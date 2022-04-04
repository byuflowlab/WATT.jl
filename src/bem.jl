struct BEM{TF, TB}
    shearexp::TF #Shear Exponent
    turbine::TB #Is a turbine (Should always be true)
    tipcorrection::TB #Include tip corrections
end

function bem(;shearexp=0.0, turbine=true, tipcorrection=true)
    return BEM(shearexp, turbine, tipcorrection)
end


function get_bem_residual(phi, rotor, section, operatingpoint) 

    ### Obtain residual 
    r, _ = CCBlade.residual(phi, rotor, section, operatingpoint) 

    if isnan(r)
        println("ccblade residual: ", r)
    end
    
    return r
end

function get_bem_residual!(outs, i, x, rotor, section, operatingpoint)
    outs[i], _ = CCBlade.residual(x[1], rotor, section, operatingpoint)
    if isnan(outs[i])
        println("ccblade residual: ", outs[i])
    end
end

function get_bem_y(phi, p, t, model::BEM, airfoil::Airfoil, env::Environment)
    # phi = x[1]
    # radius, chord, twist, pitch, rhub, rtip, hubHt = p 

    alpha = -((p[3] + p[4]) - phi) 

    Cl = airfoil.cl(alpha)
    Cd = airfoil.cd(alpha)

    Vy = env.Omega(t)*p[1]
    Vx = env.U(t)
    return [Cl, Cd, Vx, Vy]
end

function create_bemfun(model::BEM, blade::Blade, env::Environment) 
    n = length(blade.airfoils)

    function bemfun(outs, dx, x, p, t)
    
        # radius, chord, twist, pitch, rhub, rtip, hubHt = ps #Todo: I'm not sure that radius should be a differentiable variable... Like I don't think we ever really use that as a design variable... do we? 
        # Note: pitch, rhub, rtip, hubht don't change radially. 
        # Note: I'm not sure that we need rtip as a design variable.

        ### Create Rotor
        if model.tipcorrection
            rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=model.turbine)
        else
            rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=model.turbine, tip=nothing)
        end

        for i = 1:n #Iterate through the radius. 
            #Get state rates, states, and parameters for local section.
            phi = x[i]
            ps = view(p, 1+7*(i-1):7+7*(i-1)) 

            # ys = get_bem_y(xs, ps, t, model, blade.airfoils[i], env)
            Vx = env.U(0.0)
            Vy = env.Omega(0.0)*ps[1]

            U = sqrt(Vx^2 + Vy^2)

            airfoil = CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", env.rho*U*ps[2]/env.mu, U/env.a)  

            ### Create section object
            section = CCBlade.Section(ps[1], ps[2], ps[3], airfoil)

            ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
            operatingpoint = CCBlade.windturbine_op(Vx, Vy/ps[1], pitch, ps[1], 0.0, 0.0, 0.0, 0.0, ps[7], model.shearexp, env.rho, env.mu, env.a)

            #Get the residual for the local section
            get_bem_residual!(outs, i, phi, rotor, section, operatingpoint)
        end
    end
    return bemfun
end

function differentialvars(model::BEM, n)
    return fill(false, n)
end

function parsesolution(model::BEM, blade::Blade, env::Environment, p, sol)
    u = Array(sol)
    phi = u[:,end]
    n = length(phi)

    #Returnables
    cn = zeros(n)
    ct = zeros(n)
    N = zeros(n)
    T = zeros(n)

    #Working vars
    rvec = zeros(n)
    Rhub = 0.0
    Rtip = 0.0

    #Constant vars
    # B = 1
    # azimuth = 0.0
    # yaw = 0.0
    # tilt = 0.0
    # precone = 0.0 #I recognize that precone shouldn't be set to zero, but when I do the coupled solve, it should already be a part... like built into the structures... Plus... if I'm doing a steady aero, then I'd just use CCBlade. -> Althought... what happens if I want to do unsteady aero... then I'd potentially need precone.  #Todo: I should work in precone. 


    for i=1:n
        ps = p[1+7*(i-1):7+7*(i-1)]
        radius, chord, twist, pitch, rhub, rtip, hubHt = ps
        airfoil = blade.airfoils[i]

        alpha = -((twist + pitch) - phi[i])
        rvec[i] = radius 
        if i==1
            Rhub=rhub
            Rtip=rtip
        end

        cl = airfoil.cl(alpha)
        cd = airfoil.cd(alpha)

        if model.tipcorrection
            rotor = CCBlade.Rotor(rhub, rtip, 1, precone=0.0, turbine=model.turbine)
        else
            rotor = CCBlade.Rotor(rhub, rtip, 1, precone=0.0, turbine=model.turbine, tip=nothing)
        end

        Vx = env.U(0.0)
        Vy = env.Omega(0.0)*radius

        U = sqrt(Vx^2 + Vy^2) 

        alphavec = collect(-pi:0.1:pi)
        airfoil = CCBlade.AlphaAF(promote(alphavec, cl.*ones(length(alphavec)),     cd.*ones(length(alphavec)))..., "", env.rho*U*chord/env.mu, U/env.a) 

        section = CCBlade.Section(promote(radius, chord, twist)..., airfoil)

        operatingpoint = CCBlade.windturbine_op(promote(env.U(0.0), env.Omega(0.0), pitch, radius, 0.0, 0.0, 0.0, 0.0, hubHt, model.shearexp, env.rho, env.mu, env.a)...)

        _, outs = CCBlade.residual(phi[i], rotor, section, operatingpoint)
    
        N[i] = outs.Np
        T[i] = outs.Tp
        cn[i] = outs.cn
        ct[i] = outs.ct
    end

    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    rfull = [Rhub; rvec; Rtip]
    Npfull = [0.0; N; 0.0]
    Tpfull = [0.0; T; 0.0]

    # integrate Thrust and Torque (trapezoidal)
    thrust = Npfull*cos(precone)
    torque = Tpfull.*rfull*cos(precone)

    Thrust = B * FLOWMath.trapz(rfull, thrust)
    Torque = B * FLOWMath.trapz(rfull, torque)


    return phi, N, T, Thrust, Torque
end

function converge_bem_states(y, p, rotor, model::BEM, env::Environment)

    ### Extract Inputs
    Cl = y[1] 
    Cd = y[2]
    Vx = y[3]
    Vy = y[4]

    ### Extract Parameters
    radius, chord, twist, pitch, rhub, rtip, hubHt = p
    
    ### Parameters that won't change
    B = 1 # B - number of blades
    azimuth = 0.0
    yaw = 0.0
    tilt = 0.0
    precone = 0.0

    Cl, radius = promote(Cl, radius)

    omega = Vy/radius

    U = sqrt(Vx^2 + Vy^2)
    Mach = typeof(Cl)(U/env.a)
    Re = typeof(Cl)(env.rho*U*chord/env.mu)

    alphavec = collect(-pi:0.1:pi)
    airfoil = CCBlade.AlphaAF(promote(alphavec, Cl.*ones(length(alphavec)), Cd.*ones(length(alphavec)))..., "", Re, Mach)  

    ### Create section object
    section = CCBlade.Section(promote(radius, chord, twist)..., airfoil)

    ### Create OperatingPoint 
    operatingpoint = CCBlade.windturbine_op(promote(Vx, omega, pitch, radius, precone, yaw, tilt, azimuth, hubHt, model.shearexp, env.rho, env.mu, env.a)...)

    ### Converge residual 
    out = solve.(rotor, section, op)

    
    return SVector(r)
end

#Todo: Write a convenience function to go from our data structure to CCBlade structs and solve. 