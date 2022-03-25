struct BEM{TF, TB}
    shearexp::TF #Shear Exponent
    turbine::TB #Is a turbine (Should always be true)
    tipcorrection::TB #Include tip corrections
end

function bem(;shearexp=0.0, turbine=true, tipcorrection=true)
    return BEM(shearexp, turbine, tipcorrection)
end


function get_bem_residual(dx, x, y, p, t, model::BEM, env::Environment) #TODO: Get rid of inputs that aren't used.... I'm not sure that my formulation needs them. 
    ### Extract statess
    phi = x[1] 
    # println("phi: ", phi)

    ### Extract Inputs
    Cl = y[1] 
    Cd = y[2]
    Vx = y[3]
    Vy = y[4]

    # println("phi: ", phi)
    # println("ccblade y: ", y)

    ### Extract Parameters
    radius, chord, twist, pitch, rhub, rtip, hubHt = p #Todo: I'm not sure that radius should be a differentiable variable... Like I don't think we ever really use that as a design variable... do we? 

    ### Parameters that won't change
    B = 1 # B - number of blades
    azimuth = 0.0
    yaw = 0.0
    tilt = 0.0
    precone = 0.0

    Cl, radius, phi = promote(Cl, radius, phi)

    omega = Vy/radius

    if model.tipcorrection
        rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=model.turbine)
    else
        rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=model.turbine, tip=nothing)
    end

    U = sqrt(Vx^2 + Vy^2)
    Mach = typeof(Cl)(U/env.a)
    Re = typeof(Cl)(env.rho*U*chord/env.mu)
    # println(typeof(Cl))

    alphavec = collect(-pi:0.1:pi)
    airfoil = CCBlade.AlphaAF(promote(alphavec, Cl.*ones(length(alphavec)), Cd.*ones(length(alphavec)))..., "", Re, Mach)  

    ### Create section object
    section = CCBlade.Section(promote(radius, chord, twist)..., airfoil)

    ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
    # println(Vx, omega, pitch, radius, precone, yaw, tilt, azimuth, hubHt, model.shearExp, model.rho, model.mu, model.a) #I'm getting Cl and Cd as nan
    operatingpoint = CCBlade.windturbine_op(promote(Vx, omega, pitch, radius, precone, yaw, tilt, azimuth, hubHt, model.shearexp, env.rho, env.mu, env.a)...)

    ### Obtain residual 
    r, outsvec = CCBlade.residual(phi, rotor, section, operatingpoint) #TODO: This function is not compatible with dual numbers. - I did a thing in CCBlade that might mean I won't need to do promote so many times here... I don't know if it's the best... Maybe I can define a second method that won't be as fast, but won't modify the original.  - Did that, I wonder what Dr. Ning's response to that will be. -> His response was... why would you ever need to change that? I still have typing that needs to be updated in CCBlade.

    if isnan(r)
        println("ccblade residual: ", r)
    end
    
    return SVector(r)
end

function get_bem_y(dx, x, p, t, model::BEM, airfoil::Airfoil, env::Environment)
    phi = x[1]

    radius, chord, twist, pitch, rhub, rtip, hubHt = p 

    alpha = -((twist + pitch) - phi) 

    Cl = airfoil.cl(alpha)
    Cd = airfoil.cd(alpha)

    Vy = env.Omega*radius
    Vx = env.Uinf
    return [Cl, Cd, Vx, Vy]
end

function create_bemfun(model::BEM, blade::Blade, env::Environment) 
    function bemfun(outs, dx, x, p, t)
        n = length(blade.airfoils)

        for i = 1:n
            #Get state rates, states, and parameters for local section.
            dxs = [dx[i]]
            xs = [x[i]]
            ps = p[1+7*(i-1):7+7*(i-1)]
            ys = get_bem_y(dxs, xs, ps, t, model, blade.airfoils[i], env)

            #Get the residual for the local section
            outs[i] = get_bem_residual(dxs, xs, ys, ps, t, model, env)[1]
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
    B = 1
    azimuth = 0.0
    yaw = 0.0
    tilt = 0.0
    precone = 0.0 #I recognize that precone shouldn't be set to zero, but when I do the coupled solve, it should already be a part... like built into the structures... Plus... if I'm doing a steady aero, then I'd just use CCBlade. -> Althought... what happens if I want to do unsteady aero... then I'd potentially need precone.  #Todo: I should work in precone. 


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
            rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=model.turbine)
        else
            rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=model.turbine, tip=nothing)
        end

        Vx = env.Uinf
        Vy = env.Omega*radius

        U = sqrt(Vx^2 + Vy^2)
        Mach = typeof(cl)(U/env.a)
        Re = typeof(cl)(env.rho*U*chord/env.mu)

        alphavec = collect(-pi:0.1:pi)
        airfoil = CCBlade.AlphaAF(promote(alphavec, cl.*ones(length(alphavec)),     cd.*ones(length(alphavec)))..., "", Re, Mach) 

        section = CCBlade.Section(promote(radius, chord, twist)..., airfoil)

        operatingpoint = CCBlade.windturbine_op(promote(env.Uinf, env.Omega, pitch, radius, precone, yaw, tilt, azimuth, hubHt, model.shearexp, env.rho, env.mu, env.a)...)

        r_, outs = CCBlade.residual(phi[i], rotor, section, operatingpoint)
    
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