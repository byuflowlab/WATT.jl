

#= 
p = radius, chord, twist, pitch, rhub, rtip, hubHt #Same as BEM
x_section = [risostates, bemcontraints]
=#

function get_bemriso_y(x_bem, x_riso, p, t, airfoil::Airfoil, env::Environment)

    ### Extract states
    phi = x_bem[1]

    ### Extract inputs
    radius, chord, twist, pitch, _, _, _ = p

    ### Find velocities in BEM frame
    Vx = env.U(t)  #Freestream velocity
    Vy = env.Omega(t)*radius 

    Vxdot = env.Udot(t) 
    Vydot = env.Omegadot(t)*radius 

    ### Find velocities in dynamic stall frame. (Not accounting for induced velocities in Dynamic stall model.)
    u = sqrt(Vy^2 + Vx^2)  
    v = 0 

    udot = sqrt(Vxdot^2 + Vydot^2)
    vdot = 0.0

    alpha = -((twist + pitch) - phi)
    alphadot = 0.0 

    yriso = [u, udot, v, vdot, alpha, alphadot]

    ### Create Airfoil object 
    Cl, Cd = riso_coefs(x_riso, yriso, chord, airfoil)


    ybem = [Cl, Cd, Vx, Vy]
    yriso = [u, udot, v, vdot, alpha, alphadot]
    return ybem, yriso
end

function create_bemrisofun(riso::Riso, bem::BEM, blade::Blade, env::Environment)
    function bemrisofun(outs, dx, x, p, t)
        n = length(blade.airfoils)
        for i=1:n
            ps = p[1+7*(i-1):7+7*(i-1)]

            bem_idx = 4*n + i
            dxs_bem = [x[bem_idx]]
            xs_bem = [x[bem_idx]]

            riso_idx = 4*(i-1)
            dxs_riso = dx[1+riso_idx:riso_idx+4]
            xs_riso = x[1+riso_idx:riso_idx+4]
            ps_riso = ps[2]
            
            ys_bem, ys_riso = get_bemriso_y(xs_bem, xs_riso, ps, t, blade.airfoils[i], env)

            outs[1+riso_idx:riso_idx+4] = riso_residual(dxs_riso, xs_riso, ys_riso, ps_riso, t, blade.airfoils[i])

            outs[bem_idx] = get_bem_residual(dxs_bem, xs_bem, ys_bem, ps, t, bem, env)[1]

        end
    end
    return bemrisofun
end

function differentialvars(bemmodel::BEM, risomodel::Riso, n)
    return vcat(differentialvars(risomodel, n), differentialvars(bemmodel,n))
end

function parsesolution(bem::BEM, riso::Riso, blade::Blade, env::Environment, p, sol)
    u = Array(sol)'
    t = sol.t

    m = length(t)
    n = length(blade.airfoils)

    ubem = u[:, end-n+1:end]
    uds = u[:,1:end-n]

    #Returnables
    cn = zeros(m,n)
    ct = zeros(m,n)
    N = zeros(m,n)
    T = zeros(m,n)
    Torque = zeros(m)
    Thrust = zeros(m)

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

    for i = 1:m #Iterate through time
        for j=1:n
            phi = ubem[i,:]
            ps = p[1+7*(j-1):7+7*(j-1)]
            radius, chord, twist, pitch, rhub, rtip, hubHt = ps
            airfoil = blade.airfoils[j]

            alpha = -((twist + pitch) - phi[j])
            
            if i==1
                rvec[j] = radius 
                Rhub=rhub
                Rtip=rtip
            end

            cl = airfoil.cl(alpha)
            cd = airfoil.cd(alpha)

            if bem.tipcorrection
                rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bem.turbine)
            else
                rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bem.turbine, tip=nothing)
            end

            Vx = env.U(t[i])
            Vy = env.Omega(t[i])*radius

            U = sqrt(Vx^2 + Vy^2)
            Mach = typeof(cl)(U/env.a)
            Re = typeof(cl)(env.rho*U*chord/env.mu)

            alphavec = collect(-pi:0.1:pi)
            airfoil = CCBlade.AlphaAF(promote(alphavec, cl.*ones(length(alphavec)),     cd.*ones(length(alphavec)))..., "", Re, Mach) 

            section = CCBlade.Section(promote(radius, chord, twist)..., airfoil)

            operatingpoint = CCBlade.windturbine_op(promote(env.U(t[i]), env.Omega(t[i]), pitch, radius, precone, yaw, tilt, azimuth, hubHt, bem.shearexp, env.rho, env.mu, env.a)...)

            _, outs = CCBlade.residual(phi[j], rotor, section, operatingpoint)
            
            N[i,j] = outs.Np
            T[i,j] = outs.Tp
            cn[i,j] = outs.cn
            ct[i,j] = outs.ct
        end
    end

    rfull = [Rhub; rvec; Rtip]

    for k=1:m
        Npfull = [0.0; N[k,:]; 0.0]
        Tpfull = [0.0; T[k,:]; 0.0]
        
        # integrate Thrust and Torque (trapezoidal)
        thrusti = Npfull.*cos(precone)
        torquei = Tpfull.*rfull.*cos(precone)

        Thrust[k] = B * FLOWMath.trapz(rfull, thrusti)
        Torque[k] = B * FLOWMath.trapz(rfull, torquei)
    end
    return t, ubem, uds, N, T, Thrust, Torque
end