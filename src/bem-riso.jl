



#= 
p = radius, chord, twist, pitch, rhub, rtip, hubHt #Same as BEM
x_section = [risostates, bemcontraints]
=#

function get_bemriso_y(phi, x_riso, p, t, airfoil::Airfoil, env::Environment)

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

    yriso = SVector(u, udot, v, vdot, alpha, alphadot)

    ### Create Airfoil object 
    Cl, Cd = riso_coefs(x_riso, yriso, chord, airfoil)

    ybem = SVector(Cl, Cd, Vx, Vy)
    # yriso = [u, udot, v, vdot, alpha, alphadot]
    return ybem, yriso
end

function create_bemrisofun_wstates(riso::Riso, bem::BEM, blade::Blade, env::Environment)
    n = length(blade.airfoils)

    rarray = Array{CCBlade.Rotor}(undef, 1)
    afarray = Array{CCBlade.AFType}(undef,1)
    secarray = Array{CCBlade.Section}(undef,1)
    oparray = Array{CCBlade.OperatingPoint}(undef, 1)

    alphavec = collect(-pi:.1:pi)
    

    function bemrisofun(outs, dx, x, p, t)
        clvec = ones(eltype(outs), length(alphavec))
        cdvec = ones(eltype(outs), length(alphavec))

        ### Create Rotor
        if bem.tipcorrection #Question: I wonder if there is a way that I can just create these.... then update them. 
            # rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine)
            rarray[1] = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine)
        else
            # rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
            rarray[1] = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
        end

        for i=1:n
            ps = view(p, 1+7*(i-1):7+7*(i-1)) #p[1+7*(i-1):7+7*(i-1)]

            bem_idx = 4*n + i
            xbem = view(x, bem_idx)

            riso_idx = 4*(i-1)
            dxs_riso = view(dx, 1+riso_idx:riso_idx+4) #dx[1+riso_idx:riso_idx+4]
            xs_riso = view(x, 1+riso_idx:riso_idx+4)  #x[1+riso_idx:riso_idx+4]
            ps_riso = ps[2]
            
            ys_bem, ys_riso = get_bemriso_y(xbem[1], xs_riso, ps, t, blade.airfoils[i], env)

            outs[1+riso_idx:riso_idx+4] = riso_residual(dxs_riso, xs_riso, ys_riso, ps_riso, t, blade.airfoils[i])

            U = sqrt(ys_bem[3]^2 + ys_bem[4]^2)
            clvec[:] .= ys_bem[1]
            cdvec[:] .= ys_bem[2]
             
            afarray[1] = CCBlade.AlphaAF(alphavec, clvec, cdvec, "", env.rho*U*ps[2]/env.mu, U/env.a) #Note: We definitely want to leave the airfoils in the replacing array, they take a crap ton of space. 
            
            ### Create section object
            section = CCBlade.Section(ps[1], ps[2], ps[3], afarray[1])
            # secarray[1] = CCBlade.Section(ps[1], ps[2], ps[3], afarray[1])

            ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
            # operatingpoint = CCBlade.windturbine_op(ys_bem[3], ys_bem[4]/ps[1], pitch, ps[1], 0.0, 0.0, 0.0, 0.0, ps[7], bem.shearexp, env.rho, env.mu, env.a)
            oparray[1] = CCBlade.windturbine_op(ys_bem[3], ys_bem[4]/ps[1], pitch, ps[1], 0.0, 0.0, 0.0, 0.0, ps[7], bem.shearexp, env.rho, env.mu, env.a)

            # get_bem_residual!(outs, bem_idx, [x[bem_idx]], rotor, section, operatingpoint)
            get_bem_residual!(outs, bem_idx, xbem[1], rarray[1], section, oparray[1])
        end
    end
    return bemrisofun
end

function differentialvars(bemmodel::BEM, risomodel::Riso, n;withstates=false)
    if withstates
        return vcat(differentialvars(risomodel, n), differentialvars(bemmodel,n))
    else
        return differentialvars(risomodel, n)
    end
end


# function get_bemriso_residual(phi, xs_riso, ps, t, bem::BEM, airfoil::Airfoil, env::Environment)

#     ys_bem, _ = get_bemriso_y(phi, xs_riso, ps, t, airfoil, env)

#     return get_bem_residual(phi, ys_bem, ps, t, bem, env)
# end


function converge_bemriso_residual(Vx, Vy, ps, rotor, section, operatingpoint)

    radius, chord, twist, pitch, rhub, rtip, hubHt = ps

    # quadrants
    epsilon = 1e-6
    q1 = [epsilon, pi/2]
    q2 = [-pi/2, -epsilon]
    q3 = [pi/2, pi-epsilon]
    q4 = [-pi+epsilon, -pi/2]

    startfrom90 = false
    npts = 20

    # Vx = env.U(t)
    # Vy = env.Omega(t)*radius

    if Vx > 0 && Vy > 0
        # println("option 1")
        order = (q1, q2, q3, q4)
    elseif Vx < 0 && Vy > 0
        # println("option 2")
        order = (q2, q1, q4, q3)
    elseif Vx > 0 && Vy < 0
        # println("option 3")
        order = (q3, q4, q1, q2)
    else  # Vx[i] < 0 && Vy[i] < 0
        # println("Option 4")
        order = (q4, q3, q2, q1)
    end

    R(phistar) = get_bem_residual(phistar, rotor, section, operatingpoint)

    success = false
    for j = 1:length(order)  # quadrant orders.  In most cases it should find root in first quadrant searched.
        phimin, phimax = order[j]

        # check to see if it would be faster to reverse the bracket search direction
        backwardsearch = false
        if !startfrom90
            if phimin == -pi/2 || phimax == -pi/2  # q2 or q4
                backwardsearch = true
            end
        else
            if phimax == pi/2  # q1
                backwardsearch = true
            end
        end
        
        # force to dual numbers if necessary
        phimin = phimin*one(chord)
        phimax = phimax*one(chord)

        # find bracket
        success, phiL, phiU = CCBlade.firstbracket(R, phimin, phimax, npts, backwardsearch)

        # once bracket is found, solve root finding problem and compute loads
        if success
            phistar, _ = FLOWMath.brent(R, phiL, phiU)
            println("Quadrant: ", j, "   phi: ", phistar)
            return phistar #, j
        end    
    end
    @warn "Invalid data likely for this section. Radius=$radius"
end


function create_bemrisofun(riso::Riso, bem::BEM, blade::Blade, env::Environment) 
    n = length(blade.airfoils)

    function aerofun(outs, dx, x, p, t)
        #p = radius, chord, twist, pitch, rhub, rtip, hubHt

        ### Create Rotor
        if bem.tipcorrection #Question: I wonder if there is a way that I can just create these.... then update them. 
            # rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine)
            rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine)
        else
            # rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
            rotor = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
        end
        for i=1:n
            ps = view(p, 1+7*(i-1):7+7*(i-1))
            radius, chord, twist, pitch, rhub, rtip, hubHt = ps

            riso_idx = 4*(i-1)
            dxs_riso = view(dx, 1+riso_idx:riso_idx+4)
            xs_riso = view(x, 1+riso_idx:riso_idx+4)
            ps_riso = ps[2]

            Vx = env.U(t)
            Vy = env.Omega(t)*ps[1]
            Vxdot = env.Udot(t)
            Vydot = env.Omegadot(t)*ps[1]
            alphadot = 0.0

            clfun, cdfun = riso_coefs_funs(xs_riso, Vx, chord, alphadot, blade.airfoils[i])
            #Todo: I don't think that this includes twist or pitch. And I'm not sure that it goes between phi and alpha correctly. -> I'm not sure that it needs to. Right? It's just Cl and Cd as a function of alpha, then CCBlade should take care of the rest. ... 

            airfoil = CCBlade.AFfun(clfun, cdfun)

            section = CCBlade.Section(radius, chord, twist, airfoil)

            operatingpoint = CCBlade.windturbine_op(Vx, Vy, pitch, radius, 0.0, 0.0, 0.0, 0.0, hubHt, bem.shearexp, env.rho, env.mu, env.a)

            phi = converge_bemriso_residual(Vx, Vy, ps, rotor, section, operatingpoint)
            # out = CCBlade.solve(rotor, section, operatingpoint)
            
            # _, ys_riso = get_bemriso_y(phi, xs_riso, ps, t, blade.airfoils[i], env)

            u = sqrt(Vy^2 + Vx^2)  
            v = 0 

            udot = sqrt(Vxdot^2 + Vydot^2)
            vdot = 0.0

            alpha = -((twist + pitch) - phi)
            alphadot = 0.0
            ys_riso = SVector(u, udot, v, vdot, alpha, alphadot)

            outs[1+riso_idx:riso_idx+4] = riso_residual(dxs_riso, xs_riso, ys_riso, ps_riso, t, blade.airfoils[i])
            if i==n
                println(t)
                # println(phi)
            end
        end
    end
    return aerofun
end

function create_bemrisoODE(riso::Riso, bem::BEM, blade::Blade, env::Environment)
    function aerofun(x, p, t)
        n = length(blade.airfoils)
        dx = zeros(4*n)
        # println("")
        # println("t = ", t)
        for i=1:n
            
            ps = p[1+7*(i-1):7+7*(i-1)]

            riso_idx = 4*(i-1)
            xs_riso = x[1+riso_idx:riso_idx+4]
            ps_riso = ps[2]
            # println("x$i: ", xs_riso)

            phi = converge_bemriso_residual(xs_riso, ps, t, bem, blade.airfoils[i], env)
            
            _, ys_riso = get_bemriso_y([phi], xs_riso, ps, t, blade.airfoils[i], env)
            u, udot, v, vdot, theta, thetadot = ys_riso

            dx[1+riso_idx:4+riso_idx] = riso_states(xs_riso, u, udot, v, vdot, theta, thetadot, ps_riso, blade.airfoils[i])
            # if i==n
            #     println(j, "    ", phi) #Phi is flying all over the place. ALL over the place. 
            # end
        end
        return dx
    end
    return aerofun
end


function parsesolution(bem::BEM, riso::Riso, blade::Blade, env::Environment, p, sol)
    n = length(blade.airfoils)
    u = Array(sol)'
    _, m = size(u)
    i = m/n
    if i>4
        return parsesolution_wstates(bem, blade, env, p, sol)
    else
        return parsesolution_nostates(bem, blade, env, p, sol)
    end
end

function parsesolution_wstates(bem::BEM, blade::Blade, env::Environment, p, sol)
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

function parsesolution_nostates(bem::BEM, blade::Blade, env::Environment, p, sol)
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