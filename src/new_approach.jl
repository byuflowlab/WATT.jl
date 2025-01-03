#=
Re-writing Rotors, but doing everything out-of-place in the hopes that it'll be faster with Reverse-mode and so I can dodge the memory issues I'm having. -> I think it's a scoping issue. 

12/31/24 Adam Cardoza
=#

import CCBlade.afeval

function find_inittype(vars...) #todo: How does this vary from promote_type? -> That's probably a more robust approach. 
    # println("Initialization types")
    for item in vars
        # println(typeof(item))
        if isa(item, ForwardDiff.Dual)
            return typeof(item)
        elseif isa(item, ReverseDiff.TrackedReal)
            return typeof(item)
        end
    end

    return eltype(vars)
end

function afeval(af::DS.Airfoil, alpha, Re, Mach)
    return af.cl(alpha), af.cd(alpha)
end

function solve_BEM(rotor::Rotor, blade::Blade, env::Environment, idx, Vx, Vy, pitch; twist=blade.twist[idx], npts::Int=10, epsilon=1e-6) #epsilon=1.1e-3)
    
    
    
    airfoil = blade.airfoils[idx]
    # rR = blade.rR[idx]
    r = blade.r[idx]
    

    # check if we are at hub/tip
    # if isapprox(rR, 0.0, atol=1e-6) || isapprox(rR, 1.0, atol=1e-6)
    #     return Outputs()  # no loads at hub/tip
    # end

    ### unpack
    theta = twist + pitch



    # ---- determine quadrants based on case -----
    Vx_is_zero = isapprox(Vx, 0.0, atol=1e-6)
    Vy_is_zero = isapprox(Vy, 0.0, atol=1e-6)

    phi_geo = atan(Vx, Vy)
    if abs(phi_geo)<=epsilon
        warn("The geometric inflow angle is approximately the same as the solution ϵ, which may cause the BEMT to fail.")
    end

    ### quadrants
    # epsilon = 1e-6
    q1 = [epsilon, pi/2]
    q2 = [-pi/2, -epsilon]
    q3 = [pi/2, pi-epsilon]
    q4 = [-pi+epsilon, -pi/2]

    if Vx_is_zero && Vy_is_zero
        println("Vx and Vy is zero.")
        return CCBlade.Outputs()

    elseif Vx_is_zero

        startfrom90 = false  # start bracket at 0 deg.

        if Vy > 0 && theta > 0
            order = (q1, q2)
        elseif Vy > 0 && theta < 0
            order = (q2, q1)
        elseif Vy < 0 && theta > 0
            order = (q3, q4)
        else  # Vy < 0 && theta < 0
            order = (q4, q3)
        end

    elseif Vy_is_zero

        startfrom90 = true  # start bracket search from 90 deg

        if Vx > 0 && abs(theta) < pi/2
            order = (q1, q3)
        elseif Vx < 0 && abs(theta) < pi/2
            order = (q2, q4)
        elseif Vx > 0 && abs(theta) > pi/2
            order = (q3, q1)
        else  # Vx < 0 && abs(theta) > pi/2
            order = (q4, q2)
        end

    else  # normal case

        startfrom90 = false

        if Vx > 0 && Vy > 0
            order = (q1, q2, q3, q4)
        elseif Vx < 0 && Vy > 0
            order = (q2, q1, q4, q3)
        elseif Vx > 0 && Vy < 0
            order = (q3, q4, q1, q2)
        else  # Vx[i] < 0 && Vy[i] < 0
            order = (q4, q3, q2, q1)
        end

    end

    # ----- solve residual function ------
    # pull out first argument
    residual(phi, x, p) = CCBlade.residual_and_outputs(phi, x, p)[1]

    # package up variables and parameters for residual 
    xv = (r, airfoil.c, twist, blade.rhub, blade.rtip, Vx, Vy, env.rho, pitch, env.mu, env.a)
    pv = (airfoil, rotor.B, rotor.turbine, rotor.re, rotor.mach, rotor.rotation, rotor.tip)

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

        # find bracket
        # @show xv
        # if isa(xv[1], ReverseDiff.TrackedReal)
        #     @show airfoil.c.value
        #     println([xv[i].value for i in eachindex(xv)])
        # end
        success, phiL, phiU = CCBlade.firstbracket(phi -> residual(phi, xv, pv), phimin, phimax, npts, backwardsearch)

        # once bracket is found, solve root finding problem and compute loads
        if success
            function solve(x, p) #todo: Is there a more efficient way to do this instead of a closure? 
                phistar, _ = FLOWMath.brent(phi -> residual(phi, x, p), phiL, phiU)
                
                return phistar
            end
            
            phistar = IAD.implicit(solve, residual, xv, pv)
            
            _, outputs = CCBlade.residual_and_outputs(phistar, xv, pv) #TODO: Instead of creating a new function... I could just check if x0 is inbetween a and b outside of the loop. If it is, then I can check if it is a zero. If not, then replace one of the bounds into Brent's method. 
                return outputs
            end    
        end    

    # it shouldn't get to this point.  if it does it means no solution was found
    # it will return empty outputs
    # alternatively, one could increase npts and try again
    
    @warn "Invalid data (likely) for this section.  Zero loading assumed."
    return CCBlade.Outputs()
end

function get_airfoil_initial_condition(airfoil, tvec, twist, W, phi, pitch, turbine)

    theta = (twist + pitch) - phi #TODO: Isn't this alpha??? 
    
    if turbine
        theta *= -1
    end

    ys = [W, 0.0, 0.0, theta]
    
    xds0, _ = DS.initialize(airfoil, tvec, ys)
    return xds0
end

function dsmodel_initial_condition(phi, W, blade::Blade, turbine::Bool, tvec, pitch)

    na = length(blade.r)
    xds0 = [get_airfoil_initial_condition(blade.airfoils[i], tvec, blade.twist[i], W[i], phi[i], pitch, turbine) for i in 1:na]

    return xds0
end

function extract_ds_loads(airfoil::Airfoil, states, twist, phi, pitch)

    theta = (twist + pitch) - phi #TODO: Isn't this alpha??? 
    
    if turbine
        theta *= -1
    end

    ys = [W, 0.0, 0.0, theta]

    Cl, Cd, Cm = DS.get_loads(airfoil.model, airfoil, states, ys)

    sphi, cphi = sincos(phi)
    Cx = Cl*cphi + Cd*sphi
    Cy = -(Cl*sphi - Cd*cphi)

    return Cx, Cy, Cm
end

function dimensionalize(Cx, Cy, Cm, W, airfoil::Airfoil, env::Environment)
    
    q_local = 0.5*env.rho*W^2 #Local dynamic pressure
    
    Fx = Cx*q_local*airfoil.c 
    Fy = Cy*q_local*airfoil.c 
    Mx = Cm*q_local*airfoil.c^2 #The coefficient of moment is positive about the negative Z aero axis, so we need the negative of this to move it to the structural X axis. 
    return Fx, Fy, Mx
end

function get_aero_loads(airfoils, states, phi, W, pitch)

    twistvec = 1 #Todo:

    Cx, Cy, Cm = extract_ds_loads.(airfoils, states, twist, phi, pitch)

    Fx, Fy, Mx = dimensionalize.(Cx, Cy, Cm, W, airfoils, Ref(env))

    return Fx, Fy, Mx
end

function get_distributed_loads(assembly, fx, fy, mx)

    nelem = length(assembly.elements)
    
    distributed_loads = Dict(i => GXBeam.DistributedLoads(assembly_i, i; fy = (s) -> fy[i], fz = (s) -> -fx[i]) for i in 1:nelem)

    return distributed_loads
end

function initial_condition(rotor::Rotors.Rotor, blade, mesh, env::Environment, tvec; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, g=9.81, runtimeflag::Bool=false, runtimeiter::Int=speakiter, runtime = (aerostates, gxhistory, i) ->nothing, gxflag=nothing, prepp=nothing, p=nothing, azimuth0=0.0)

    # println("unpacking...")
    ### unpack the data structures. 
    @unpack assembly, system, prescribed_conditions, distributed_loads, point_masses, linear_velocity, xpfunc, pfunc, structural_damping, two_dimensional, linear = mesh

    @unpack azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds = aerostates #Todo: I want these states to be created, rather than done in-place. 

    na = length(blade.r)
    nt = length(tvec)

    ### Initial Condition analysis. 
    t0 = tvec[1] 


    if verbose
        println("Calculating initial condition...")
    end

    initial_condition_checks(gxflag)

    # println("Initializing...")
    ### Initialize BEM solution 
    Vx, Vy = get_aero_velocities.(Ref(rotor), Ref(blade), Ref(env), Ref(t0), 1:na, Ref(azimuth0))
    ccout = solve_BEM(Ref(rotor), Ref(blade), Ref(env), 1:na, Vx, Vy, Ref(pitch))

    phi0 = ccout.phi #todo: I don't know if I need to extract it. 
    alpha0 = ccout.alpha
    W0 = ccout.W

    

    #Note: I don't use take_aero_step because these functions are different. 
    xds0 = dsmodel_initial_condition(phi0, W, blade, turbine, tvec, pitch)

    
    # extract_ds_loads!(blade.airfoils, xds0, mesh.xds_idxs, phi0, mesh.p_ds, cx0, cy0, cm0)


    # dimensionalize!(fx0, fy0, mx0, cx0, cy0, cm0, blade::Blade, env::Environment, W0) 

    
    update_forces!(distributed_loads, fx0, fy0, mx0, blade, assembly) #Todo: Doesn't matter, because Taylor's adjoint dooesn't use distributed loads, it recreates it. -> Will it matter if this is run without the pfunc? -> I think it will. -> I'll put this in the following if statement. 


    ### GXBeam initial solution 
    if isnothing(prepp)
        p = nothing
        distributed_loads
    else
        prepp(p, fx0, fy0, mx0)
    end

    Omega0 = SVector(0.0, 0.0, -env.RS(t0))
    gravity0 = SVector(-g*cos(azimuth0), -g*sin(azimuth0), 0.0)

    # #todo: Only works with initial response (not steady state or spinning solution. )
    # system, gxhistory[1], converged = GXBeam.initial_condition_analysis!(system, assembly, t0; prescribed_conditions, distributed_loads, angular_velocity=Omega0, gravity=gravity0, steady_state=false, structural_damping, linear, pfunc, p, show_trace=false)

    # @show typeof(prescribed_conditions)
    # @show typeof(distributed_loads)
    # @show typeof(assembly)
    # @show typeof(gravity0), typeof(Omega0)
    # @show typeof(system)
    # @show typeof(p)
    # @show p

    system, gxstate, constants, paug, xgx, converged = GXBeam.initialize_system!(system, assembly, tvec; prescribed_conditions, distributed_loads, gravity=gravity0, angular_velocity=Omega0, structural_damping, reset_state=true, pfunc, p) #todo: This has extra allocations that I don't need.
    
    # @show typeof(gxstate), length(gxstate)
    # println("Finished initializing GXBeam...")
    # @show typeof(paug)

    #Todo. What does gxstate look like? -> A vector of assembly states (but that's because it's reallocating the history vector. )

    gxhistory[1] = gxstate[1]

    if !converged
        @warn("The initial condition structural analysis failed to converge.")
        return 
    end



    ### Update mesh transfer variables
    update_mesh!(blade, mesh, assembly, gxhistory[1], env, t0, na)


    
    ### Take the first step then set that as the first value. -> This is what OpenFAST does. 
    # azimuth = env.RS(t)*dt + azimuth0
    
    if isa(xds0[1], ReverseDiff.TrackedReal)
        # ns = DS.numberofstates_total(airfoils)
        # xds_old = Array{inittype, 1}(undef, ns) 
        # xds_old = deepcopy(xds)
        inittype = eltype(xds0)
        ns = DS.numberofstates_total(blade.airfoils)
        xds_old = Array{inittype, 1}(undef, ns)

        #Todo. How do I copy the derivative? (Because I want to preserve the pointers to xds.)
        for i in eachindex(xds0) 
            xds_old[i] = xds0[i]
        end
        
    elseif isa(xds[1], ForwardDiff.Dual)
        # xds_old = deepcopy(xds) #Derivatives are getting nuked
        inittype = eltype(xds0)
        ns = DS.numberofstates_total(blade.airfoils)
        xds_old = Array{inittype, 1}(undef, ns)

        #Todo. How do I copy the derivative? (Because I want to preserve the pointers to xds.)
        for i in eachindex(xds0) 
            xds_old[i] = xds0[i]
        end

        # @show xds[1].partials
        # @show xds_old[1].partials
    else
        xds_old = deepcopy(xds0)
    end

    dt = tvec[2] - tvec[1] #Note: this passing in phi0 and replacing phi0 might break things... 
    # take_aero_step!(phi0, alpha0, W0, xds0, cx0, cy0, cm0, fx0, fy0, mx0, phi0, xds_old, azimuth0, t0, dt, pitch, mesh, rotor, blade, env; solver, pfunc, prepp, p)
    take_aero_step!(phi0, alpha0, W0, xds0, cx0, cy0, cm0, fx0, fy0, mx0, xds_old, azimuth0, t0, dt, pitch, mesh, rotor, blade, env; solver)
end



function take_step(rotor::Rotors.Rotor, blade, mesh, env::Environment, tvec; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, g=9.81, runtimeflag::Bool=false, runtimeiter::Int=speakiter, runtime = (aerostates, gxhistory, i) ->nothing, gxflag=nothing, prepp=nothing, p=nothing, azimuth0=0.0)




    # println("unpacking...")
    ### unpack the data structures. 
    @unpack assembly, system, prescribed_conditions, distributed_loads, point_masses, linear_velocity, xpfunc, pfunc, structural_damping, two_dimensional, linear = mesh

    @unpack azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds = aerostates #Todo: I want these states to be created, rather than done in-place. 

    na = length(blade.r)
    nt = length(tvec)


    # println("Beginning time loop...")
    for i in 2:nt

        ### Unpack
        phi_i = view(phi, i, :)
        alpha_i = view(alpha, i, :)
        W_i = view(W, i, :)
        cx_i = view(Cx, i, :)
        cy_i = view(Cy, i, :)
        cm_i = view(Cm, i, :)
        fx_i = view(Fx, i, :)
        fy_i = view(Fy, i, :)
        mx_i = view(Mx, i, :)
        xds_i = view(xds, i, :)
        # phi_im1 = view(phi, i-1, :)
        xds_im1 = view(xds, i-1, :)
        
        t = tvec[i]
        tprev = tvec[i-1]
        dt = t - tprev
    
        
        if dt<0
            error("Time step is negative")
        end

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. 
        #todo: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

        if azimuth[i]<azimuth[i-1]
            @warn("Blade moved backwards")
        end

        #Note: Coupling from structural velocities doesn't appear to kick in until the third time step. 
        # take_aero_step!(phi_i, alpha_i, W_i, xds_i, cx_i, cy_i, cm_i, fx_i, fy_i, mx_i, phi_im1, xds_im1, azimuth[i], t, dt, pitch, mesh, rotor, blade, env; solver, pfunc, prepp, p)
        take_aero_step!(phi_i, alpha_i, W_i, xds_i, cx_i, cy_i, cm_i, fx_i, fy_i, mx_i, xds_im1, azimuth[i], t, dt, pitch, mesh, rotor, blade, env; solver)

        # @show fx_i[end]
        # println("")


        ### Update GXBeam loads 
        # @show distributed_loads
        update_forces!(distributed_loads, fx_i, fy_i, mx_i, blade, assembly) 
        # println("") #Note: Changing, but going super negative? 
        # @show distributed_loads

        Omega = SVector(0.0, 0.0, -env.RS(t)) 
        # gravity = SVector(-g*cos(aerostates.azimuth[i-1]), -g*sin(aerostates.azimuth[i-1]), 0.0) #TODO: I need to include tilt, and precone here. 

        # @show typeof(aerostates.azimuth)

        if isa(azimuth[i], ForwardDiff.Dual)||isa(azimuth[i], ReverseDiff.TrackedReal) #todo: Does azimuth need to be a state?
            a0 = azimuth[i].value
            a1 = azimuth[i-1].value
        else
            a0 = azimuth[i]
            a1 = azimuth[i-1]
        end

        # if isnothing(p)
        #     a0 = azimuth
        #     a1 = azimuth0
        # else
        #     a0 = azimuth.value
        #     a1 = azimuth0.value
        # end

        #TODO: Figure out if I need to use this, or if I can just calculate it at a single point. 
        gravity = (tee) -> SVector(-g*cos((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), -g*sin((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), 0.0) ##Todo t = tvec[i].... So this be way wrong. Oh... this is an inline function so the solver can linearly interpolate the gravity vector across time. But... I think the time domain analysis only analyzes at the given time steps... which means that this function doesn't get called really... I dunno. -> TODO: This would only make a difference if the solver used intermediate steps (using DifferentialEquations to solve GXBeam.) -> I should get rid of this. 

        # @show typeof(gravity2)

        #Note: Taylor applies the gravitational load by C'*mass*C*gvec

        ### Solve GXBeam for time step 
        #todo: This function is taking a lot of time. -> I might be able to save time by branching his code and writing another function, but most of the time is spent in nlsolve. I think all of the time spent is just time solving, not really inside of Taylor's code, but of course, if I make his code faster, then I make the solve faster. 
        # @show eltype(system)

        if isnothing(prepp)  
            p = nothing
            # println("Got here")
        else
            prepp(p, fx_i, fy_i, mx_i)
        end #Todo: derivatives passing needs to be included. 

    

        system, gxhistory[i], constants, paug, xgx, convergedi = GXBeam.step_system!(system, paug, xgx, constants, gxhistory[i-1], assembly, tvec, i; prescribed_conditions, distributed_loads, structural_damping, gravity, angular_velocity=Omega)

        if !convergedi
            @warn("GXBeam failed to converge on the $i th time step.")
            break
        end




        ### Update aero inputs from structures.
        update_mesh!(blade, mesh, assembly, gxhistory[i], env, t, na)



        if verbose & (mod(i-1, speakiter)==0) #todo: remove the dependence on i (move to just a verbose and runtime flag)
            println("")
            println("Simulation time: ", t)
        end



        if runtimeflag & (mod(i-1, runtimeiter)==0) 
            runtime(aerostates, gxhistory_new, i) 
        end
    end
end

function simulate()
end