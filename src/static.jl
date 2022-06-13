



function fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations=1000, tolerance=1e-3, g=0.0, verbose=false)


    ### Seperate the parameters and read in any pertainent data
    n = gxmodel.n
    pa = view(p, 1:(7*n))
    ps = view(p, (7*n+1):length(p))

    ### Constant problem dependent parameters
    pitch = pa[4]
    rhub = pa[5]
    rtip = pa[6]
    hubHt = pa[7]

    ### Environmental variables
    Omega = SVector(0.0, 0.0, -env.Omega(0.0))
    grav = SVector(0.0, -g, 0.0)
    U = env.U(0.0)


    ### Base vectors for turbine. 
    radii_idx = 1:7:(7*n)
    chord_idx = 2:7:(7*n)
    twist_idx = 3:7:(7*n)
    rvec = view(pa, radii_idx) #Todo: This isn't what I want. Currently the BEM node and the GXBeam node aren't co-located. 
    chordvec = view(pa, chord_idx)
    twistvec = view(pa, twist_idx)

    ### Unchanging constants
    t = 0.0
    precone = 0.0 #Precone. -> I'm going to change the radius manually. -> I could have a base precone, then allow for changes on top of that? 
    B = 1 #Number of blades
    yaw = 0.0
    tilt = 0.0
    azimuth = 0.0





    ### Create storage structures for solution
    resids = zeros(maxiterations, 3) #columns of resids is deflection y, force y, force z





    #### Create CCBlade inputs
    ### Create Rotor
    if bemmodel.tipcorrection
        rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bemmodel.turbine)
    else
        rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bemmodel.turbine, tip=nothing)
    end

    # @show rvec
    # @show chordvec
    # @show twistvec

    ### Create Airfoils
    airfoils = [CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", env.rho*U*rvec[i]/env.mu, U/env.a) for i in 1:n]

    ### Create Section
    sections = CCBlade.Section.(rvec, chordvec, twistvec, airfoils)

    ### Create Operating Point
    operatingpoints = CCBlade.windturbine_op.(U, env.Omega(t), pitch, rvec, precone, yaw, tilt, azimuth, hubHt, bemmodel.shearexp, env.rho)



    #### Create GXBeam inputs
    ### Create Assembly
    assembly = create_gxbeam_assembly(gxmodel, ps)
    elements = view(assembly.elements, :) 


    ### Create Prescribed Conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))



    ##### Run initial solution
    ### Run CCBlade
    outs = CCBlade.solve.(Ref(rotor), sections, operatingpoints) 

    ### Extract CCBlade Loads
    Fz = -outs.Np #Todo: I need to interpolate the loads onto the GXBeam nodes. 
    Fy = outs.Tp

    ### Create distributed_load
    distributed_load = Dict{Int64, GXBeam.DistributedLoads{eltype(p)}}()
    f = @SVector zeros(3)
    m = @SVector zeros(3)
    m_follower = @SVector zeros(3)
    for i = 1:n #Iterate through the elements and apply the distributed load at every element. 
        f_follower = SVector(0.0, Fy[i]/elements[i].L, Fz[i]/elements[i].L) #Dividing by the length of the element so the force is distributed across the element. 
        distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end

    ### Run GXBeam
    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

    ### Obtain the deflections
    def_x = [state.elements[ielem].u[1] for ielem = 1:n]
    def_y = [state.elements[ielem].u[2] for ielem = 1:n]
    def_z = [state.elements[ielem].u[3] for ielem = 1:n]
    def_thetax = [state.elements[ielem].theta[1] for ielem = 1:n]
    Vx = [env.U(t) for ielem = 1:n]
    Vy = [-state.elements[ielem].V[2] for ielem = 1:n]


    ### Store old solution 
    Fz_new = deepcopy(Fz)
    Fz_old = deepcopy(Fz)

    Fy_new = deepcopy(Fy)
    Fy_old = deepcopy(Fy)

    deflection_new = deepcopy(def_y)
    deflection_old = deepcopy(def_y)

    # @show Fz_new
    # @show deflection_new

    ### Create solution flags
    converged_flag = false
    def_flag = false
    y_flag = false
    z_flag = false
    iter = 1

    ##### Iterate to solve
    while !converged_flag

        ### Iterate through the aerodynamic stations and solve
        ## Create Section
        r_displaced = rvec .+ def_x #Todo: Why is there deflection in the positive x? #Todo: The rvec location and the gxbeam location are not the same. 
        twist_displaced = twistvec .- def_thetax
        sects = CCBlade.Section.(r_displaced, chordvec, twist_displaced, airfoils) #Todo: Try adding the velocities. 

        ## Create Operating Point
        # ops = CCBlade.windturbine_op.(U, env.Omega(t), pitch, r_displaced, precone, yaw, tilt, azimuth, hubHt, bemmodel.shearexp, env.rho)
        for i = 1:n #Iterate throught the nodes and update the operating points
            operatingpoints[i] = CCBlade.OperatingPoint(Vx[i], Vy[i], env.rho, pitch, env.mu, env.a)
        end
        # @show operatingpoints[end].Vx operatingpoints[end].Vy

        ## Run CCBlade
        outs = CCBlade.solve.(Ref(rotor), sects, operatingpoints)

        ## Extract CCBlade Loads
        Fz_new .= -outs.Np
        Fy_new .= outs.Tp

        # @show Fz_new

        ### Solve GXBeam
        ## Update GXBeam Loads
        for i = 1:n #Iterate through the elements and apply the distributed load at every element. 
            f_follower = SVector(0.0, Fy_new[i]/elements[i].L, Fz_new[i]/elements[i].L) #Dividing by the length of the element so the force is distributed across the element. 
            distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
        end

        ## Run GXBeam
        system, converged = steady_state_analysis!(system, assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

        state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

        ## Obtain the deflections
        def_x .= [state.elements[ielem].u[1] for ielem = 1:n]
        def_y .= [state.elements[ielem].u[2] for ielem = 1:n]
        def_z .= [state.elements[ielem].u[3] for ielem = 1:n]
        def_thetax .= [state.elements[ielem].theta[1] for ielem = 1:n]

        ## Update the velocities
        for i = 1:n
            V_elem = state.elements[i].V # Element linear velocity
            Vx[i] = env.U(t) + V_elem[3] #Freestream velocity  
            Vy[i] = - V_elem[2]
        end


        ### Update solution
        deflection_new .= def_y


        ### Calculate residuals
        resids[iter, 1] = maximum(abs.(deflection_old.-deflection_new))
        resids[iter, 2] = maximum(abs.(Fy_old.-Fy_new))
        resids[iter, 3] = maximum(abs.(Fz_old.-Fz_new))

        ### Check if the solution has converged 
        if resids[iter,1]<=tolerance
            if verbose*!def_flag
                println("Y deflection converged after $iter iterations.")
            end
            def_flag = true
        else
            def_flag = false
        end

        if resids[iter,2]<=tolerance
            if verbose*!y_flag
                println("Y loadings converged after $iter iterations.")
            end
            y_flag = true
        else
            y_flag = false
        end

        if resids[iter,3]<=tolerance
            if verbose*!z_flag
                println("Z loadings converged after $iter iterations.")
            end
            z_flag = true
        else
            z_flag = false
        end

        if def_flag*y_flag*z_flag
            converged_flag = true
        end

        if iter>=maxiterations
            break
        end
        iter += 1


        ### Update old vectors
        deflection_old .= deflection_new
        Fy_old .= Fy_new
        Fz_old .= Fz_new
    end
    # iter -= 1

    return outs, state, system, assembly, prescribed_conditions, converged, iter, resids[1:iter,:] #Todo: Is the state and system that are getting passed out the state and system that were most recently calculated... or the one before the loop. I'm guessing it's the one before the loop. ... but then how/why does my dynamic solution look like it is converging to the fixed-point solution? 
end




function static_solve(fun, x0, p, t, lowbounds, upbounds; iterations=1000, dx0=zero(x0))

    solveme! = function(resids, x)
        fun(resids, dx0, x, p, t)
    end


    mcpsolve(solveme!, lowbounds, upbounds, x0, autodiff = :forward; iterations=iterations)
end