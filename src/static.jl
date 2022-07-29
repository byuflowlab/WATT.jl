

"""
    fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations=1000, tolerance=1e-3, g=0.0, verbose=false)

Iterates between the BEM (CCBlade) and structural model (GXBeam) until the steady state solution converges to a given tolerance. The convergence criteria is Y deflections and the Y and Z loadings. 

### Inputs
- bemmodel::BEM - The bem model. Doesn't do much right now because I don't allow for propeller analysis yet. 
- gxmodel::gxbeam - The structural model. Holds the x location for the points and elements. 
- env::Environment - the environment model
- blade::Blade - the blade model
- p::Array{TF, 1} - Array of parameters for the aerodynamic and structural models. It goes [p_a, p_s]. There are 7*n aero parameters (where n is the number of nodes). The rest are structural nodes. 
- maxiterations::TI - The maximum number of iterations that the solver will take. 
- tolerance::TF - The accuracy of decimals that the user requires. 
- g::TF - The amount of gravity applied. Currently defaults to 0.0, but one day will default to 9.817
- verbose::Bool - Whether or not to print out statements showing the convergence progress. 

### Outputs
- outs::CCBlade.Outputs - The final CCBlade output struct
- state::GXBeam.State - The final GXBeam state
- system::GXBeam.System
- assembly::GXBeam.Assembly
- converged::Bool - Whether or not the solver converged within the maximum number of iterations.
- iterations::TI - The number of iterations that the solver required to converge the residuals. 
- residuals::Array{TF, 2} - a 3 x j matrix of the residuals at every iteration. Column 1 is the y deflection residuals, Column 2 is the y force residuals, and Column 3 is the z force residuals. j is the number of iterations taken by the solver. 


### Notes
- Currently yaw, tilt, and azimuth aren't included. 
"""
function fixedpoint(bemmodel, gxmodel, env, blade, p, B; maxiterations=1000, tolerance=1e-3, g=0.0, verbose=false)


    ### Seperate the parameters and read in any pertainent data
    n = gxmodel.ne
    pa = view(p, 1:(7*n))
    ps = view(p, (7*n+1):length(p))

    ### Constant problem dependent parameters
    pitch = pa[4]
    rhub = pa[5]
    rtip = pa[6]
    hubHt = pa[7]

    ### Environmental variables
    Omega = SVector(0.0, 0.0, -env.RS(0.0))
    grav = SVector(0.0, -g, 0.0)
    U = env.Vinf(0.0)


    ### Base vectors for turbine. 
    radii_idx = 1:7:(7*n)
    chord_idx = 2:7:(7*n)
    twist_idx = 3:7:(7*n)
    rvec = view(pa, radii_idx)  
    chordvec = view(pa, chord_idx)
    twistvec = view(pa, twist_idx)

    ### Unchanging constants
    t = 0.0
    precone = 0.0 #Precone. -> I'm going to change the radius manually. -> I could have a base precone, then allow for changes on top of that? 
    # B = 1 #Number of blades #Note: You can't just default to 1. B affects the Prandtl Tip/hub correction. 
    yaw = 0.0 #TODO: Need to include these. 
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


    ### Create Section
    sections = [Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i in 1:n] 

    ### Create Operating Point
    operatingpoints = [CCBlade.OperatingPoint(U, env.RS(t)*rvec[i], env.rho, pitch, env.mu, env.a) for i in 1:n] 


    #### Create GXBeam inputs
    ### Create Assembly
    assembly = create_gxbeam_assembly(gxmodel, ps)
    # elements = view(assembly.elements, :) 


    ### Create Prescribed Conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))



    ##### Run initial solution
    ### Run CCBlade
    outs = CCBlade.solve.(Ref(rotor), sections, operatingpoints)  
    outshistory = [outs]

    ### Extract CCBlade Loads
    Fz = -outs.Np  
    Fy = outs.Tp

    ### Create distributed_load
    distributed_load = Dict{Int64, GXBeam.DistributedLoads{eltype(p)}}()
    f_follower = @SVector zeros(3)
    m = @SVector zeros(3)
    m_follower = @SVector zeros(3)
    for i = 1:n #Iterate through the elements and apply the distributed load at every element.
        f = SVector(0.0, Fy[i], Fz[i]) # Update the distributed dead load (as opposed to a follower load). 
        distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end

    ### Run GXBeam
    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

    ### Create vectors to update across iterations so that I can pass the solution out of the iteration vector. 
    history = [state]
    syshistory = [system]

    ### Obtain the deflections
    def_x = [history[1].elements[ielem].u[1] for ielem = 1:n]
    def_y = [history[1].elements[ielem].u[2] for ielem = 1:n]
    def_z = [history[1].elements[ielem].u[3] for ielem = 1:n]
    def_thetax = [history[1].elements[ielem].theta[1] for ielem = 1:n]


    Vx = [env.Vinf(t) for ielem = 1:n]
    Vy = [-history[1].elements[ielem].V[2] for ielem = 1:n]


    ### Store old solution 
    Fz_new = deepcopy(Fz)
    Fz_old = deepcopy(Fz)

    Fy_new = deepcopy(Fy)
    Fy_old = deepcopy(Fy)

    deflection_new = deepcopy(def_y)
    deflection_old = deepcopy(def_y)



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
        twist_displaced = twistvec .- def_thetax
        sects = CCBlade.Section.(rvec, chordvec, twist_displaced, blade.airfoils) #Update the section. 
        #I don't displace the radial vector to avoid issues with the hub and tip corrections. The effects of their displacement is felt through the X and Y velocities. 
        

        ## Create Operating Point
        for i = 1:n #Iterate throught the nodes and update the operating points
            operatingpoints[i] = CCBlade.OperatingPoint(Vx[i], Vy[i], env.rho, pitch, env.mu, env.a)
        end

        ## Run CCBlade
        outs = CCBlade.solve.(Ref(rotor), sects, operatingpoints)

        ## Extract CCBlade Loads
        Fz_new .= -outs.Np #TODO: Does the dot operator do anything in this case? 
        Fy_new .= outs.Tp

        outshistory[1] = outs


        ### Solve GXBeam
        ## Update GXBeam Loads
        for i = 1:n #Iterate through the elements and apply the distributed load at every element. 
            f = SVector(0.0, Fy_new[i], Fz_new[i]) #Update the distributed dead load
            distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
        end

        ## Run GXBeam
        syshistory[1], converged = steady_state_analysis!(system, assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

        history[1] = AssemblyState(syshistory[1], assembly; prescribed_conditions = prescribed_conditions)

        ## Obtain the deflections
        def_x .= [history[1].elements[ielem].u[1] for ielem = 1:n]
        def_y .= [history[1].elements[ielem].u[2] for ielem = 1:n]
        def_z .= [history[1].elements[ielem].u[3] for ielem = 1:n]
        def_thetax .= [history[1].elements[ielem].theta[1] for ielem = 1:n]

        ## Update the velocities
        for i = 1:n
            V_elem = history[1].elements[i].V # Element linear velocity
            Vx[i] = env.Vinf(t) + V_elem[3] #Freestream velocity  #Todo: This assumes that the structural and aerodynamic meshs are co-located. I need to switch to an interpolated method. 
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
            break
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

    return outshistory[1], history[1], syshistory[1], assembly, converged, iter, resids[1:iter,:] 
end




function static_solve(fun, x0, p, t, lowbounds, upbounds; iterations=1000, dx0=zero(x0))

    solveme! = function(resids, x)
        fun(resids, dx0, x, p, t)
    end


    mcpsolve(solveme!, lowbounds, upbounds, x0, autodiff = :forward; iterations=iterations)
end