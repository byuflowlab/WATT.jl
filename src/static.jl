#=
Adam Cardoza 

A fixed point solver that iterates back and forth between CCBlade and GXBeam until it converges on a solution. A dynamic stall model isn't included, so the solver assumes that the blade is at some sort of non-oscillitory steady state. 

- 6/?/22 - Initial creation

- 7/29/22 - Reworked to have the assembly and design vectors passed in rather than be created inside the function. (We can always do that inside another function for optimization.) Also working on having the loads be interpolated. Note that it appears that the deflections drastically decreased



TODO: 
- Add asymmetric environment variables. 
        - add tilt
        - add yaw
        - add precone
        - add azimuth

=#

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
function fixedpoint(rvec, chordvec, twistvec, pitch, rhub, rtip, assembly::GXBeam.Assembly, env, blade, B;t=0.0, maxiterations=1000, tolerance=1e-3, g=0.0, verbose=false, speakiterations=10, turbine::Bool=true, tip::CCBlade.TipCorrection=PrandtlTipHub(), precone=0.0, yaw=0.0, tilt=0.0, azimuth=0.0)

    ############ Initial read in and preparations 
    if verbose
        println("Preparing Solution")
    end

    ### Extract parameters
    na = length(rvec) #Number of Aerodynamic nodes
    ne = length(assembly.elements) #Number of structural nodes

    rvec_gxbeam = [assembly.elements[i].x[1] for i in 1:ne]


    ### Environmental variables
    Omega = SVector(0.0, 0.0, -env.RS(t))
    grav = SVector(0.0, -g, 0.0)
    U = env.Vinf(t)


    ### Create storage structures for solution
    resids = zeros(maxiterations, 3) #columns of resids is deflection y, force y, force z






    ############## Run CCBlade ############
    ### Create Rotor
    rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=turbine, tip=tip)


    ### Create Section
    sections = [Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i in 1:na] 

    ### Create Operating Point
    operatingpoints = [CCBlade.OperatingPoint(U, env.RS(t)*rvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na] 

    ### Run CCBlade
    outs = CCBlade.solve.(Ref(rotor), sections, operatingpoints)  
    outshistory = [outs]







    ############## Run GXBeam  ##################

    ### Create Prescribed Conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Fixed root condition. 

    
    ### Extract CCBlade Loads
    Fzfit = Akima(rvec, -outs.Np) 
    Fyfit = Akima(rvec, outs.Tp)

    ### Create distributed_load
    distributed_load = Dict{Int64, GXBeam.DistributedLoads{eltype(rvec)}}()
    f_follower = @SVector zeros(3)
    m = @SVector zeros(3)
    m_follower = @SVector zeros(3)
    for i = 1:ne #Iterate through the elements and apply the distributed load at every element.
        f = SVector(0.0, Fyfit(rvec_gxbeam[i]), Fzfit(rvec_gxbeam[i])) # Update the distributed dead load (as opposed to a follower load). 
        distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end

    ### Run GXBeam
    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)



    ############## Prepare to iterate between CCBlade and GXBeam ################
    ### Create vectors to update across iterations so that I can pass the solution out of the iteration vector. 
    history = [state]
    syshistory = [system]

    ### Obtain the deflections
    def_x = [history[1].elements[ielem].u[1] for ielem = 1:ne]  #Linear deflection in the freestream direction
    def_y = [history[1].elements[ielem].u[2] for ielem = 1:ne]  #Linear deflection in the leadlag direction
    def_z = [history[1].elements[ielem].u[3] for ielem = 1:ne]  #Linear deflection in the flapwise direction
    def_thetax = [history[1].elements[ielem].theta[1] for ielem = 1:ne] #Pitching deflection


    interpolationpoints = create_interpolationpoints(assembly, rvec) #Create points to interpolate the velocities from
    aeroV = [interpolate_velocity(interpolationpoints[j], assembly, state) for j = 1:na] #Interpolate the velocities from the elements to the aerodynamic nodes
    aero_thetax =  [interpolate_angle(interpolationpoints[j], assembly, state) for j=1:na]

    Vx = [env.Vinf(t) + aeroV[iaero][3] for iaero = 1:na]  
    Vy = [-aeroV[iaero][2] for iaero = 1:na]


    ### Store old solution 
    Fz_new = deepcopy(-outs.Np) 
    Fz_old = deepcopy(-outs.Np)

    Fy_new = deepcopy(outs.Tp)
    Fy_old = deepcopy(outs.Tp)

    deflection_new = deepcopy(def_y)
    deflection_old = deepcopy(def_y)



    ### Create solution flags
    converged_flag = false
    def_flag = false
    y_flag = false
    z_flag = false
    iter = 1

    ##### Iterate to solve 
    if verbose
        println("Beginning Iterations")
    end
    while !converged_flag
        if verbose&&(mod(iter,speakiterations)==0)
            println("Beginning Iteration $iter")
        end

        ############### Iterate through the aerodynamic stations and solve
        ## Create Section
        twist_displaced = twistvec .- aero_thetax #deform the twist distribution.  
        sects = CCBlade.Section.(rvec, chordvec, twist_displaced, blade.airfoils) #Update the section. 
        #I don't displace the radial vector to avoid issues with the hub and tip corrections. The effects of their displacement is felt through the X and Y velocities. 
        

        ## Create Operating Point
        for i = 1:na #Iterate throught the nodes and update the operating points
            operatingpoints[i] = CCBlade.OperatingPoint(Vx[i], Vy[i], env.rho, pitch, env.mu, env.a)
        end

        ## Run CCBlade
        outs = CCBlade.solve.(Ref(rotor), sects, operatingpoints)

        ## Extract CCBlade Loads
        Fz_new .= -outs.Np #TODO: Does the dot operator do anything in this case? 
        Fy_new .= outs.Tp

        outshistory[1] = outs





        ##########  Solve GXBeam
        ## Update GXBeam Loads
        Fyfit_local = Akima(rvec, Fy_new)
        Fzfit_local = Akima(rvec, Fz_new)
        for i = 1:ne #Iterate through the elements and apply the distributed load at every element. 
            f = SVector(0.0, Fyfit_local(rvec_gxbeam[i]), Fzfit_local(rvec_gxbeam[i])) #Update the distributed dead load
            distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
        end

        ## Run GXBeam
        syshistory[1], converged = steady_state_analysis!(system, assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

        history[1] = AssemblyState(syshistory[1], assembly; prescribed_conditions = prescribed_conditions)






        ########## Update the outputs and prepare for the next iteration
        ## Update the deflections
        def_x .= [history[1].elements[ielem].u[1] for ielem = 1:ne]
        def_y .= [history[1].elements[ielem].u[2] for ielem = 1:ne]
        def_z .= [history[1].elements[ielem].u[3] for ielem = 1:ne]
        def_thetax .= [history[1].elements[ielem].theta[1] for ielem = 1:ne]

        ## Update the velocities
        for iaero = 1:na
            aeroV[iaero] = interpolate_velocity(interpolationpoints[iaero], assembly, history[1])
            aero_thetax[iaero] =  interpolate_angle(interpolationpoints[iaero], assembly, history[1])
            Vx[iaero] = env.Vinf(t) + aeroV[iaero][3] #Freestream velocity  
            Vy[iaero] = -aeroV[iaero][2]
        end


        ### Update solution
        deflection_new .= def_y



        ############ Calculate residuals
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

    return outshistory[1], history[1], syshistory[1], converged_flag, iter, resids[1:iter,:] 
end




function static_solve(fun, x0, p, t, lowbounds, upbounds; iterations=1000, dx0=zero(x0))

    solveme! = function(resids, x)
        fun(resids, dx0, x, p, t)
    end


    mcpsolve(solveme!, lowbounds, upbounds, x0, autodiff = :forward; iterations=iterations)
end