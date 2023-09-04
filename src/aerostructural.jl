#=


=#

export initialize, initial_condition!, take_step!, simulate, simulate!






function checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    if minimum(rvec)<=rhub
        @warn("A member of rvec is less than or equal to rhub. This will cause problems with the hub corrections.")
    end
    if maximum(rvec)>=rtip
        @warn("A member of rvec is greater than or equal to rtip. This will cause problems with the tip corrections.")
    end
    if pitch>(pi/4)
        @warn("The pitch value is greater than 45 degrees. Ensure that the pitch is in radians and not degrees.")
    end
    if precone>(pi/4)
        @warn("The precone value is greater than 45 degrees. Ensure that the precone is in radians and not degrees.")
    end
    if tilt>(pi/4)
        @warn("The tilt value is greater than 45 degrees. Ensure that the tilt is in radians and not degrees.")
    end
    if yaw>(pi/4)
        @warn("The yaw value is greater than 45 degrees. Ensure that the yaw is in radians and not degrees.")
    end
    if maximum(twistvec)>(pi/2)
        @warn("The maximum twist value is greater than 90 degrees. Ensure that the twist distribution is given in radians.")
    end
end

#TODO: It might be a good idea to make a version that is completely in place. (Pass in data storage and riso ode.)

function find_inittype(vars...)
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


"""
    initialize()

Prepare data structures for the simulation. 


"""
function initialize(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false)

    if verbose
        println("Rotors.jl initializing solution...")
    end

    # if warnings
    #     checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    # end

    
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types.

    ### Initialization information
    na = length(blade.rR)
    nt = length(tvec)


    # inittype = eltype(blade.airfoils.c)
    inittype = find_inittype(blade.airfoils[1].c, blade.twist[1])
    @show typeof(blade.airfoils[1].c)
    @show typeof(blade.twist[1])
    @show inittype
    # @show typeof(blade.twist), typeof(blade.twist[1])



    ### Prepare data storage
    azimuth = Array{inittype}(undef, nt)
    phi = Array{inittype}(undef,(nt, na))
    alpha = Array{inittype}(undef,(nt, na))
    W = Array{inittype}(undef,(nt, na))

    Cx = Array{inittype}(undef,(nt, na))
    Cy = Array{inittype}(undef,(nt, na))
    Cm = Array{inittype}(undef,(nt, na))

    Fx = Array{inittype}(undef,(nt, na))
    Fy = Array{inittype}(undef,(nt, na))
    Mx = Array{inittype}(undef,(nt, na))

     
    # cchistory = Vector{CCBlade.Outputs{inittype}}(undef, na) 
    cchistory = Vector{CCBlade.Outputs{inittype}}(undef, na) 
    # @show inittype
    xcc = Vector{inittype}(undef, 11)

    gxhistory = Array{GXBeam.AssemblyState{inittype, Vector{GXBeam.PointState{inittype}}, Vector{GXBeam.ElementState{inittype}}}}(undef, nt)



    ### Initialize DS solution
    xds, xds_idxs, p_ds = initialize_ds_model(blade.airfoils, nt; inittype)  

     
    
    ### Store everything in the aerostates 
    aerostates = AeroStates(azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)

    

    ### Allocate the distributed load 
    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{inittype}}()
    

    ### Allocate the prescribed conditions 
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)) # root section is fixed



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 


    delta = Vector{SVector{3, inittype}}(undef, na)
    def_theta = Vector{SVector{3, inittype}}(undef, na)
    aeroV = Vector{SVector{3, inittype}}(undef, na)

    mesh = Mesh(interpolationpoints, prescribed_conditions, distributed_loads, delta, def_theta, aeroV, cchistory, xcc, xds_idxs, p_ds)

    return aerostates, gxhistory, mesh
end

function initial_condition_checks(gxflag)

    if !in(gxflag, [nothing, :steady, :spinning])
        error("The flag $gxflag isn't offered for GXBeam initialization.")
    end
end

function initial_condition!(rotor::Rotor, blade::Blade, assembly::GXBeam.Assembly, env::Environment, aerostates::AeroStates, gxstates, mesh::Mesh, tvec, azimuth0, pitch; g=9.81, gxflag=nothing, structural_damping::Bool=true, linear::Bool=false, verbose::Bool=false, solver::Solver=RK4(), pfunc = (p,t) -> (;), prepp=nothing, p=nothing)
    #TODO: Maybe include the coupled aerostructural solution as an option? 

    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. -> These might best go in the blade constructor. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). -> I'm not entirely sure what I meant by that. 

    t0 = tvec[1]

    aerostates.azimuth[1] = azimuth0

    na = length(blade.r)

    if verbose
        println("Calculating initial condition...")
    end

    # println("checks")
    initial_condition_checks(gxflag)

    ### Initialize BEM solution 
    # println("BEM solution")
    for j = 1:na
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth0)

        mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, mesh, 1, j)
    end

    # println("dsmodel initial condition")
    dsmodel_initial_condition(aerostates, mesh, blade, rotor.turbine, t0, pitch)

    # println("extract ds loads")
    extract_ds_loads!(blade.airfoils, view(aerostates.xds, 1, :), mesh.xds_idxs, view(aerostates.phi, 1, :), mesh.p_ds, view(aerostates.cx, 1, :), view(aerostates.cy, 1, :), view(aerostates.cm, 1, :))

    # println("dimensionalize loades")
    dimensionalize!(view(aerostates.fx, 1, :), view(aerostates.fy, 1, :), view(aerostates.mx, 1, :), view(aerostates.cx, 1, :), view(aerostates.cy, 1, :), view(aerostates.cm, 1, :), blade::Blade, env::Environment, view(aerostates.W, 1, :)) 

    # println("update distributed loads")
    update_forces!(mesh.distributed_loads, view(aerostates.fx, 1,:), view(aerostates.fy, 1,:), view(aerostates.mx, 1,:), blade, assembly)


    ### GXBeam initial solution 
    # println("GXBeam initial conditions")
    if isnothing(prepp)
        p = nothing
    else
        prepp(p, aerostates, 1)
    end
    gxstates[1], system = gxbeam_initial_conditions(env, assembly, mesh.prescribed_conditions, mesh.distributed_loads, t0, azimuth0, g, structural_damping, linear, gxflag, pfunc, p)

    ### Update mesh transfer variables
    # println("Update mesh")
    update_mesh!(blade, mesh, assembly, gxstates[1], env, t0, na)


    # @show length(tvec)
    # ### Trying to take the first step then set that as the first value. 
    # system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, 2, pitch; verbose, structural_damping, linear, g, solver)
    # tvi = [tvec[1], tvec[1]] #Returned nlsolve error on later time step. 
    take_aero_step!(aerostates, mesh, rotor, blade, env, tvec, 2, pitch; solver) #Todo: I think there is some overall problem in how the initial states are handled. Rather, maybe not a problem, but a difference. 

    for j = 1:na
        aerostates.W[1,j] = aerostates.W[2,j]
        aerostates.phi[1,j] = aerostates.phi[2,j]
        aerostates.xds[1,j] = aerostates.xds[2,j]
        
        aerostates.cx[1,j] = aerostates.cx[2,j]
        aerostates.cy[1,j] = aerostates.cy[2,j]
        aerostates.cm[1,j] = aerostates.cm[2,j]

        aerostates.fx[1,j] = aerostates.fx[2,j]
        aerostates.fy[1,j] = aerostates.fy[2,j]
        aerostates.mx[1,j] = aerostates.mx[2,j]
    end

    # update_mesh!(blade, mesh, assembly, gxstates[1], env, t0, na)

    # gxstates[1] = gxstates[2] #Todo: Do I need to reset the system for initialization? Do I even want to be updating? If I'm never updating the structural? 

    # # println("update distributed loads")
    # update_forces!(mesh.distributed_loads, view(aerostates.fx, 1,:), view(aerostates.fy, 1,:), view(aerostates.mx, 1,:), blade, assembly)


    # ### GXBeam initial solution 
    # # println("GXBeam initial conditions")
    # gxstates[1], system = gxbeam_initial_conditions(env, assembly, mesh.prescribed_conditions, mesh.distributed_loads, t0, azimuth0, g, structural_damping, linear, gxflag)

    # ### Update mesh transfer variables
    # # println("Update mesh")
    # update_mesh!(blade, mesh, assembly, gxstates[1], env, t0, na)



    return system
end

function take_aero_step!(aerostates::AeroStates, mesh::Mesh, rotor::Rotor, blade::Blade, env::Environment, tvec, i, pitch; solver::Solver=RK4(), pfunc = (p,t) -> (;), prepp=nothing, p=nothing)
    na = length(blade.r)

    t = tvec[i]
    dt = tvec[i] - tvec[i-1]
    
    
    if dt<0
        error("Time step is negative")
    end

    #update azimuthal position
    aerostates.azimuth[i] = env.RS(t)*dt + aerostates.azimuth[i-1] #Euler step for azimuthal position. #TODO: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

    if aerostates.azimuth[i]<aerostates.azimuth[i-1]
        @warn("Blade moved backwards")
    end

    ### Update BEM inputs and solve
    for j = 1:na
        ### Update base inflow velocities
        Vx, Vy = Rotors.get_aerostructural_velocities(rotor, blade, env, t, j, aerostates.azimuth[i], mesh.delta[j], mesh.def_theta[j], mesh.aerov[j])
        
        #TODO: Write a solver that is initialized with the previous inflow angle.
        mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, mesh, i, j)
    end
    
    ### Update Dynamic Stall model inputs 
    update_ds_inputs!(blade.airfoils, view(mesh.p_ds, :), view(aerostates.W, i, :), view(aerostates.phi, i, :), blade.twist, pitch, dt, rotor.turbine)
    
    ### Integrate Dynamic Stall model
    update_ds_states!(solver, blade.airfoils, view(aerostates.xds, i-1, :), view(aerostates.xds, i, :), mesh.xds_idxs, mesh.p_ds, t, dt)

    ### Extract loads 
    extract_ds_loads!(blade.airfoils, view(aerostates.xds, i, :), mesh.xds_idxs, view(aerostates.phi, i, :), mesh.p_ds, view(aerostates.cx, i, :), view(aerostates.cy, i, :), view(aerostates.cm, i, :))
 
    
    dimensionalize!(view(aerostates.fx, i, :), view(aerostates.fy, i, :), view(aerostates.mx, i, :), view(aerostates.cx, i, :), view(aerostates.cy, i, :), view(aerostates.cm, i, :), blade::Blade, env::Environment, view(aerostates.W, i, :))
    #These loads do not need to be rotated because they will be applied in the deflected frame (a follower load). This should also be true for things like precone, tilt, and yaw if they are defined correctly in GXBeam. 
end



#TODO: Function headers
function take_step!(aerostates::AeroStates, gxstates, mesh::Mesh, rotor::Rotor, blade::Blade, assembly::GXBeam.Assembly, env::Environment, system::GXBeam.DynamicSystem, tvec, i, pitch; verbose::Bool=false, speakiter::Int=100, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter, solver::Solver=RK4(), pfunc = (p,t) -> (;), prepp=nothing, p=nothing)

    na = length(blade.r)

    t = tvec[i]
    dt = tvec[i] - tvec[i-1]
    
    
    if dt<0
        error("Time step is negative")
    end

    #update azimuthal position
    aerostates.azimuth[i] = env.RS(t)*dt + aerostates.azimuth[i-1] #Euler step for azimuthal position. #TODO: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

    if aerostates.azimuth[i]<aerostates.azimuth[i-1]
        @warn("Blade moved backwards")
    end

    ### Update BEM inputs and solve
    for j = 1:na
        ### Update base inflow velocities
        Vx, Vy = Rotors.get_aerostructural_velocities(rotor, blade, env, t, j, aerostates.azimuth[i], mesh.delta[j], mesh.def_theta[j], mesh.aerov[j])
        
        #TODO: Write a solver that is initialized with the previous inflow angle.
        mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, mesh, i, j)
    end
    
    ### Update Dynamic Stall model inputs 
    update_ds_inputs!(blade.airfoils, view(mesh.p_ds, :), view(aerostates.W, i, :), view(aerostates.phi, i, :), blade.twist, pitch, dt, rotor.turbine)
    
    ### Integrate Dynamic Stall model
    update_ds_states!(solver, blade.airfoils, view(aerostates.xds, i-1, :), view(aerostates.xds, i, :), mesh.xds_idxs, mesh.p_ds, t, dt)

    ### Extract loads 
    extract_ds_loads!(blade.airfoils, view(aerostates.xds, i, :), mesh.xds_idxs, view(aerostates.phi, i, :), mesh.p_ds, view(aerostates.cx, i, :), view(aerostates.cy, i, :), view(aerostates.cm, i, :))
 
    
    dimensionalize!(view(aerostates.fx, i, :), view(aerostates.fy, i, :), view(aerostates.mx, i, :), view(aerostates.cx, i, :), view(aerostates.cy, i, :), view(aerostates.cm, i, :), blade::Blade, env::Environment, view(aerostates.W, i, :))
    #These loads do not need to be rotated because they will be applied in the deflected frame (a follower load). This should also be true for things like precone, tilt, and yaw if they are defined correctly in GXBeam. 
    
    
    ### Update GXBeam loads  
    update_forces!(mesh.distributed_loads, view(aerostates.fx, i-1, :), view(aerostates.fy, i-1, :), view(aerostates.mx, i-1, :), blade, assembly) 

    Omega = SVector(0.0, 0.0, -env.RS(t))
    # gravity = SVector(-g*cos(aerostates.azimuth[i-1]), -g*sin(aerostates.azimuth[i-1]), 0.0) #TODO: I need to include tilt, and precone here. 

    # @show typeof(aerostates.azimuth)

    if isnothing(p)
        a0 = aerostates.azimuth[i-1]
        a1 = aerostates.azimuth[i]
    else
        a0 = aerostates.azimuth[i-1].value
        a1 = aerostates.azimuth[i].value
    end

    gravity2 = (tee) -> SVector(-g*cos((a0*(tvec[i]-tee) + a1*(tee-tvec[i-1]))/(tvec[i]-tvec[i-1])), -g*sin((a0*(tvec[i]-tee) + a1*(tee-tvec[i-1]))/(tvec[i]-tvec[i-1])), 0.0) ##Todo t = tvec[i].... So this be way wrong. Oh... this is an inline function so the solver can linearly interpolate the gravity vector across time. But... I think the time domain analysis only analyzes at the given time steps... which means that this function doesn't get called really... I dunno. 

    # @show typeof(gravity2)

    #Note: Taylor applies the gravitational load by C'*mass*C*gvec

    ### Solve GXBeam for time step #TODO: This function is taking a lot of time. -> I might be able to save time by branching his code and writing another function, but most of the time is spent in nlsolve. I think all of the time spent is just time solving, not really inside of Taylor's code, but of course, if I make his code faster, then I make the solve faster. 
    # @show eltype(system)

    if isnothing(prepp)
        p = nothing
    else
        prepp(p, aerostates, i-1)
    end
    system, localhistory, _ = GXBeam.time_domain_analysis!(system, assembly, tvec[i-1:i]; prescribed_conditions=mesh.prescribed_conditions, distributed_loads=mesh.distributed_loads, linear, angular_velocity = Omega, reset_state=false, initial_state=gxstates[i-1], structural_damping, gravity=gravity2, pfunc=pfunc, p=p, show_trace=false) #TODO: I feel like there is a faster way to accomplish this. Like, do I really need to reallocate Omega and gravity every time step? -> Is this really a time cost though? 
     

    #TODO: Can I save memory by directly allocating to the gxstates vector? -> I can probably save allocations by augmenting the time_domain_analysis!() function to already have the results allocated. 

    ### Extract GXBeam outputs
    gxstates[i] = localhistory[end] 

    ### Update aero inputs from structures.
    update_mesh!(blade, mesh, assembly, gxstates[i], env, t, na)



    if verbose & (mod(i-1, speakiter)==0)
        println("")
        println("Simulation time: ", t)
    end



    if plotbool & (mod(i-1, plotiter)==0)
        # tipdef_x = [gxhistory[k].points[end].u[1] for k in eachindex(tvec[1:i])]
        # tipdef_y = [gxhistory[k].points[end].u[2] for k in eachindex(tvec[1:i])]
        # tipdef_z = [gxhistory[k].points[end].u[3] for k in eachindex(tvec[1:i])]
        # plt = plot(xaxis="Time (s)", yaxis="Tip Deflection", legend=:outerright)
        # plot!(tvec[1:i], tipdef_x, lab="X deflection")
        # plot!(tvec[1:i], tipdef_y, lab="Y deflection")
        # plot!(tvec[1:i], tipdef_z, lab="Z deflection")
        # display(plt)

        thetamat = zeros(i, 3)

        for k in 1:i
            thetamat[k,:] = gxhistory[k].points[end].theta
        end

        plt = plot(xaxis="Time (s)", yaxis="theta def", legend=:outerright)
        plot!(tvec[1:i], thetamat[:,1], lab="X deflection")
        plot!(tvec[1:i], thetamat[:,2], lab="Y deflection")
        plot!(tvec[1:i], thetamat[:,3], lab="Z deflection")
        display(plt)
    end

    return system #Todo: Should this be returned? (For scoping and passing of data)
end

function simulate(rotor::Rotors.Rotor, blade::Blade, env::Environment, assembly::GXBeam.Assembly, tvec; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, warnings::Bool=true, azimuth0=0.0, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter)

    nt = length(tvec)

    aerostates, gxstates, mesh = initialize(blade, assembly, tvec; verbose=true)

    
    system = initial_condition!(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch; verbose)

    for i = 2:nt
        system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch; verbose, speakiter, plotiter, plotbool, structural_damping, linear, g, solver)
    end

    return aerostates, gxstates
end


function simulate!(rotor::Rotors.Rotor, blade::Blade, env::Environment, assembly::GXBeam.Assembly, tvec, aerostates::AeroStates, gxstates, mesh::Mesh; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, warnings::Bool=true, azimuth0=0.0, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter) #Todo: Move g to the environment struct. 

    nt = length(tvec)

    system = initial_condition!(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch; verbose)

    for i = 2:nt
        system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch; verbose, speakiter, plotiter, plotbool, structural_damping, linear, g, solver) #Todo: I want to change take_step! to rely on t and dt... but I'm not sure how to access the data structures then... because I've condensed the data into the aerostates, gxstates and they require an integer to index. 
    end

    return aerostates, gxstates
end

#TODO: Make a memory efficient take_step!() function that only saves certain time indices. 

