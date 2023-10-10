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
function initialize(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false, p=nothing, pfunc = (p,t) -> (;), xpfunc=nothing)

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

    t0 = first(tvec)


    # inittype = eltype(blade.airfoils.c)
    inittype = find_inittype(blade.airfoils[1].c, blade.twist[1])
    # @show typeof(blade.airfoils[1].c)
    # @show typeof(blade.twist[1])
    # @show inittype
    # @show typeof(blade.twist), typeof(blade.twist[1])



    ### ----- Prepare data storage for aerodynamic models ----- ###

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

    # A vector that CCBlade uses for solving. 
    xcc = Vector{inittype}(undef, 11)

    # Initialize DS solution
    xds, xds_idxs, p_ds = initialize_ds_model(blade.airfoils, nt; inittype)  

    # Store everything in the aerostates 
    aerostates = AeroStates(azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)

    

    ### ----- Allocate the GXBeam Data ----- ###

    system = GXBeam.DynamicSystem(assembly)
    
    gxhistory = Array{GXBeam.AssemblyState{inittype, 
        Vector{GXBeam.PointState{inittype}},
        Vector{GXBeam.ElementState{inittype}}}}(undef, nt)

    # Distributed Loads #Todo: I don't know if I'll need this. I mean the forces will be added, but they'll go through the p/pfunc portion
    if inittype <: ReverseDiff.TrackedReal
        distributed_loads = Dict{Int64, GXBeam.DistributedLoads{ReverseDiff.TrackedReal{Float64,Float64,Nothing}}}()
    else
        distributed_loads = Dict{Int64, GXBeam.DistributedLoads{inittype}}()
    end

    # Allocate the prescribed conditions 
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)) # root section is fixed

    # Get the P_augment matrix -> TODO: What exactly is the Paug matrix
    if isnothing(p)
        # parameter vector is non-existant
        paug = zeros(eltype(inittype), 12*length(assembly.points)) #Note: This was pulling type from the state vector, not the design variable vector... but I think by this point the state vector had duals in it. 
    else 
        paug = zeros(eltype(inittype), 12*length(assembly.points) + length(p))
        # copy parameter values to the augmented parameter vector
        paug[12*length(assembly.points) + 1 : 12*length(assembly.points) + length(p)] .= p
    end

    
    

    # The currently un-used GXBeam constants
    point_masses=Dict{Int,PointMass{Float64}}()
    linear_velocity=(@SVector zeros(3))
    angular_velocity=(@SVector zeros(3))
    linear_acceleration=(@SVector zeros(3)) #Todo: I'm not sure this is used. 
    angular_acceleration=(@SVector zeros(3)) #Todo: I'm not sure this is used. 
    gravity=(@SVector zeros(3))
    structural_damping=true #Todo: I might add this in later. 
    two_dimensional=false
    show_trace=false
    method=:newton
    linesearch=LineSearches.BackTracking(maxstep=1e6)
    ftol=1e-9
    iterations=1000
    xpfunc = nothing
    converged = Ref(true)

    constants = (; assembly, two_dimensional, structural_damping,
            x=system.x, resid=system.r, jacob=system.K, force_scaling = system.force_scaling,
            indices = system.indices, xpfunc, pfunc, p, t=t0, converged=converged,
            prescribed_conditions, point_masses, gravity,
            linear_velocity, angular_velocity,
            show_trace, method, linesearch, ftol, iterations)



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 


    delta = Vector{SVector{3, inittype}}(undef, na) #Todo: What is this used for? 
    def_theta = Vector{SVector{3, inittype}}(undef, na) #Todo: What is this used for
    #The structural velocities interpolated to the aerodynamic nodes.
    aeroV = Vector{SVector{3, inittype}}(undef, na)

    mesh = Mesh(interpolationpoints, delta, def_theta, aeroV, xcc, xds_idxs, p_ds, 
                system, prescribed_conditions, distributed_loads, paug, constants)

    return aerostates, gxhistory, mesh
end




"""
    initial_condition_checks(gxflag)

A function to double check that the inputs are compatible with the solution method.  
"""
function initial_condition_checks(gxflag)

    if !in(gxflag, [nothing, :steady, :spinning])
        error("The flag $gxflag isn't offered for GXBeam initialization.")
    end
end




"""
    initial_condition!()

Calculate the initial condition of the rotor. 
"""
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

        # mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc) #TODO: This shouldn't have an exclamation point in the name because it doesn't actually appear to be inplace. And it returns a struct, not a vector, so it's not an allocation... Or shouldn't be. 

        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, ccout, 1, j)
    end

    # println("dsmodel initial condition")
    dsmodel_initial_condition(aerostates, mesh, blade, rotor.turbine, t0, pitch)

    # println("extract ds loads")
    extract_ds_loads!(blade.airfoils, view(aerostates.xds, 1, :), mesh.xds_idxs, view(aerostates.phi, 1, :), mesh.p_ds, view(aerostates.cx, 1, :), view(aerostates.cy, 1, :), view(aerostates.cm, 1, :))

    # println("dimensionalize loades")
    dimensionalize!(view(aerostates.fx, 1, :), view(aerostates.fy, 1, :), view(aerostates.mx, 1, :), view(aerostates.cx, 1, :), view(aerostates.cy, 1, :), view(aerostates.cm, 1, :), blade::Blade, env::Environment, view(aerostates.W, 1, :)) 

    # println("update distributed loads")
    # @show typeof(mesh.distributed_loads)
    ReverseDiff.@skip update_forces!(mesh.distributed_loads, view(aerostates.fx, 1,:), view(aerostates.fy, 1,:), view(aerostates.mx, 1,:), blade, assembly) #Todo: Do I need the ReverseDiff skip here? I fee like this could be ... bad. I dunno. I should have notated what I was doing. 


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
    #Todo: Why am I taking an aero step after I've done the initial initialization and the GX step (I must have a reason. I'm getting pretty good results.)

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


    #Todo: I don't know why I have all these comments here. 
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


function initial_condition2(rotor::Rotor, blade::Blade, assembly::GXBeam.Assembly, env::Environment, aerostates::AeroStates, gxstates, mesh::Mesh, tvec, azimuth0, pitch; g=9.81, gxflag=nothing, structural_damping::Bool=true, linear::Bool=false, verbose::Bool=false, solver::Solver=RK4(), pfunc = (p,t) -> (;), prepp=nothing, p=nothing)
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

        # mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc) #TODO: This shouldn't have an exclamation point in the name because it doesn't actually appear to be inplace. And it returns a struct, not a vector, so it's not an allocation... Or shouldn't be. 

        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, ccout, 1, j)
    end

    # println("dsmodel initial condition")
    dsmodel_initial_condition(aerostates, mesh, blade, rotor.turbine, t0, pitch)

    # println("extract ds loads")
    extract_ds_loads!(blade.airfoils, view(aerostates.xds, 1, :), mesh.xds_idxs, view(aerostates.phi, 1, :), mesh.p_ds, view(aerostates.cx, 1, :), view(aerostates.cy, 1, :), view(aerostates.cm, 1, :))

    # println("dimensionalize loades")
    dimensionalize!(view(aerostates.fx, 1, :), view(aerostates.fy, 1, :), view(aerostates.mx, 1, :), view(aerostates.cx, 1, :), view(aerostates.cy, 1, :), view(aerostates.cm, 1, :), blade::Blade, env::Environment, view(aerostates.W, 1, :)) 

    # println("update distributed loads")
    # @show typeof(mesh.distributed_loads)
    ReverseDiff.@skip update_forces!(mesh.distributed_loads, view(aerostates.fx, 1,:), view(aerostates.fy, 1,:), view(aerostates.mx, 1,:), blade, assembly) #Todo: Do I need the ReverseDiff skip here? I fee like this could be ... bad. I dunno. I should have notated what I was doing. 


    ### GXBeam initial solution 
    # println("GXBeam initial conditions")
    if isnothing(prepp)
        p = nothing
    else
        prepp(p, aerostates, 1)
    end
    gxstates[1], system = gxbeam_initial_conditions!(env, mesh.system, assembly, mesh.prescribed_conditions, mesh.distributed_loads, t0, azimuth0, g, structural_damping, linear, gxflag, pfunc, p)

    ### Update mesh transfer variables
    # println("Update mesh")
    update_mesh!(blade, mesh, assembly, gxstates[1], env, t0, na)


    # @show length(tvec)
    # ### Trying to take the first step then set that as the first value. 
    # system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, 2, pitch; verbose, structural_damping, linear, g, solver)
    # tvi = [tvec[1], tvec[1]] #Returned nlsolve error on later time step. 
    take_aero_step!(aerostates, mesh, rotor, blade, env, tvec, 2, pitch; solver) #Todo: I think there is some overall problem in how the initial states are handled. Rather, maybe not a problem, but a difference. 
    #Todo: Why am I taking an aero step after I've done the initial initialization and the GX step (I must have a reason. I'm getting pretty good results.)

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


    #Todo: I don't know why I have all these comments here. 
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





#Todo: Shouldn't this be in the aero_only file? 
"""
    take_aero_step!()

Take an in-place step of the aerodynamic models. 

"""
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
        # mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)
        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, ccout, i, j)
        # update_aerostates!(aerostates, mesh, i, j)
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
#Todo: I need to change this function to take the old state and new state... which I think means I want to change aerostates to a vector of structs. -> What's wrong with it taking the full state vector and an index? It's faster. You can make it work with ImplicitAD onestep by putting an indexer in the parameters tuple. 
"""
    take_step!()

Take a step in the aerostructural simulation of the rotor. 
"""
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
        # mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)
        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, ccout, i, j)
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

    # @show system.x

    system, localhistory, _ = GXBeam.time_domain_analysis!(system, assembly, tvec[i-1:i]; prescribed_conditions=mesh.prescribed_conditions, distributed_loads=mesh.distributed_loads, linear, angular_velocity = Omega, reset_state=false, initial_state=gxstates[i-1], structural_damping, gravity=gravity2, pfunc=pfunc, p=p, show_trace=false) #TODO: I feel like there is a faster way to accomplish this. Like, do I really need to reallocate Omega and gravity every time step? -> Is this really a time cost though? 

    # @show system.x
     

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


function take_step2(aerostates::AeroStates, gxstates, mesh::Mesh, rotor::Rotor, blade::Blade, assembly::GXBeam.Assembly, env::Environment, system::GXBeam.DynamicSystem, tvec, i, pitch; verbose::Bool=false, speakiter::Int=100, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter, solver::Solver=RK4(), pfunc = (p,t) -> (;), prepp=nothing, p=nothing)

    na = length(blade.r)

    t = tvec[i]
    dt = tvec[i] - tvec[i-1]
    if i == 2
        dtprev = 0.
    else
        dtprev = tvec[i-1] - tvec[i-2]
    end
    
    
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
        # mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)
        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        # update_aerostates!(aerostates, mesh, i, j)
        update_aerostates!(aerostates, ccout, i, j)
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
    ReverseDiff.@skip update_forces!(mesh.distributed_loads, view(aerostates.fx, i-1, :), view(aerostates.fy, i-1, :), view(aerostates.mx, i-1, :), blade, assembly) 

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

    ### Solve GXBeam for time step 
    #TODO: This function is taking a lot of time. -> I might be able to save time by branching his code and writing another function, but most of the time is spent in nlsolve. I think all of the time spent is just time solving, not really inside of Taylor's code, but of course, if I make his code faster, then I make the solve faster. 
    # @show eltype(system)

    if isnothing(prepp) #Todo: I need to figure out how I'm passing my derivatives around. 
        p = nothing
    else
        prepp(p, aerostates, i-1)
    end
    # system, localhistory, _ = GXBeam.time_domain_analysis!(system, assembly, tvec[i-1:i]; prescribed_conditions=mesh.prescribed_conditions, distributed_loads=mesh.distributed_loads, linear, angular_velocity = Omega, reset_state=false, initial_state=gxstates[i-1], structural_damping, gravity=gravity2, pfunc=pfunc, p=p, show_trace=false) #TODO: I feel like there is a faster way to accomplish this. Like, do I really need to reallocate Omega and gravity every time step? -> Is this really a time cost though? 

    # current_state = gxstates[i-1]
    
    constants = (;mesh.constants..., t, dt, dtprev, x=system.x, p, distributed_loads=mesh.distributed_loads, gravity=gravity2, angular_velocity=Omega)

    #Todo: Do I need to update p and paug? 
    paug = mesh.paug

    # @show system.x #Matched for the first time step (after initialization). but not during the second solution. 

    system, xgx, converged = GXBeam.take_step(system, gxstates[i-1], constants, paug;) 

    # @show system.x #Doesn't match after the first time step. #Todo: 

    gxstates[i] = GXBeam.newmark_output(system, xgx, paug, constants)
     

    #TODO: Can I save memory by directly allocating to the gxstates vector? -> I can probably save allocations by augmenting the time_domain_analysis!() function to already have the results allocated. 

    ### Extract GXBeam outputs
    # gxstates[i] = localhistory[end] 

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

    aerostates, gxstates, mesh = initialize(blade, assembly, tvec; verbose)

    
    system = initial_condition!(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch; verbose)

    for i = 2:nt
        system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch; verbose, speakiter, plotiter, plotbool, structural_damping, linear, g, solver)
    end

    return aerostates, gxstates
end

function simulate2(rotor::Rotors.Rotor, blade::Blade, env::Environment, assembly::GXBeam.Assembly, tvec; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, warnings::Bool=true, azimuth0=0.0, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter)

    nt = length(tvec)

    aerostates, gxstates, mesh = initialize(blade, assembly, tvec; verbose)

    
    system = initial_condition2(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch; verbose)

    for i = 2:nt
        system = take_step2(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch; verbose, speakiter, plotiter, plotbool, structural_damping, linear, g, solver)
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

