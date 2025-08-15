#=


=#

export initialize, initial_condition, take_step!, simulate, simulate!






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
function initialize(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false, p=nothing, pfunc = (p,t) -> (;), xpfunc=nothing, structural_damping::Bool=true, linear::Bool=false)

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

    inittype = find_inittype(blade.c[1], blade.twist[1])



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
    xds, xds_idxs, y_ds = initialize_ds_model(blade.airfoils, nt; inittype)  

    # Store everything in the aerostates 
    # aerostates = AeroStates(azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)
    aerostates = (;azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)

    

    ### ----- Allocate the GXBeam Data ----- ###

    system = GXBeam.DynamicSystem(assembly)

    ngx = length(system.x)

    gxstates = Array{inittype}(undef, (nt, 2*ngx))
    
    gxhistory = Array{GXBeam.AssemblyState{inittype, 
        Vector{GXBeam.PointState{inittype}},
        Vector{GXBeam.ElementState{inittype}}}}(undef, nt)

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{Float64}}()

    # Allocate the prescribed conditions 
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)) # root section is fixed

    
    

    # The currently un-used GXBeam constants
    point_masses=Dict{Int,PointMass{Float64}}()
    linear_velocity=(@SVector zeros(3))
    angular_velocity=(@SVector zeros(3))
    # structural_damping=true #todo: I might add this in later. 
    two_dimensional=false
    xpfunc = nothing
    



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 


    delta = Vector{SVector{3, inittype}}(undef, na) #todo: What is this used for? 
    def_theta = Vector{SVector{3, inittype}}(undef, na) #todo. What is this used for
    #The structural velocities interpolated to the aerodynamic nodes.
    aerov = Vector{SVector{3, inittype}}(undef, na)

    
    mesh = (; interpolationpoints, delta, def_theta, aerov, xcc, xds_idxs, y_ds,
                assembly, system, prescribed_conditions, distributed_loads,
                point_masses, linear_velocity, angular_velocity,
                xpfunc, pfunc, two_dimensional, structural_damping, linear)

    return aerostates, gxstates, gxhistory, mesh
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
function initial_condition!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, gxstates, azimuth0, rotor::Rotor, env::Environment, blade, mesh, tvec, pitch; g=9.81, gxflag=nothing, verbose::Bool=false, solver::Solver=RK4(), prepp=nothing, p=nothing)
    #TODO: Maybe include the coupled aerostructural solution as an option? 

    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. -> These might best go in the blade constructor. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). -> I'm not entirely sure what I meant by that. 

    @unpack assembly, system, structural_damping, linear, pfunc, distributed_loads, prescribed_conditions = mesh

    t0 = tvec[1] #todo: Why not just pass in t0? Do I need other values? -> I use it later when I'm doing a psuedo step at the beginning. I don't know why I, because I redo that step later down the line. But I think that I was planning on putting that step in the initial condition. 


    na = length(blade.r)

    if verbose
        println("Calculating initial condition...")
    end

    initial_condition_checks(gxflag)


    ### Initialize BEM solution 
    for j = 1:na
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth0)

        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        phi[j] = ccout.phi
        alpha[j] = ccout.alpha
        W[j] = ccout.W
    end

    #Note: I don't use take_aero_step because these functions are different. 
    dsmodel_initial_condition!(xds, phi, W, mesh, blade, rotor.turbine, t0, pitch) #

    
    extract_ds_loads!(blade.airfoils, xds, mesh.xds_idxs, phi, mesh.y_ds, cx, cy, cm)


    dimensionalize!(fx, fy, mx, cx, cy, cm, blade, env, W) 

    
    update_forces!(distributed_loads, fx, fy, mx, blade, assembly)


    ### GXBeam initial solution 
    if isnothing(prepp)
        p = nothing
    else
        prepp(p, fx, fy, mx)
    end

    y0, dy0, gxhistory, _ = gxbeam_initial_conditions!(env, system, assembly, prescribed_conditions, distributed_loads, t0, azimuth0, g, structural_damping, linear, gxflag, pfunc, p) 

    ngx = length(system.x)
    gxstates[1:ngx] = y0
    gxstates[ngx+1:end] = dy0

    # gxstates .= vcat(y0, dy0)
# 

    ### Update mesh transfer variables
    update_mesh!(blade, mesh, assembly, gxhistory, env, t0, na)


    
    ### Trying to take the first step then set that as the first value. 
    # azimuth = env.RS(t)*dt + azimuth0
    
    if isa(xds[1], ReverseDiff.TrackedReal)
        # ns = DS.numberofstates_total(airfoils)
        # xds_old = Array{inittype, 1}(undef, ns) 
        # xds_old = deepcopy(xds)
        inittype = eltype(xds)
        ns = DS.numberofstates_total(blade.airfoils)
        xds_old = Array{inittype, 1}(undef, ns)

        #Todo. How do I copy the derivative? (Because I want to preserve the pointers to xds.)
        for i in eachindex(xds) 
            xds_old[i] = xds[i]
        end
        
    elseif isa(xds[1], ForwardDiff.Dual)
        # xds_old = deepcopy(xds) #Derivatives are getting nuked
        inittype = eltype(xds)
        ns = DS.numberofstates_total(blade.airfoils)
        xds_old = Array{inittype, 1}(undef, ns)

        #Todo. How do I copy the derivative? (Because I want to preserve the pointers to xds.)
        for i in eachindex(xds) 
            xds_old[i] = xds[i]
        end

        # @show xds[1].partials
        # @show xds_old[1].partials
    else
        xds_old = deepcopy(xds)
    end

    dt = tvec[2] - tvec[1]
    take_aero_step!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, xds_old, azimuth0, t0, dt, pitch, mesh, rotor, blade, env; solver, pfunc, prepp, p)
    #todo: I think there is some overall problem in how the initial states are handled. Rather, maybe not a problem, but a difference. 
    #todo: Why am I taking an aero step after I've done the initial initialization and the GX step (I must have a reason. I'm getting pretty good results.) -> I was trying to recreate the initial condition like what OpenFAST does. That makes sense with what I was doing below. As I stick the second solution in the initial index, then the second condition will get overwritten. 

    # for j = 1:na
    #     aerostates.W[1,j] = aerostates.W[2,j]
    #     aerostates.phi[1,j] = aerostates.phi[2,j]
    #     aerostates.xds[1,j] = aerostates.xds[2,j]
        
    #     aerostates.cx[1,j] = aerostates.cx[2,j]
    #     aerostates.cy[1,j] = aerostates.cy[2,j]
    #     aerostates.cm[1,j] = aerostates.cm[2,j]

    #     aerostates.fx[1,j] = aerostates.fx[2,j]
    #     aerostates.fy[1,j] = aerostates.fy[2,j]
    #     aerostates.mx[1,j] = aerostates.mx[2,j]
    # end

    return gxhistory
end

"""
    take_step!()

Take a step based on the previous states. 
"""
function take_step!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, gxstates, t, tprev, azimuth0, xds_old, gxstates_old, blade, mesh, rotor::Rotor, env::Environment, tvec, i, pitch; verbose::Bool=false, speakiter::Int=100, g=9.81, runtimeflag::Bool=false, runtimeiter::Int=speakiter, runtime = (aerostates, gxstates, gxhistory, i)->nothing, solver::Solver=RK4(), prepp=nothing, p=nothing)

    @unpack assembly, system, prescribed_conditions, distributed_loads, point_masses, linear_velocity, xpfunc, pfunc, structural_damping, two_dimensional = mesh

    na = length(blade.r)

    dt = t - tprev
    # if i == 2
    #     dtprev = 0.
    # else
    #     dtprev = tvec[i-1] - tvec[i-2]
    # end
    
    
    if dt<0
        error("Time step is negative")
    end

    #update azimuthal position
    azimuth = env.RS(t)*dt + azimuth0 #Euler step for azimuthal position. 
    #TODO: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

    if azimuth<azimuth0
        @warn("Blade moved backwards")
    end

    #Note: Coupling from structural velocities doesn't appear to kick in until the third time step. 
    take_aero_step!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, xds_old, azimuth, t, dt, pitch, mesh, rotor, blade, env; solver, pfunc, prepp, p)
    
    
    ### Update GXBeam loads 
    update_forces!(distributed_loads, fx, fy, mx, blade, assembly) 

    Omega = SVector(0.0, 0.0, -env.RS(t)) 
    # gravity = SVector(-g*cos(aerostates.azimuth[i-1]), -g*sin(aerostates.azimuth[i-1]), 0.0) #TODO: I need to include tilt, and precone here. 

    # @show typeof(aerostates.azimuth)

    if isa(azimuth0, ForwardDiff.Dual)||isa(azimuth0, ReverseDiff.TrackedReal) #todo: Does azimuth need to be a state?
        a0 = azimuth.value
        a1 = azimuth0.value
    else
        a0 = azimuth
        a1 = azimuth0
    end

    # if isnothing(p)
    #     a0 = azimuth
    #     a1 = azimuth0
    # else
    #     a0 = azimuth.value
    #     a1 = azimuth0.value
    # end

    #TODO: Figure out if I need to use this, or if I can just calculate it at a single point. 
    gravity = (tee) -> SVector(-g*cos((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), -g*sin((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), 0.0) ##Todo t = tvec[i].... So this be way wrong. Oh... this is an inline function so the solver can linearly interpolate the gravity vector across time. But... I think the time domain analysis only analyzes at the given time steps... which means that this function doesn't get called really... I dunno. 

    # @show typeof(gravity2)

    #Note: Taylor applies the gravitational load by C'*mass*C*gvec

    ### Solve GXBeam for time step 
    #TODO: This function is taking a lot of time. -> I might be able to save time by branching his code and writing another function, but most of the time is spent in nlsolve. I think all of the time spent is just time solving, not really inside of Taylor's code, but of course, if I make his code faster, then I make the solve faster. 
    # @show eltype(system)

    if isnothing(prepp)  
        p = nothing
        # println("Got here")
    else
        prepp(p, fx, fy, mx)
    end
    

    ngx = length(system.x)

    xp = view(gxstates_old, 1:ngx)
    dxp = view(gxstates_old, ngx+1:2ngx)

    # @show xp
    # @show dxp

    # println("Running GXBeam...")
    xi, dxi = GXBeam.take_step(xp, dxp, system, assembly, t, tprev, prescribed_conditions, distributed_loads, point_masses, gravity, linear_velocity, Omega, xpfunc, pfunc, p, structural_damping, two_dimensional) 

    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(xi, p, t)
    pcond = get(parameters, :prescribed_conditions, prescribed_conditions)
    pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)

    gxstates[1:ngx] = xi
    gxstates[ngx+1:end] = dxi
    # vcat(xi, dxi) #TODO: There is probably a way to make it so I don't split and recombine the states and state rates quite so much. 
    gxhistory_new = GXBeam.AssemblyState(dxi, xi, system, assembly; prescribed_conditions=pcond)
     

    #TODO: Can I save memory by directly allocating to the gxstates vector? -> I can probably save allocations by augmenting the time_domain_analysis!() function to already have the results allocated. -> This would require a significant restructure of GXBeam.


    ### Update aero inputs from structures.
    update_mesh!(blade, mesh, assembly, gxhistory_new, env, t, na)



    if verbose & (mod(i-1, speakiter)==0) #TODO: remove the dependence on i (move to just a verbose and runtime flag)
        println("")
        println("Simulation time: ", t)
    end



    if runtimeflag & (mod(i-1, runtimeiter)==0) #TODO: Turn this into a runtime function. 
        runtime(aerostates, gxstates_new, gxhistory_new, i) #TODO: Update this to work with the actual aerostates. 
    end

    return azimuth, gxhistory_new
end


function simulate(rotor::Rotors.Rotor, blade::Blade, env::Environment, assembly::GXBeam.Assembly, tvec; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, warnings::Bool=true, azimuth0=0.0, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter)

    nt = length(tvec)

    aerostates, gxstates, mesh = initialize(blade, assembly, tvec; verbose)

    
    system = initial_condition(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch; verbose)

    for i = 2:nt
        system = take_step(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch; verbose, speakiter, plotiter, plotbool, structural_damping, linear, g, solver)
    end

    return aerostates, gxstates
end


function simulate!(rotor::Rotors.Rotor, env::Environment, tvec, aerostates, gxstates, gxhistory, blade, mesh; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, warnings::Bool=true, azimuth0=0.0, g=9.81, runtimeflag::Bool=false, runtimeiter::Int=speakiter, runtime=(aerostates, gxstates, gxhistory, i)->nothing, gxflag=nothing, prepp=nothing, p=nothing) 
    #TODO: Move g to the environment struct. And be evaluated with time (maybe)

    nt = length(tvec)

    #Todo: I need extract the states from aerostates. 

    @unpack azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds = aerostates

    phi0 = view(phi, 1, :)
    alpha0 = view(alpha, 1, :)
    W0 = view(W, 1, :)
    Cx0 = view(Cx, 1, :)
    Cy0 = view(Cy, 1, :)
    Cm0 = view(Cm, 1, :)
    Fx0 = view(Fx, 1, :)
    Fy0 = view(Fy, 1, :)
    Mx0 = view(Mx, 1, :)
    xds0 = view(xds, 1, :)
    gxstates0 = view(gxstates, 1, :)

    gxhistory[1] = initial_condition!(phi0, alpha0, W0, xds0, Cx0, Cy0, Cm0, Fx0, Fy0, Mx0, gxstates0, azimuth0, rotor, env, blade, mesh, tvec, pitch; g, gxflag, verbose, solver, prepp, p)

    # println("Finished initial conditions")
    # @show gxstates[1,:]

    for i = 2:nt
        if verbose
            println("Running step $i ...")
        end

        t = tvec[i]
        tprev = tvec[i-1]

        phii = view(phi, i, :)
        alphai = view(alpha, i, :)
        Wi = view(W, i, :)
        Cxi = view(Cx, i, :)
        Cyi = view(Cy, i, :)
        Cmi = view(Cm, i, :)
        Fxi = view(Fx, i, :)
        Fyi = view(Fy, i, :)
        Mxi = view(Mx, i, :)
        xdsi = view(xds, i, :)
        gxstatesi = view(gxstates, i, :)


        azimuth_old = azimuth[i]
        xds_old = view(xds, i-1, :)
        gxstates_old = view(gxstates, i-1, :)

        azimuth[i], gxhistory[i] = take_step!(phii, alphai, Wi, xdsi, Cxi, Cyi, Cmi, Fxi, Fyi, Mxi, gxstatesi, t, tprev, azimuth_old, xds_old, gxstates_old, blade, mesh, rotor, env, tvec, i, pitch; verbose, speakiter, g, runtimeflag, runtimeiter, runtime, solver, prepp, p)
    end

    return aerostates, gxstates, gxhistory
end

#TODO: Make a memory efficient take_step!() function that only saves certain time indices. 


"""
initialize_sim()

Initialize the data structures for a simulation. Specifically for the run_sim set of functions. 
"""
function initialize_sim(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false, p=nothing, pfunc = (p,t) -> (;), xpfunc=nothing, structural_damping::Bool=true, linear::Bool=false)
    if verbose
        println("Rotors.jl initializing simulation...")
    end

    #Created testing file
    #TODO: Create tests

    #TODO: Look into whether or not I really need this function. If I can just use the initialize function... then why would I need duplicate functions? 

    

    # if warnings
    #     checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    # end

    
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types.

    ### Initialization information
    na = length(blade.rR)
    nt = length(tvec)

    t0 = first(tvec)

    inittype = find_inittype(blade.c[1], blade.twist[1])



    ### ----- Prepare data storage for aerodynamic models ----- ###

    azimuth = Array{inittype}(undef, nt)
    # azimuth = Array{Float64}(undef, nt) #Note: I don't know of a situation that the azimuth would be a dual number.
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
    xds, xds_idxs, y_ds, p_ds = initialize_ds_model(blade, nt, inittype)  

    # Store everything in the aerostates 
    # aerostates = AeroStates(azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)
    aerostates = (;azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)

    

    ### ----- Allocate the GXBeam Data ----- ###

    system = GXBeam.DynamicSystem(assembly)
    # @show typeof(system)

    gxhistory = Array{GXBeam.AssemblyState{inittype, 
        Vector{GXBeam.PointState{inittype}},
        Vector{GXBeam.ElementState{inittype}}}}(undef, nt)

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{Float64}}()

    # Allocate the prescribed conditions 
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)) # root section is fixed

    
    

    # The currently un-used GXBeam constants
    point_masses=Dict{Int,PointMass{Float64}}()
    linear_velocity=(@SVector zeros(3))
    angular_velocity=(@SVector zeros(3))
    # structural_damping=true #todo: I might add this in later. 
    two_dimensional=false
    xpfunc = nothing
    



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 


    delta = Vector{SVector{3, inittype}}(undef, na) #todo: What is this used for? 
    def_theta = Vector{SVector{3, inittype}}(undef, na) #todo. What is this used for
    #The structural velocities interpolated to the aerodynamic nodes.
    aerov = Vector{SVector{3, inittype}}(undef, na)

    
    mesh = (; interpolationpoints, delta, def_theta, aerov, xcc, xds_idxs, y_ds, p_ds,
                assembly, system, prescribed_conditions, distributed_loads,
                point_masses, linear_velocity, angular_velocity,
                xpfunc, pfunc, two_dimensional, structural_damping, linear)

    return aerostates, gxhistory, mesh
end




"""
run_sim!()

The pre-allocated version of run_sim(). 


"""
function run_sim!(rotor::Rotors.Rotor, blade, mesh, env::Environment, tvec, aerostates, gxhistory; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, g=9.81, runtimeflag::Bool=false, runtimeiter::Int=speakiter, runtime = (aerostates, gxhistory, i) ->nothing, gxflag=nothing, prepp=nothing, p=nothing, azimuth0=0.0)


    ### unpack the data structures. 
    @unpack assembly, system, prescribed_conditions, distributed_loads, point_masses, linear_velocity, xpfunc, pfunc, structural_damping, two_dimensional, linear = mesh

    @unpack azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds = aerostates

    na = length(blade.r)
    nt = length(tvec)



    ### Initial Condition analysis. 
    t0 = tvec[1] 

    phi0 = view(phi, 1, :)
    alpha0 = view(alpha, 1, :)
    W0 = view(W, 1, :)
    cx0 = view(Cx, 1, :)
    cy0 = view(Cy, 1, :)
    cm0 = view(Cm, 1, :)
    fx0 = view(Fx, 1, :)
    fy0 = view(Fy, 1, :)
    mx0 = view(Mx, 1, :)
    xds0 = view(xds, 1, :)

    if verbose
        println("Calculating initial condition...")
    end

    initial_condition_checks(gxflag)


    # println("Initializing...")
    ### Initialize BEM solution 
    for j = 1:na
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth0)


        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)
        # ccout = solve_BEM!(rotor, blade, env, 0.0, j, Vx, Vy, pitch, mesh.xcc; newbounds=false)

        phi0[j] = ccout.phi
        alpha0[j] = ccout.alpha
        W0[j] = ccout.W
    end

    #Note: I don't use take_aero_step because these functions are different. 
    dsmodel_initial_condition!(xds0, phi0, W0, mesh, blade, rotor.turbine, t0, pitch)

    
    extract_ds_loads!(blade.airfoils, xds0, mesh.xds_idxs, phi0, mesh.y_ds, mesh.p_ds, cx0, cy0, cm0)


    dimensionalize!(fx0, fy0, mx0, cx0, cy0, cm0, blade::Blade, env::Environment, W0) 

    
    update_forces!(distributed_loads, fx0, fy0, mx0, blade, assembly)


    ### GXBeam initial solution 
    if isnothing(prepp)
        p = nothing
    else
        prepp(p, fx0, fy0, mx0)
    end

    Omega0 = SVector(0.0, 0.0, -env.RS(t0))
    gravity0 = SVector(-g*cos(azimuth0), -g*sin(azimuth0), 0.0)

    # #todo: Only works with initial response (not steady state or spinning solution. )
    # system, gxhistory[1], converged = GXBeam.initial_condition_analysis!(system, assembly, t0; prescribed_conditions, distributed_loads, angular_velocity=Omega0, gravity=gravity0, steady_state=false, structural_damping, linear, pfunc, p, show_trace=false)


    system, gxstate, constants, paug, xgx, converged = GXBeam.initialize_system!(system, assembly, tvec; prescribed_conditions, distributed_loads, gravity=gravity0, angular_velocity=Omega0, structural_damping, reset_state=true, pfunc, p) #todo: This has extra allocations that I don't need.
    


    gxhistory[1] = gxstate[1]

    if !converged
        @warn("The initial condition structural analysis failed to converge.")
        return 
    end



    ### Update mesh transfer variables
    update_mesh!(blade, mesh, assembly, gxhistory[1], env, t0, na)


    
    ### Take the first step then set that as the first value. -> This is what OpenFAST does. 
    # azimuth = env.RS(t)*dt + azimuth0
    
    # if isa(xds0[1], ReverseDiff.TrackedReal) #todo. Move this into a function. 
   
    #     inittype = eltype(xds0)
    #     ns = DS.numberofstates_total(blade.airfoils)
    #     xds_old = Array{inittype, 1}(undef, ns)

    #     #Todo. How do I copy the derivative? (Because I want to preserve the pointers to xds.)
    #     for i in eachindex(xds0) 
    #         xds_old[i] = xds0[i]
    #     end
        
    # elseif isa(xds0[1], ForwardDiff.Dual)
    #     # xds_old = deepcopy(xds) #Derivatives are getting nuked
    #     inittype = eltype(xds0)
    #     ns = DS.numberofstates_total(blade.airfoils)
    #     xds_old = Array{inittype, 1}(undef, ns)

    #     #Todo. How do I copy the derivative? (Because I want to preserve the pointers to xds.)
    #     for i in eachindex(xds0) 
    #         xds_old[i] = xds0[i]
    #     end

    #     # @show xds[1].partials
    #     # @show xds_old[1].partials
    # else
    #     xds_old = deepcopy(xds0)
    # end
    xds_old = dualcopy(xds0) #todo: I feel like there is a better way to do this. 

    dt = tvec[2] - tvec[1] #Note: this passing in phi0 and replacing phi0 might break things... 
    # take_aero_step!(phi0, alpha0, W0, xds0, cx0, cy0, cm0, fx0, fy0, mx0, phi0, xds_old, azimuth0, t0, dt, pitch, mesh, rotor, blade, env; solver, pfunc, prepp, p)
    take_aero_step!(phi0, alpha0, W0, xds0, cx0, cy0, cm0, fx0, fy0, mx0, xds_old, azimuth0, t0, dt, pitch, mesh, rotor, blade, env; solver)





    # println("Beginning time loop...")
    for i in 2:nt
        # println("i = $i")

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


"""
run_sim()

Simulate the physical response of a rotor blade. 
"""
function run_sim()
end

