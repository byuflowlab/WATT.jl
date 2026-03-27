


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
    initial_condition_checks(gxflag)

A function to double check that the inputs are compatible with the solution method.  
"""
function initial_condition_checks(gxflag)

    if !in(gxflag, [nothing, :steady, :spinning])
        error("The flag $gxflag isn't offered for GXBeam initialization.")
    end
end






"""
initialize_sim()

Initialize the data structures for a simulation. Specifically for the run_sim set of functions. 
"""
function initialize_sim(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false, p=nothing, pfunc = (p,t) -> (;), xpfunc=nothing, structural_damping::Bool=true, linear::Bool=false)
    if verbose
        println("WATT.jl initializing simulation...")
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
function run_sim!(rotor::Rotor, blade, mesh, env::Environment, tvec, aerostates, gxhistory; pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, g=9.81, runtimeflag::Bool=false, runtimeiter::Int=speakiter, runtime = (aerostates, gxhistory, i) ->nothing, gxflag=nothing, prepp=nothing, p=nothing, azimuth0=0.0)


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


    ### Initialize BEM solution 
    for j = 1:na
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth0)

        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

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

    system, gxstate, constants, paug, xgx, converged = GXBeam.initialize_system!(system, assembly, tvec; prescribed_conditions, distributed_loads, gravity=gravity0, angular_velocity=Omega0, structural_damping, reset_state=true, pfunc, p) #todo: This has extra allocations that I don't need.


    gxhistory[1] = gxstate[1]

    if !converged
        @warn("The initial condition structural analysis failed to converge.")
        return 
    end



    ### Update mesh transfer variables
    update_mesh!(blade, mesh, assembly, gxhistory[1], env, t0, na)


    
    ### Take the first step then set that as the first value. -> This is what OpenFAST does. 
    xds_old = dualcopy(xds0) #todo: I feel like there is a better way to do this. 

    dt = tvec[2] - tvec[1] #Note: this passing in phi0 and replacing phi0 might break things... 
    take_aero_step!(phi0, alpha0, W0, xds0, cx0, cy0, cm0, fx0, fy0, mx0, xds_old, azimuth0, t0, dt, pitch, mesh, rotor, blade, env; solver)

    azimuth[1] = azimuth0


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
        take_aero_step!(phi_i, alpha_i, W_i, xds_i, cx_i, cy_i, cm_i, fx_i, fy_i, mx_i, xds_im1, azimuth[i], t, dt, pitch, mesh, rotor, blade, env; solver)



        ### Update GXBeam loads 
        update_forces!(distributed_loads, fx_i, fy_i, mx_i, blade, assembly) 

        Omega = SVector(0.0, 0.0, -env.RS(t)) 

        if isa(azimuth[i], ForwardDiff.Dual)||isa(azimuth[i], ReverseDiff.TrackedReal) #todo: Does azimuth need to be a state?
            a0 = azimuth[i].value
            a1 = azimuth[i-1].value
        else
            a0 = azimuth[i]
            a1 = azimuth[i-1]
        end

        #TODO: Figure out if I need to use this, or if I can just calculate it at a single point. 
        gravity = (tee) -> SVector(-g*cos((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), -g*sin((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), 0.0) ##Todo t = tvec[i].... So this be way wrong. Oh... this is an inline function so the solver can linearly interpolate the gravity vector across time. But... I think the time domain analysis only analyzes at the given time steps... which means that this function doesn't get called really... I dunno. -> TODO: This would only make a difference if the solver used intermediate steps (using DifferentialEquations to solve GXBeam.) -> I should get rid of this. 

        #Note: Taylor applies the gravitational load by C'*mass*C*gvec

        ### Solve GXBeam for time step 
        #todo: This function is taking a lot of time. -> I might be able to save time by branching his code and writing another function, but most of the time is spent in nlsolve. I think all of the time spent is just time solving, not really inside of Taylor's code, but of course, if I make his code faster, then I make the solve faster. 

        if isnothing(prepp)  
            p = nothing
        else
            prepp(p, fx_i, fy_i, mx_i)
        end 

        

    
        system, gxhistory[i], constants, paug, xgx, convergedi = GXBeam.step_system!(system, paug, xgx, constants, gxhistory[i-1], assembly, tvec, i; prescribed_conditions, distributed_loads, structural_damping, gravity, angular_velocity=Omega, pfunc=pfunc, p=p) 


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

Simulate the physical response of a rotor blade. Preallocate the data structures using initialize_sim() and then run the time stepping loop using run_sim!().
"""
function run_sim()
end

