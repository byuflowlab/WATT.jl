


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

"""
    find_inittype(vars...) -> Type

Pick the element type that pre-allocated simulation buffers should carry.

Scans `vars` in order and returns the type of the first `ForwardDiff.Dual` or
`ReverseDiff.TrackedReal` it finds; if none is present, returns `eltype(vars)`.

**Arguments**
- `vars...`: Representative scalars from the problem definition — typically
  `blade.c[1]` and `blade.twist[1]`.

**Returns**
- `Type`: The element type to allocate with.

**Notes**
Internal. This is what makes the `initialize_*` functions AD-transparent: pass a
dual-valued chord or twist and every downstream buffer is allocated wide enough
to carry duals, so no retyping is needed mid-simulation. Used by
[`initialize_sim`](@ref), [`initialize_aero`](@ref), and
[`initialize_static`](@ref).
"""
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
    initialize_sim(blade, assembly, tvec; structural_damping=true, linear=false, kwargs...) -> aerostates, gxhistory, mesh

Pre-allocate every buffer the coupled aero-structural transient solver needs.
The returned tuple is the canonical input to [`run_sim!`](@ref).

**Arguments**
- `blade::Blade`
- `assembly::GXBeam.Assembly`
- `tvec::AbstractVector`: Time vector; `length(tvec)` sets the history depth.

**Keyword Arguments**
- `structural_damping::Bool = true`: Passed through to GXBeam.
- `linear::Bool = false`: Use linear structural response.
- `pfunc`, `xpfunc`: Sensitivity parameter callbacks (advanced). Xpfunc currently unused. 
- `p::Vector{<:Real} = nothing`: The sensitivity parameter vector to pass to the pfunc and xpfunc functions (advanced).
- `verbose::Bool = false`: Print progress.

**Returns**
- `aerostates::AeroStates`: Time-indexed aero state history.
- `gxhistory::Vector{GXBeam.AssemblyState}`: Per-step structural state.
- `mesh::SimMesh`: Coupling buffers reused at every time step.

**Notes**
The element type of the allocated buffers is inferred from `blade.c[1]` and
`blade.twist[1]` via [`find_inittype`](@ref), so passing `ForwardDiff.Dual`
chord/twist propagates duals through every buffer.
"""
function initialize_sim(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false, p=nothing, pfunc = (p,t) -> (;), xpfunc=nothing, structural_damping::Bool=true, linear::Bool=false)
    if verbose
        println("WATT.jl initializing simulation...")
    end

    #TODO: Why don't I store p and prepp in the mesh? -> It's not that big of a deal to pass them around... so might not be worth the effort. 

    # if warnings
    #     checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    # end #TODO: Why did I comment out these checks? -> What if I just have the user check the inputs? Like check for warnings could be a function that the user can call if they want, but I don't want to call it by default.
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types.

    ### Initialization information
    na = length(blade.rR)
    nt = length(tvec)

    t0 = first(tvec)

    inittype = find_inittype(blade.c[1], blade.twist[1])


    ### ----- Prepare data storage for aerodynamic models ----- ###

    # Initialize DS state vector first so we know its width (ns).
    xds, xds_idxs, y_ds, p_ds = initialize_ds_model(blade, nt, inittype)
    ns = size(xds, 2)

    aerostates = AeroStates{inittype}(undef, nt, na, ns)
    # The DS state history lives in `aerostates.xds`; share the buffer that
    # `initialize_ds_model` already populated so call sites that pass `xds`
    # around continue to work.
    copyto!(aerostates.xds, xds) #Copying the initialized xds values into the aerostates buffer. This is a little bit redundant, but it allows me to keep the xds buffer that I initialized for the DS model, and also have it be part of the aerostates struct. -> Could update initialize_ds_model to accept the aerostates buffer and write directly into that, but this is probably fine for now.
    xds = aerostates.xds #Aliasing the buffer... but it 

    # CCBlade scratch vector
    xcc = Vector{inittype}(undef, 11)

    

    ### ----- Allocate the GXBeam Data ----- ###

    system = GXBeam.DynamicSystem(assembly)

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
    two_dimensional=false
    xpfunc = nothing
    



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 


    delta = Vector{SVector{3, inittype}}(undef, na) #Linear deflection of the aerodynamic nodes relative to the reference configuration.
    def_theta = Vector{SVector{3, inittype}}(undef, na) #Angular deflection of the aerodynamic nodes relative to the reference configuration.
    #The structural velocities interpolated to the aerodynamic nodes.
    aerov = Vector{SVector{3, inittype}}(undef, na)

    
    mesh = SimMesh(interpolationpoints, delta, def_theta, aerov, xcc,
                xds_idxs, y_ds, p_ds,
                assembly, system, prescribed_conditions, distributed_loads,
                point_masses, linear_velocity, angular_velocity,
                xpfunc, pfunc, two_dimensional, structural_damping, linear)

    return aerostates, gxhistory, mesh
end




"""
    run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; pitch=0.0, solver=RK4(), kwargs...)

Pre-allocated transient aero-structural simulation. Mutates `aerostates`,
`gxhistory`, and the mutable fields of `mesh` in place. The first step is
treated as an initial condition (BEM + DS init + GXBeam `initialize_system!`),
the remaining `length(tvec)-1` steps advance with `take_aero_step!` and
`GXBeam.step_system!`.

**Arguments**
- `rotor::Rotor`, `blade::Blade`, `mesh::SimMesh`, `env::Environment`
- `tvec::AbstractVector`: Time vector.
- `aerostates::AeroStates`, `gxhistory`: From [`initialize_sim`](@ref).

**Keyword Arguments**
- `pitch::Real = 0.0`: Blade pitch (rad).
- `solver::Solver = RK4()`: DS state integrator.
- `g::Real = 9.81`: Gravity magnitude.
- `azimuth0::Real = 0.0`: Initial azimuth.
- `gxflag::Union{Nothing,Symbol} = nothing`: A flag to initialize the structures in a particular way. Options are `nothing`, `:steady`, and `:spinning`. Currently does nothing. 
- `verbose::Bool = false`, `speakiter::Int = 100`: Progress printing.
- `runtimeflag::Bool = false`, `runtimeiter::Int`, `runtime`: Custom per-step callback.
- `prepp::Function = nothing`: A function to update the loads within the parameter vector. 
- `p`::Vector{<:Real} = nothing: The sensitivity parameter vector to pass to the pfunc and prepp functions.


**Notes**
AD through this function relies on the parameter prep function and p. See examples for how to use this for AD.
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

        ccout = solve_BEMT(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

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


    #Note: Coupling from structural velocities doesn't appear to kick in until the third time step.
    #Note: the single-step logic lives in `step_solution!`; the loop below just threads the Newmark
    #carriers (xgx/paug/constants), the aero history views, and the structural history.
    for i in 2:nt

        ### Unpack single-step aero output views (writes land directly in the history arrays).
        aero = (; phi   = view(phi,   i, :),
                  alpha = view(alpha, i, :),
                  W     = view(W,     i, :),
                  Cx    = view(Cx,    i, :),
                  Cy    = view(Cy,    i, :),
                  Cm    = view(Cm,    i, :),
                  Fx    = view(Fx,    i, :),
                  Fy    = view(Fy,    i, :),
                  Mx    = view(Mx,    i, :))
        xds_i   = view(xds, i, :)
        xds_im1 = view(xds, i-1, :)

        ### Advance one coupled Newmark step (see step_solution!).
        gxhistory[i], xgx, paug, constants, azimuth[i], convergedi = step_solution!(
            gxhistory[i-1], xgx, paug, constants, aero, xds_i, xds_im1, azimuth[i-1],
            mesh, rotor, blade, env, tvec, i;
            pitch, solver, g, prepp, p)

        if !convergedi
            @warn("GXBeam failed to converge on the $i th time step.")
            break
        end

        t = tvec[i]
        if verbose & (mod(i-1, speakiter)==0) #todo: remove the dependence on i (move to just a verbose and runtime flag)
            println("")
            println("Simulation time: ", t)
        end

        if runtimeflag & (mod(i-1, runtimeiter)==0)
            runtime(aerostates, gxhistory[i], i)
        end
    end
end


"""
    run_sim(rotor, blade, assembly, env, tvec; kwargs...) -> aerostates, gxhistory, mesh

Allocating wrapper for [`run_sim!`](@ref). Equivalent to calling
[`initialize_sim`](@ref) and then `run_sim!`, returning the populated
[`AeroStates`](@ref), GXBeam history, and [`SimMesh`](@ref).

**Arguments**
- `rotor::Rotor`
- `blade::Blade`
- `assembly::GXBeam.Assembly`
- `env::Environment`
- `tvec::AbstractVector`: Time vector.

**Keyword Arguments**
Forwarded to `run_sim!`. See its docstring for the full list.

**Notes**
For repeated solves (e.g. in an optimization outer loop), prefer the
in-place pair `initialize_sim` + `run_sim!` to reuse the buffers.
"""
function run_sim(rotor::Rotor, blade::Blade, assembly::GXBeam.Assembly, env::Environment, tvec; kwargs...)
    aerostates, gxhistory, mesh = initialize_sim(blade, assembly, tvec)
    run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; kwargs...)
    return aerostates, gxhistory, mesh
end


# Strip a ForwardDiff/ReverseDiff dual to its primal value; passthrough otherwise.
# Used to build the (non-dual) gravity-interpolation angles inside a coupled step.
_azimuth_value(x) = (isa(x, ForwardDiff.Dual) || isa(x, ReverseDiff.TrackedReal)) ? x.value : x


"""
    step_solution!(state_prev, xgx, paug, constants, aero, xds_new, xds_old, azimuth_prev,
                   mesh, rotor, blade, env, tvec, i;
                   pitch=0.0, solver=RK4(), g=9.81, prepp=nothing, p=nothing)
        -> (state_new, xgx, paug, constants, azimuth, converged)

Advance the coupled aero-structural solution **exactly one Newmark-β time step**, from the state at
`tvec[i-1]` to `tvec[i]`. This is the per-step primitive factored out of [`run_sim!`](@ref)'s time
loop (`run_sim!` calls it), and the building block for windowed frozen-start sensitivity sweeps via
[`initialize_from_state`](@ref) / [`run_from_state!`](@ref). It makes **no rest-IC assumption**: the
starting state is an input.

The sequence mirrors `run_sim!`'s loop body:
1. Integrate the azimuth (`azimuth = env.RS(t)*dt + azimuth_prev`).
2. Aero step (`take_aero_step!`) from `xds_old` → `xds_new` + sectional loads (written into `aero`).
3. `update_forces!` to push aero loads into the GXBeam distributed loads.
4. Structural step via `GXBeam.step_system!` → new flat state `xgx`, updated `paug`/`constants`, and
   a new `AssemblyState`.
5. `update_mesh!` (only when converged) to feed the new deflections/velocities back to the aero inputs.

**Arguments**
- `state_prev::GXBeam.AssemblyState`: structural state at `tvec[i-1]` (carries the Newmark rates).
- `xgx`, `paug`, `constants`: GXBeam Newmark solver carriers (flat state, augmented rate-init/parameter
  vector, and step constants) — as produced by [`initialize_from_state`](@ref) or a prior step.
- `aero`: single-step aero scratch exposing `.phi/.alpha/.W/.Cx/.Cy/.Cm/.Fx/.Fy/.Mx` (a
  [`StaticAeroStates`](@ref) or a `NamedTuple` of length-`na` vectors; written in place).
- `xds_new`, `xds_old`: DS state buffers (written / read).
- `azimuth_prev::Real`: azimuth at `tvec[i-1]`.
- `mesh::AbstractSimMesh`, `rotor::Rotor`, `blade::Blade`, `env::Environment`
- `tvec::AbstractVector`, `i::Int`: the step advances `tvec[i-1] → tvec[i]`; `step_system!` also reads
  `tvec[i-2]` for `i>2` (the Newmark rate recursion), so pass the same time grid you initialized with.

**Keyword Arguments**
- `pitch::Real = 0.0`, `solver::Solver = RK4()`, `g::Real = 9.81`
- `prepp::Function = nothing`, `p = nothing`: sensitivity-parameter prep/vector (see [`run_sim!`](@ref)).

**Returns**
`(state_new, xgx, paug, constants, azimuth, converged)`. `xds_new` is written in place.
"""
function step_solution!(state_prev, xgx, paug, constants, aero, xds_new, xds_old, azimuth_prev,
                        mesh, rotor::Rotor, blade::Blade, env::Environment, tvec, i;
                        pitch=0.0, solver::Solver=RK4(), g=9.81, prepp=nothing, p=nothing)

    @unpack assembly, system, prescribed_conditions, distributed_loads, pfunc, structural_damping = mesh

    na = length(blade.r)

    t = tvec[i]
    tprev = tvec[i-1]
    dt = t - tprev
    if dt < 0
        error("Time step is negative")
    end

    # 1. Azimuth (Euler step for the azimuthal position).
    azimuth = env.RS(t)*dt + azimuth_prev
    if azimuth < azimuth_prev
        @warn("Blade moved backwards")
    end

    # 2. Aero step (BEM + DS) from the previous DS state → new DS state + sectional loads.
    take_aero_step!(aero.phi, aero.alpha, aero.W, xds_new, aero.Cx, aero.Cy, aero.Cm,
                    aero.Fx, aero.Fy, aero.Mx, xds_old, azimuth, t, dt, pitch,
                    mesh, rotor, blade, env; solver)

    # 3. Push aero loads into the GXBeam distributed loads.
    update_forces!(distributed_loads, aero.Fx, aero.Fy, aero.Mx, blade, assembly)

    Omega = SVector(0.0, 0.0, -env.RS(t))

    # Gravity swings with the (linearly interpolated) azimuth over the step; build it from the
    # primal azimuth values so the closure stays cheap. Note: Taylor applies gravity via C'*mass*C*gvec.
    a0 = _azimuth_value(azimuth)
    a1 = _azimuth_value(azimuth_prev)
    gravity = (tee) -> SVector(-g*cos((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)),
                               -g*sin((a0*(t-tee) + a1*(tee-tprev))/(t-tprev)), 0.0)

    # 4. Update the sensitivity parameter vector with the new aero loads, if requested.
    if isnothing(prepp)
        p = nothing
    else
        prepp(p, aero.Fx, aero.Fy, aero.Mx)
    end

    # Structural step (Newmark-β) for this time step.
    _, state_new, constants, paug, xgx, converged = GXBeam.step_system!(
        system, paug, xgx, constants, state_prev, assembly, tvec, i;
        prescribed_conditions, distributed_loads, structural_damping,
        gravity, angular_velocity=Omega, pfunc=pfunc, p=p)

    # 5. Feed the new structural deflections/velocities back to the aero inputs (skip on failure,
    # matching run_sim!'s break-before-update behavior).
    if converged
        update_mesh!(blade, mesh, assembly, state_new, env, t, na)
    end

    return state_new, xgx, paug, constants, azimuth, converged
end


"""
    initialize_from_state(state0, xds0, azimuth0, mesh, blade, env, tvec_window, t0;
                          prescribed_conditions=mesh.prescribed_conditions,
                          distributed_loads=mesh.distributed_loads,
                          gravity=nothing, angular_velocity=nothing,
                          structural_damping=mesh.structural_damping, g=9.81,
                          pfunc=mesh.pfunc, p=nothing)
        -> (state0, xgx, paug, constants, xds0, azimuth0)

Make a caller-supplied `mesh`/`system` consistent with an **injected** structural state so that
[`step_solution!`](@ref) / [`run_from_state!`](@ref) can march from it — **without** any rest-IC
solve. Reuses `GXBeam.initialize_system!(...; initial_state=state0)`, which skips
`initial_condition_analysis!` entirely and only builds the Newmark carriers (`constants`, `paug`,
`xgx`) from the injected state.

The returned tuple is exactly the argument prefix [`run_from_state!`](@ref) expects (splat it).

**Caller constraints**
- `mesh` must be built by [`initialize_sim`](@ref) (allocation-only). For an AD sweep, allocate it for
  the marched eltype (e.g. pass a `Dual`-promoted `blade` to `initialize_sim`) so `update_mesh!` can
  write dual-typed deflections.
- `state0::GXBeam.AssemblyState` must carry the Newmark rates (`udot/θdot/Vdot/Ωdot`), i.e. it must be a
  state produced by a prior `run_sim!`/`run_from_state!` (`gxhistory[s]`), **not** a hand-built
  displacement-only state — the rates seed the first window step.
- **Frozen-start AD:** pass a `Float64` `state0`/`xds0` but a `Dual` `p`; GXBeam promotes `xgx/paug`
  through `pfunc`/`p`, so the injected state enters with **zero partials** (derivative seeded only by `p`).

**Arguments**
- `state0`, `xds0`, `azimuth0`: the snapshot triple (structural state, DS state, azimuth) at step `s`.
- `mesh`, `blade`, `env`: as in [`run_sim!`](@ref).
- `tvec_window`: the window time grid `[t_s, t_{s+1}, …, t_{s+N}]`; `first(tvec_window)` should equal `t0`.
- `t0`: the snapshot time (window start).

**Keyword Arguments**
- `gravity`, `angular_velocity`: defaults consistent with `azimuth0`/`env.RS(t0)`. These are only stored
  as `constants` defaults (each `step_solution!` overrides them per step), so they are largely inert.
- others as in [`run_sim!`](@ref).
"""
function initialize_from_state(state0, xds0, azimuth0, mesh, blade::Blade, env::Environment, tvec_window, t0;
                               prescribed_conditions=mesh.prescribed_conditions,
                               distributed_loads=mesh.distributed_loads,
                               gravity=nothing, angular_velocity=nothing,
                               structural_damping::Bool=mesh.structural_damping, g=9.81,
                               pfunc=mesh.pfunc, p=nothing)

    @unpack assembly, system = mesh
    na = length(blade.r)

    a0 = _azimuth_value(azimuth0)
    if isnothing(angular_velocity)
        angular_velocity = SVector(0.0, 0.0, -env.RS(t0))
    end
    if isnothing(gravity)
        gravity = SVector(-g*cos(a0), -g*sin(a0), 0.0)
    end

    # initial_state provided → initialize_system! skips the rest-IC analysis and just builds the
    # Newmark carriers from the injected state (promoting through p for AD).
    _, _, constants, paug, xgx, _ = GXBeam.initialize_system!(
        system, assembly, tvec_window; initial_state=state0, reset_state=false,
        prescribed_conditions, distributed_loads, gravity, angular_velocity,
        structural_damping, pfunc, p)

    # Populate the coupling buffers from the snapshot so the first take_aero_step! sees deflections
    # consistent with state0 (mirrors run_sim!'s post-IC update_mesh!).
    update_mesh!(blade, mesh, assembly, state0, env, t0, na)

    return state0, xgx, paug, constants, xds0, azimuth0
end


"""
    run_from_state!(state0, xgx, paug, constants, xds0, azimuth0, mesh, rotor, blade, env, tvec_window;
                    pitch=0.0, solver=RK4(), g=9.81, prepp=nothing, p=nothing, out=nothing)
        -> (state_final, xgx, paug, constants, xds, azimuth, out)

March [`step_solution!`](@ref) over a window, starting from the carriers returned by
[`initialize_from_state`](@ref) (splat that tuple into the leading arguments). Single-step aero/DS
buffers are allocated once here (at `eltype(xgx)`, so they follow `Dual` under AD) and reused; the
caller owns any history storage via the `out` hook.

**Arguments**
- `state0, xgx, paug, constants, xds0, azimuth0`: the tuple from [`initialize_from_state`](@ref).
- `mesh`, `rotor`, `blade`, `env`, `tvec_window`: as in [`step_solution!`](@ref). Marches local index
  `i = 2:length(tvec_window)`.

**Keyword Arguments**
- `pitch`, `solver`, `g`, `prepp`, `p`: forwarded to [`step_solution!`](@ref).
- `out`: optional per-step sink. If a `Function`, called as `out(state, i)` after each converged step;
  if an `AbstractVector`, assigns `out[i] = state`. Use it to capture the deflection history (e.g. for a
  `ForwardDiff` Jacobian of `∂u/∂p`).

**Returns**
`(state_final, xgx, paug, constants, xds, azimuth, out)`.
"""
function run_from_state!(state0, xgx, paug, constants, xds0, azimuth0,
                         mesh, rotor::Rotor, blade::Blade, env::Environment, tvec_window;
                         pitch=0.0, solver::Solver=RK4(), g=9.81, prepp=nothing, p=nothing, out=nothing)

    nt = length(tvec_window)
    na = length(blade.r)

    TF = eltype(xgx)

    # Single-step scratch, reused across the window. Promote the DS state to the marched eltype so it
    # enters with zero partials under AD (frozen start).
    aero = StaticAeroStates{TF}(undef, na)
    xds_old = TF.(xds0)
    xds_new = similar(xds_old)

    state = state0
    azimuth = azimuth0

    for i in 2:nt
        state, xgx, paug, constants, azimuth, converged = step_solution!(
            state, xgx, paug, constants, aero, xds_new, xds_old, azimuth,
            mesh, rotor, blade, env, tvec_window, i;
            pitch, solver, g, prepp, p)

        if !converged
            @warn("GXBeam failed to converge on window step $i.")
            break
        end

        if isa(out, Function)
            out(state, i)
        elseif isa(out, AbstractVector)
            out[i] = state
        end

        # Ping-pong the DS buffers: this step's output becomes next step's input.
        xds_old, xds_new = xds_new, xds_old
    end

    return state, xgx, paug, constants, xds_old, azimuth, out
end

