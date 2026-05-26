#=
Surrogate-driven aerostructural transient solver. Mirrors aerostructural.jl,
replacing the GXBeam structural advance with three pure calls into a
user-supplied AbstractStructuralSurrogate: encode_initial, step, decode.
The aero side (BEMT + dynamic stall) is unchanged.
=#


"""
    SurrogateMesh{TF, TIP, TXDS, TYDS, TPDS, TA}

Scratch state for the surrogate-driven transient simulation. Carries the
aero-coupling fields shared with [`SimMesh`](@ref) and [`AeroMesh`](@ref)
plus an `assembly` handle for the structural-point interpolation, but no
GXBeam `DynamicSystem` / `distributed_loads` / `prescribed_conditions`.

**Fields**
- `interpolationpoints::TIP`: Cached structural→aero interpolation handles.
- `delta::Vector{SVector{3,TF}}`: Per-aero-station structural displacement (updated each step).
- `def_theta::Vector{SVector{3,TF}}`: Per-aero-station structural rotation (updated each step).
- `aerov::Vector{SVector{3,TF}}`: Per-aero-station structural velocity (updated each step).
- `xcc::Vector{TF}`: CCBlade scratch vector.
- `xds_idxs::TXDS`, `y_ds::TYDS`, `p_ds::TPDS`: DS model scratch.
- `assembly::TA`: Structural assembly (reference configuration; used by `update_mesh!`).
"""
struct SurrogateMesh{TF, TIP, TXDS, TYDS, TPDS, TA} <: AbstractSimMesh
    interpolationpoints::TIP
    delta::Vector{SVector{3, TF}}
    def_theta::Vector{SVector{3, TF}}
    aerov::Vector{SVector{3, TF}}
    xcc::Vector{TF}
    xds_idxs::TXDS
    y_ds::TYDS
    p_ds::TPDS
    assembly::TA
end


"""
    build_surrogate_loads(fx, fy, mx, assembly, g, azimuth) -> Matrix

Assemble the per-element 6-column load matrix the structural surrogate
was trained on. Element length `L` and element mass density `m` are read
from `assembly.elements[j]`. Pure / allocating so dual types in `fx`,
`fy`, `mx`, or `azimuth` propagate cleanly.

Column convention: `[Fx Fy Fz Mx My Mz]`. `Mx, My, Mz` are zeroed (the
surrogate's training data set the BEMT moment to zero on the load side).

**Arguments**
- `fx, fy, mx::AbstractVector`: Sectional aerodynamic loads at the
  `nelem` aero stations (matched 1:1 with structural elements).
- `assembly::GXBeam.Assembly`: Structural assembly.
- `g::Real`: Gravity magnitude.
- `azimuth::Real`: Current azimuthal position (rad).
"""
function build_surrogate_loads(fx::AbstractVector, fy::AbstractVector, mx::AbstractVector,
                                assembly::GXBeam.Assembly, g, azimuth)
    nelem = length(assembly.elements)
    TF = promote_type(eltype(fx), eltype(fy), eltype(mx), typeof(g), typeof(azimuth))
    f = zeros(TF, nelem, 6)
    ca = cos(azimuth)
    sa = sin(azimuth)
    for j in 1:nelem
        L = assembly.elements[j].L
        m = assembly.elements[j].mass[1, 1] * L
        f[j, 1] = -g * m * ca
        f[j, 2] = fy[j] * L / 2 - g * m * sa
        f[j, 3] = -fx[j] * L / 2
        # f[j, 4:6] left at zero
    end
    return f
end


"""
    initialize_sim_surrogate(blade, assembly, tvec; verbose=false) -> (aerostates, surr_history, mesh)

Pre-allocate every buffer the surrogate-driven transient solver needs.
Mirrors [`initialize_sim`](@ref) but skips the GXBeam allocations and
returns a `SurrogateAssemblyState` history in place of `gxhistory`.

**Arguments**
- `blade::Blade`
- `assembly::GXBeam.Assembly`
- `tvec::AbstractVector`: Time vector.

**Returns**
- `aerostates::AeroStates`
- `surr_history::Vector{SurrogateAssemblyState}`: One entry per time step,
  fully populated after [`run_sim_surrogate!`](@ref).
- `mesh::SurrogateMesh`
"""
function initialize_sim_surrogate(blade::Blade, assembly::GXBeam.Assembly, tvec; verbose::Bool=false)
    if verbose
        println("WATT.jl initializing surrogate simulation...")
    end

    na = length(blade.rR)
    nt = length(tvec)
    np = length(assembly.points)

    inittype = find_inittype(blade.c[1], blade.twist[1])

    ### Aero buffers (identical to initialize_sim)
    xds, xds_idxs, y_ds, p_ds = initialize_ds_model(blade, nt, inittype)
    ns = size(xds, 2)

    aerostates = AeroStates{inittype}(undef, nt, na, ns)
    copyto!(aerostates.xds, xds)

    xcc = Vector{inittype}(undef, 11)

    ### Surrogate history
    surr_history = Vector{SurrogateAssemblyState{inittype}}(undef, nt)
    for i in 1:nt
        surr_history[i] = zero_surrogate_state(inittype, np)
    end

    ### Structural→aero interpolation
    interpolationpoints = create_interpolationpoints(assembly, blade)

    delta     = Vector{SVector{3, inittype}}(undef, na)
    def_theta = Vector{SVector{3, inittype}}(undef, na)
    aerov     = Vector{SVector{3, inittype}}(undef, na)
    for j in 1:na
        delta[j]     = SVector{3, inittype}(0.0, 0.0, 0.0)
        def_theta[j] = SVector{3, inittype}(0.0, 0.0, 0.0)
        aerov[j]     = SVector{3, inittype}(0.0, 0.0, 0.0)
    end

    mesh = SurrogateMesh(interpolationpoints, delta, def_theta, aerov, xcc,
                         xds_idxs, y_ds, p_ds, assembly)

    return aerostates, surr_history, mesh
end


"""
    run_sim_surrogate!(rotor, blade, mesh, env, tvec, aerostates, surr_history, surr;
                      pitch=0.0, solver=RK4(), g=9.81, azimuth0=0.0,
                      verbose=false, speakiter=100,
                      runtimeflag=false, runtimeiter=speakiter,
                      runtime=(aerostates, surr_history, i) -> nothing)

Run the transient aerostructural simulation using `surr` (an
[`AbstractStructuralSurrogate`](@ref)) in place of GXBeam for the
structural advance. The first step is treated as an initial condition
(BEM + DS init + `encode_initial(surr, u0)`); each subsequent step
solves BEM/DS, assembles the surrogate's per-element load matrix via
[`build_surrogate_loads`](@ref), advances the latent state via
`step_latent(surr, z, f)`, and decodes onto the structural-point grid via
`decode(surr, z)`.

**Surrogate interface (pure / AD-safe).** `surr` must implement:
- `z0     = encode_initial(surr, u0::SurrogateAssemblyState)`
- `z_next = step_latent(surr, z, f::AbstractMatrix)` — `f` is `(nelem, 6)`
- `state  = decode(surr, z)::SurrogateAssemblyState`

No method may mutate `surr` or `z`. The latent state is threaded as a
local variable inside the loop, so re-binding is AD-safe.

**Arguments**
- `rotor::Rotor`, `blade::Blade`, `mesh::SurrogateMesh`, `env::Environment`
- `tvec::AbstractVector`
- `aerostates::AeroStates`, `surr_history`: From [`initialize_sim_surrogate`](@ref).
- `surr::AbstractStructuralSurrogate`: Pre-conditioned surrogate (FiLM
  / `x` already baked in by the caller).

**Keyword Arguments**
- `u0_struct::Union{Nothing, SurrogateAssemblyState} = nothing`: Initial
  structural state passed to `encode_initial`. Default `nothing` uses an
  at-rest zero state; pass a precomputed structural equilibrium (e.g.
  the first slot of a GXBeam baseline run, converted via the user's
  state-extraction convention) when the surrogate was trained on data
  whose first frame is non-zero.
- `pitch::Real = 0.0`
- `solver::Solver = RK4()`
- `g::Real = 9.81`
- `azimuth0::Real = 0.0`
- `verbose::Bool = false`, `speakiter::Int = 100`
- `runtimeflag::Bool = false`, `runtimeiter::Int`, `runtime`: per-step callback.
"""
function run_sim_surrogate!(rotor::Rotor, blade::Blade, mesh::SurrogateMesh,
                            env::Environment, tvec,
                            aerostates::AeroStates, surr_history,
                            surr::AbstractStructuralSurrogate;
                            u0_struct=nothing,
                            pitch=0.0, solver::Solver=RK4(), g=9.81, azimuth0=0.0,
                            verbose::Bool=false, speakiter::Int=100,
                            runtimeflag::Bool=false, runtimeiter::Int=speakiter,
                            runtime=(aerostates, surr_history, i) -> nothing)

    assembly = mesh.assembly

    @unpack azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds = aerostates

    na = length(blade.r)
    nt = length(tvec)

    ### --- Initial condition ---
    t0 = tvec[1]

    phi0   = view(phi, 1, :)
    alpha0 = view(alpha, 1, :)
    W0     = view(W, 1, :)
    cx0    = view(Cx, 1, :)
    cy0    = view(Cy, 1, :)
    cm0    = view(Cm, 1, :)
    fx0    = view(Fx, 1, :)
    fy0    = view(Fy, 1, :)
    mx0    = view(Mx, 1, :)
    xds0   = view(xds, 1, :)

    if verbose
        println("Calculating initial surrogate condition...")
    end

    for j = 1:na
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth0)
        ccout = solve_BEMT(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)
        phi0[j]   = ccout.phi
        alpha0[j] = ccout.alpha
        W0[j]     = ccout.W
    end

    dsmodel_initial_condition!(xds0, phi0, W0, mesh, blade, rotor.turbine, t0, pitch)
    extract_ds_loads!(blade.airfoils, xds0, mesh.xds_idxs, phi0, mesh.y_ds, mesh.p_ds, cx0, cy0, cm0)
    dimensionalize!(fx0, fy0, mx0, cx0, cy0, cm0, blade, env, W0)

    azimuth[1] = azimuth0

    ### Encode the IC into the latent space.
    # If the caller didn't supply one, default to the at-rest zero state.
    # Surrogates trained on data whose first frame is a non-zero structural
    # equilibrium (e.g. GXBeam's initialize_system! output) should be given
    # that equilibrium here via the `u0_struct` kwarg.
    np = length(assembly.points)
    inittype = find_inittype(blade.c[1], blade.twist[1])
    u0_ic = u0_struct === nothing ? zero_surrogate_state(inittype, np) : u0_struct
    z = encode_initial(surr, u0_ic) #
    surr_history[1] = decode(surr, z) #Todo: The zero struct doesn't map to a zero surrogate state after decode, at least with what's currently in the example. 

    update_mesh!(blade, mesh, assembly, surr_history[1], env, t0, na)

    ### Warm-up aero step (matches the GXBeam path; lines 286–292 in aerostructural.jl).
    xds_old = dualcopy(xds0)
    dt0 = tvec[2] - tvec[1]
    take_aero_step!(phi0, alpha0, W0, xds0, cx0, cy0, cm0, fx0, fy0, mx0,
                    xds_old, azimuth0, t0, dt0, pitch, mesh, rotor, blade, env; solver)

    ### --- Main loop ---
    for i in 2:nt
        phi_i   = view(phi, i, :)
        alpha_i = view(alpha, i, :)
        W_i     = view(W, i, :)
        cx_i    = view(Cx, i, :)
        cy_i    = view(Cy, i, :)
        cm_i    = view(Cm, i, :)
        fx_i    = view(Fx, i, :)
        fy_i    = view(Fy, i, :)
        mx_i    = view(Mx, i, :)
        xds_i   = view(xds, i, :)
        xds_im1 = view(xds, i-1, :)

        t     = tvec[i]
        tprev = tvec[i-1]
        dt    = t - tprev

        if dt < 0
            error("Time step is negative")
        end

        azimuth[i] = env.RS(t) * dt + azimuth[i-1]

        if azimuth[i] < azimuth[i-1]
            @warn("Blade moved backwards")
        end

        take_aero_step!(phi_i, alpha_i, W_i, xds_i, cx_i, cy_i, cm_i,
                        fx_i, fy_i, mx_i, xds_im1, azimuth[i], t, dt, pitch,
                        mesh, rotor, blade, env; solver)

        ### Build the surrogate force matrix and advance the latent state.
        f_mat = build_surrogate_loads(fx_i, fy_i, mx_i, assembly, g, azimuth[i])
        z = step_latent(surr, z, f_mat)
        surr_history[i] = decode(surr, z)

        update_mesh!(blade, mesh, assembly, surr_history[i], env, t, na)

        if verbose & (mod(i-1, speakiter) == 0)
            println("")
            println("Simulation time: ", t)
        end

        if runtimeflag & (mod(i-1, runtimeiter) == 0)
            runtime(aerostates, surr_history, i)
        end
    end

    return nothing
end


"""
    run_sim_surrogate(rotor, blade, assembly, env, tvec, surr; kwargs...) -> (aerostates, surr_history, mesh)

Allocating wrapper for [`run_sim_surrogate!`](@ref). Equivalent to
calling [`initialize_sim_surrogate`](@ref) and then `run_sim_surrogate!`.

For repeated solves prefer the in-place pair to reuse buffers.
"""
function run_sim_surrogate(rotor::Rotor, blade::Blade, assembly::GXBeam.Assembly,
                           env::Environment, tvec, surr::AbstractStructuralSurrogate; kwargs...)
    aerostates, surr_history, mesh = initialize_sim_surrogate(blade, assembly, tvec)
    run_sim_surrogate!(rotor, blade, mesh, env, tvec, aerostates, surr_history, surr; kwargs...)
    return aerostates, surr_history, mesh
end
