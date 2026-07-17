


"""
    Rotor{TI, TF, TB, T1, T2, T3, T4}

Rotor-level properties. Modeled after `CCBlade.Rotor`.

**Fields**
- `B::Int`: number of blades.
- `hubht::Float`: hub height in the global reference frame (m).
- `tilt::Float`: tilt angle (rad). Positive about the negative X axis (tilts the turbine up).
- `yaw::Float`: yaw angle (rad). Positive about the negative Y axis (turns the turbine clockwise).
- `turbine::Bool`: `true` for wind-turbine sign conventions; `false` otherwise.
- `mach::CCBlade.MachCorrection`: Mach-number correction model, or `nothing`.
- `re::CCBlade.ReCorrection`: Reynolds-number correction model, or `nothing`.
- `rotation::CCBlade.RotationCorrection`: rotation correction model, or `nothing`.
- `tip::CCBlade.TipCorrection`: hub/tip-loss correction model, or `nothing`.
"""
struct Rotor{TI, TF, TB, T1 <: Union{Nothing, CCBlade.MachCorrection}, T2 <: Union{Nothing, CCBlade.ReCorrection}, T3 <: Union{Nothing, CCBlade.RotationCorrection}, T4 <: Union{Nothing, CCBlade.TipCorrection}}
    B::TI
    hubht::TF
    tilt::TF
    yaw::TF
    turbine::TB
    mach::T1
    re::T2
    rotation::T3
    tip::T4
end

"""
    Rotor(B, hubht, turbine; tilt=0.0, yaw=0.0, mach=nothing, re=nothing, rotation=nothing, tip=nothing) -> Rotor

Convenience constructor for [`Rotor`](@ref) with all geometry/correction
options as keyword arguments and zero-default tilt/yaw.

**Arguments**
- `B::Int`: number of blades.
- `hubht::Float`: hub height (m).
- `turbine::Bool`: `true` for wind-turbine sign conventions.

**Keyword Arguments**
- `tilt::Float = 0.0`: tilt angle (rad).
- `yaw::Float = 0.0`: yaw angle (rad).
- `mach::CCBlade.MachCorrection`: Mach-number correction model, or `nothing`.
- `re::CCBlade.ReCorrection`: Reynolds-number correction model, or `nothing`.
- `rotation::CCBlade.RotationCorrection`: rotation correction model, or `nothing`.
- `tip::CCBlade.TipCorrection`: hub/tip-loss correction model, or `nothing`.
"""
function Rotor(B, hubht, turbine; tilt=0.0, yaw=0.0, mach=nothing, re=nothing, rotation=nothing, tip=nothing)
    return Rotor(B, hubht, tilt, yaw, turbine, mach, re, rotation, tip)
end


"""
    Blade{TF, TF2}

Per-blade geometry, twist, and airfoil definitions. All aerodynamic-node
quantities are aligned with the blade reference frame (precone, sweep,
and curve applied; tilt and yaw are *not*).

**Fields**
- `rhub::TF`: radial distance from the rotor center to the hub edge (m).
- `rtip::TF`: radial distance from the rotor center to the blade tip (m).
- `rx::AbstractVector{<:TF}`: x-position of the aero nodes (lead-lag / freestream direction, includes curve).
- `ry::AbstractVector{<:TF}`: y-position of the aero nodes (flapwise direction, includes sweep).
- `rz::AbstractVector{<:TF}`: z-position of the aero nodes (radial direction).
- `r::AbstractVector{<:TF}`: vector-norm radial position `‖(rx, ry, rz)‖` of each aero node.
- `rR::AbstractVector{<:TF}`: fractional blade-length location of each aero node (`r/rtip`).
- `c::AbstractVector{<:TF2}`: chord length at each aero node.
- `thetax::AbstractVector{<:TF}`: sweep angle of each local element about the X axis (rad).
- `thetay::AbstractVector{<:TF}`: curve angle of each local element about the Y axis (rad).
- `twist::AbstractVector{<:TF2}`: twist distribution (rad).
- `precone::TF`: blade precone angle (rad).
- `xcp::AbstractVector{<:TF}`: per-section center-of-pressure location.
- `airfoils::AbstractVector{<:DS.Airfoil}`: airfoil instance at each aero node.
"""
struct Blade{TF, TF2}
    rhub::TF
    rtip::TF
    rx::AbstractVector{<:TF} #Lead-lag direction (freestream) curve value
    ry::AbstractVector{<:TF} #Flapwise direction (Sweep value)
    rz::AbstractVector{<:TF} #Radial direction
    r::AbstractVector{<:TF} #Todo: Decide if I want this in here, or if I'll just use the norm of the vectorized version. 
    rR::AbstractVector{<:TF}
    c::AbstractVector{<:TF2} #Chord length
    thetax::AbstractVector{<:TF} #Sweep angle
    thetay::AbstractVector{<:TF} #Curve angle
    twist::AbstractVector{<:TF2}
    precone::TF
    xcp::AbstractVector{<:TF} 
    airfoils::AbstractVector{<:DS.Airfoil}
end

"""
    Blade(span, chord, twist, xcp, airfoils; rhub=span[1], rtip=span[end], precone=0.0, sweep=0.0, curve=0.0, rx=zero(span), ry=zero(span)) -> Blade

Convenience constructor for [`Blade`](@ref). Issues warnings when input
angles look like they were passed in degrees instead of radians, and
broadcasts scalar `sweep` / `curve` to per-node vectors.

**Arguments**
- `span::AbstractVector`: radial-direction distances of the aerodynamic nodes (m).
- `chord::AbstractVector`: chord length at each node.
- `twist::AbstractVector`: twist at each node (rad).
- `xcp::AbstractVector`: center-of-pressure location at each node.
- `airfoils::AbstractVector{<:DS.Airfoil}`: one airfoil instance per node.

**Keyword Arguments**
- `rhub`: radial distance to the hub. Defaults to `span[1]`.
- `rtip`: blade tip radius. Defaults to `span[end]`.
- `precone`: blade precone angle (rad).
- `sweep`: per-node sweep angle (rad), or a scalar broadcast to all nodes.
- `curve`: per-node curve angle (rad), or a scalar broadcast to all nodes.
- `rx`: x-direction distances (curve offsets) of the aerodynamic nodes. Defaults to zeros.
- `ry`: y-direction distances (sweep offsets) of the aerodynamic nodes. Defaults to zeros.
"""
function Blade(span, chord, twist, xcp, airfoils::AbstractVector{<:DS.Airfoil}; rhub=span[1], rtip=span[end], precone=0.0, sweep=0.0, curve=0.0, rx=zero(span), ry=zero(span))
    n = length(airfoils)

    if length(span)!=length(twist)!=length(chord)!=length(xcp)!=n
        error("Blade(): The number of airfoils and nodes (span) but be the same.")
    end

    if length(sweep)==1
        sweep = sweep.*ones(n)
    end

    if length(curve)==1
        curve = curve.*ones(n)
    end

    if length(sweep)!=length(curve)!=n
        error("Blade(): The length of the sweep and curve vectors must be as long as the radial node vector.")
    end

    pi2 = pi/2

    if precone>pi2
        @warn("Blade(): The precone angle you provided appears to be degrees, not radians.")
    elseif any(i->i>pi2, sweep)
        @warn("Blade(): The sweep angle(s) you provided appears to be degrees, not radians.")
    elseif any(i->i>pi2, curve)
        @warn("Blade(): The curve angle(s) you provided appears to be degrees, not radians.")
    end

    if any(i->i>pi2, twist)
        @warn("Blade(): The twist angle appears to be in degrees, not radians.")
    end

    rvec = @. sqrt(rx^2 + ry^2 + span^2)
    rR = span./rtip

    return Blade(rhub, rtip, rx, ry, span, rvec, rR, chord, sweep, curve, twist, precone, xcp, airfoils)
end



"""
    AeroStates{TF}

Aerodynamic state history of a transient simulation. All fields are mutated in
place by the time-stepping loops; the struct itself is immutable because no
field reference is reassigned.

**Fields**
- `azimuth::Vector{TF}`: Azimuthal position at each of the `nt` time steps.
- `phi::Matrix{TF}`: Inflow angle `nt × na` (rad).
- `alpha::Matrix{TF}`: Angle of attack `nt × na` (rad).
- `W::Matrix{TF}`: Inflow velocity magnitude `nt × na` (m/s).
- `Cx::Matrix{TF}`, `Cy::Matrix{TF}`, `Cm::Matrix{TF}`: Force/moment coefficients `nt × na`.
- `Fx::Matrix{TF}`, `Fy::Matrix{TF}`, `Mx::Matrix{TF}`: Dimensional sectional loads `nt × na`.
- `xds::Matrix{TF}`: Dynamic stall state history `nt × ns`.
"""
struct AeroStates{TF}
    azimuth::Vector{TF}
    phi::Matrix{TF}
    alpha::Matrix{TF}
    W::Matrix{TF}
    Cx::Matrix{TF}
    Cy::Matrix{TF}
    Cm::Matrix{TF}
    Fx::Matrix{TF}
    Fy::Matrix{TF}
    Mx::Matrix{TF}
    xds::Matrix{TF}
end

Base.eltype(::AeroStates{TF}) where {TF} = TF
Base.eltype(::Type{AeroStates{TF}}) where {TF} = TF

"""
    AeroStates{TF}(undef, nt, na, ns) -> aerostates

Allocate an uninitialized `AeroStates{TF}` with the standard shapes for a
transient simulation: `nt` time steps, `na` aerodynamic nodes, `ns` total DS
states across the blade.
"""
function AeroStates{TF}(::UndefInitializer, nt::Integer, na::Integer, ns::Integer) where {TF}
    AeroStates{TF}(
        Vector{TF}(undef, nt),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, na),
        Matrix{TF}(undef, nt, ns),
    )
end

"""
    StaticAeroStates{TF}

Aerodynamic state at a single steady-state instant — no time dimension, no DS
states. Used by the static fixed-point solver.

**Fields**
- `phi::Vector{TF}`, `alpha::Vector{TF}`, `W::Vector{TF}`: BEM solution per node.
- `Cx::Vector{TF}`, `Cy::Vector{TF}`, `Cm::Vector{TF}`: Force/moment coefficients per node.
- `Fx::Vector{TF}`, `Fy::Vector{TF}`, `Mx::Vector{TF}`: Dimensional sectional loads per node.
"""
struct StaticAeroStates{TF}
    phi::Vector{TF}
    alpha::Vector{TF}
    W::Vector{TF}
    Cx::Vector{TF}
    Cy::Vector{TF}
    Cm::Vector{TF}
    Fx::Vector{TF}
    Fy::Vector{TF}
    Mx::Vector{TF}
end

Base.eltype(::StaticAeroStates{TF}) where {TF} = TF
Base.eltype(::Type{StaticAeroStates{TF}}) where {TF} = TF

"""
    StaticAeroStates{TF}(undef, na) -> aerostates

Allocate an uninitialized `StaticAeroStates{TF}` for `na` aerodynamic nodes.
"""
function StaticAeroStates{TF}(::UndefInitializer, na::Integer) where {TF}
    StaticAeroStates{TF}(
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
        Vector{TF}(undef, na),
    )
end


"""
    SurrogatePointState{TF}

Per-structural-point state predicted by a structural surrogate.

Duck-compatible with `GXBeam.PointState` on the three fields the aero
coupling reads (`u`, `theta`, `V`); also carries `Omega`, `F`, `M` for
downstream analyses (fatigue, failure checks).

**Fields**
- `u::SVector{3,TF}`     — linear displacement (aero coupling)
- `theta::SVector{3,TF}` — Wiener–Milenkovic rotation (aero coupling)
- `V::SVector{3,TF}`     — linear velocity (aero coupling)
- `Omega::SVector{3,TF}` — angular velocity (post-processing)
- `F::SVector{3,TF}`     — internal force (post-processing)
- `M::SVector{3,TF}`     — internal moment (post-processing)
"""
struct SurrogatePointState{TF}
    u::SVector{3,TF}
    theta::SVector{3,TF}
    V::SVector{3,TF}
    Omega::SVector{3,TF}
    F::SVector{3,TF}
    M::SVector{3,TF}
end

"""
    SurrogateAssemblyState{TF}

Per-step assembly state produced by a structural surrogate's `decode`.
Holds one [`SurrogatePointState`](@ref) per structural point.

**Fields**
- `points::Vector{SurrogatePointState{TF}}`
"""
struct SurrogateAssemblyState{TF}
    points::Vector{SurrogatePointState{TF}}
end

Base.eltype(::SurrogateAssemblyState{TF}) where {TF} = TF
Base.eltype(::Type{SurrogateAssemblyState{TF}}) where {TF} = TF

"""
    zero_surrogate_state(TF, np) -> SurrogateAssemblyState{TF}

Allocate a `SurrogateAssemblyState{TF}` with `np` points, all fields zero.
"""
function zero_surrogate_state(::Type{TF}, np::Integer) where {TF}
    z = SVector{3,TF}(zero(TF), zero(TF), zero(TF))
    pts = [SurrogatePointState{TF}(z, z, z, z, z, z) for _ in 1:np]
    return SurrogateAssemblyState{TF}(pts)
end

"""
    AbstractStructuralSurrogate

Supertype for user-supplied structural surrogates plugged into
[`run_sim_surrogate!`](@ref). Concrete subtypes must implement the three
pure methods `encode_initial`, `step`, and `decode`. See the
`run_sim_surrogate!` docstring for the interface contract.
"""
abstract type AbstractStructuralSurrogate end

"""
    encode_initial(surr::AbstractStructuralSurrogate, u0_struct::SurrogateAssemblyState)

Return the initial latent state `z0` from a physical-state IC on the
structural points. Pure — must not mutate `surr`.
"""
function encode_initial end

"""
    step_latent(surr::AbstractStructuralSurrogate, z, f_per_element::AbstractMatrix)

Advance the latent state one step under the per-element load matrix
`f_per_element` (size `nelem × 6`, columns are `Fx, Fy, Fz, Mx, My, Mz`).
Pure — returns the new latent state, must not mutate `surr` or `z`.

Named `step_latent` (not `step`) to avoid shadowing `Base.step`.
"""
function step_latent end

"""
    decode(surr::AbstractStructuralSurrogate, z) -> SurrogateAssemblyState

Decode the latent state `z` into a fresh [`SurrogateAssemblyState`](@ref)
on the structural-point grid. Pure.
"""
function decode end

"""
    decode!(surr, z, u, theta, V, Omega, F, M)

Device / batched decode used by [`run_sim_surrogate_gpu!`](@ref). Writes the
decoded structural point states in place into six preallocated `(3, np, ns)`
buffers (`u`, `theta`, `V`, `Omega`, `F`, `M`) from the batched latent state
`z` `(nlatent, ns)`. Companion to the scalar/host [`decode`](@ref); kept
separate so the GPU path can stay allocation-free per step and never leaves the
device. Concrete GPU surrogates implement this method.
"""
function decode! end
