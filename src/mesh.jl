#=
functions to transfer loads to the structural mesh and interpolate structural displacements/velocities to the aerodynamic mesh. Also includes the SimMesh structs that hold the coupling state and scratch buffers for the different simulation modes (aero-only transient, coupled transient, steady aero-structural).
=#

"""
    find_point_indices(rvec, r)

Find the indices of the two nearest neighboring nodes in a 1D mesh of a given position. Note that this function assumes that the mesh is monatonically increasing.

**Arguments**
- `rvec::Vector{Number}`: The 1D mesh of interest.
- r - the position of the point of interest

**Returns**
- pair - a tuple containing the indices of the point before and after the aerodynamic point (before, after)
"""
function find_point_indices(rvec, r)

    n1 = length(rvec)

    if r<rvec[1]
        @warn("find_point_indices(): The point was not found in the 1D mesh.")
        return 1, 2
    elseif r>rvec[end]
        @warn("find_point_indices(): The point was not found in the 1D mesh.")
        return (n1-1, n1)
    end

    ### Iterate through the points 
    for i = 1:n1-1
        if rvec[i]>r
            return (i-1, i)
        end
    end

    return (n1-1, n1)
end


"""
    find_interpolation_percent(rvec, pair, r) -> Number

Fractional position of `r` between the two `rvec` nodes identified by `pair`.

**Arguments**
- `rvec`: 1D mesh containing the neighboring nodes.
- `pair::Tuple{Int, Int}`: indices of the bracketing neighbors.
- `r`: position of interest on the same axis as `rvec`.

**Returns**
- Fractional location in `[0, 1]` (0 at `rvec[pair[1]]`, 1 at `rvec[pair[2]]`).
"""
function find_interpolation_percent(rvec, pair, r)

    a = r - rvec[pair[1]]
    L = rvec[pair[2]]- rvec[pair[1]]
    return a/L
end

"""
    InterpolationPoint{TF}(pair, percent)

A precomputed 1D linear-interpolation handle: the two source-mesh indices
that bracket a target point, plus the fractional position between them.
Cached so the structural→aero interpolation doesn't repeat the bracket
search every time step.

**Fields**
- `pair::Tuple{Int64, Int64}`: indices of the bracketing source-mesh nodes.
- `percent::TF`: fractional position in `[0, 1]` between `pair[1]` and `pair[2]`.
"""
struct InterpolationPoint{TF}
    pair::Tuple{Int64, Int64}
    percent::TF
end

#=
Note: Maybe down the line I can switch InterpolationPoint to a type, and make different structs of linear, quadratic, cubic, or maybe even more general of polynomial, Lagrange, Legendre and interpolate that way.

=#


"""
    AbstractSimMesh

Supertype for the simulation mesh structs that hold the per-step scratch
buffers and coupling state between aero, structures, and dynamic stall.
Concrete subtypes: [`SimMesh`](@ref) (full coupled transient), [`AeroMesh`](@ref)
(aero-only transient), [`StaticMesh`](@ref) (steady-state fixed point).
"""
abstract type AbstractSimMesh end

"""
    AeroMesh{TF, TXDS, TYDS, TPDS}

Scratch state for an aero-only transient simulation. No structural coupling
fields — `delta`/`def_theta`/`aerov` exist solely so `take_aero_step!` can
reuse the same code path as the coupled solver; they stay zero.

**Fields**
- `delta::Vector{SVector{3,TF}}`: Always zero in aero-only mode.
- `def_theta::Vector{SVector{3,TF}}`: Always zero in aero-only mode.
- `aerov::Vector{SVector{3,TF}}`: Always zero in aero-only mode.
- `xcc::Vector{TF}`: CCBlade scratch vector.
- `xds_idxs::TXDS`: DS state-vector slice indices, per airfoil.
- `y_ds::TYDS`: DS environment inputs (W, Udot, alpha, alphadot per node).
- `p_ds::TPDS`: DS parameters (chord, x_cp per node).
"""
struct AeroMesh{TF, TXDS, TYDS, TPDS} <: AbstractSimMesh
    delta::Vector{SVector{3, TF}}
    def_theta::Vector{SVector{3, TF}}
    aerov::Vector{SVector{3, TF}}
    xcc::Vector{TF}
    xds_idxs::TXDS
    y_ds::TYDS
    p_ds::TPDS
end

"""
    SimMesh{TF, TIP, TXDS, TYDS, TPDS, TA, TS, TPC, TDL, TPM, TXP, TPF}

Scratch state for the full coupled aero-structural transient simulation.
Carries the aero-side fields of [`AeroMesh`](@ref) plus the GXBeam coupling
state.

**Fields**
- `interpolationpoints`, `delta`, `def_theta`, `aerov`, `xcc`, `xds_idxs`,
  `y_ds`, `p_ds`: As in [`AeroMesh`](@ref).
- `assembly::TA`: GXBeam assembly.
- `system::TS`: GXBeam `DynamicSystem` reused across time steps.
- `prescribed_conditions::TPC`: GXBeam boundary conditions.
- `distributed_loads::TDL`: Updated each step from BEM/DS loads.
- `point_masses::TPM`: GXBeam point masses (empty in standard use).
- `linear_velocity`, `angular_velocity::SVector{3,Float64}`: Base body motion.
- `xpfunc::TXP`, `pfunc::TPF`: GXBeam parameter callbacks.
- `two_dimensional::Bool`, `structural_damping::Bool`, `linear::Bool`: GXBeam flags.
"""
struct SimMesh{TF, TIP, TXDS, TYDS, TPDS, TA, TS, TPC, TDL, TPM, TXP, TPF} <: AbstractSimMesh
    interpolationpoints::TIP
    delta::Vector{SVector{3, TF}}
    def_theta::Vector{SVector{3, TF}}
    aerov::Vector{SVector{3, TF}}
    xcc::Vector{TF}
    xds_idxs::TXDS
    y_ds::TYDS
    p_ds::TPDS
    assembly::TA
    system::TS
    prescribed_conditions::TPC
    distributed_loads::TDL
    point_masses::TPM
    linear_velocity::SVector{3, Float64}
    angular_velocity::SVector{3, Float64}
    xpfunc::TXP
    pfunc::TPF
    two_dimensional::Bool
    structural_damping::Bool
    linear::Bool
end

"""
    StaticMesh{TF, TIP, TA, TS, TPC, TDL, TPM, TXP, TPF}

Scratch state for the steady-state fixed-point aero-structural solver. Like
[`SimMesh`](@ref) but with no dynamic-stall fields (the static solver has no
DS history).
"""
struct StaticMesh{TF, TIP, TA, TS, TPC, TDL, TPM, TXP, TPF} <: AbstractSimMesh
    interpolationpoints::TIP
    delta::Vector{SVector{3, TF}}
    def_theta::Vector{SVector{3, TF}}
    aerov::Vector{SVector{3, TF}}
    xcc::Vector{TF}
    assembly::TA
    system::TS
    prescribed_conditions::TPC
    distributed_loads::TDL
    point_masses::TPM
    linear_velocity::SVector{3, Float64}
    angular_velocity::SVector{3, Float64}
    xpfunc::TXP
    pfunc::TPF
    two_dimensional::Bool
    structural_damping::Bool
    linear::Bool
end

"""
    create_interpolationpoints(assembly, blade) -> Vector{InterpolationPoint}

Build a cache of [`InterpolationPoint`](@ref)s that maps each aerodynamic
station in `blade.r` onto the structural arc-length of `assembly`. Used
once during simulation setup so per-step calls to
[`interpolate_position`](@ref) / [`interpolate_deflection`](@ref) etc. can
skip the bracket search.

Assumes both the structural assembly points and `blade.r` are ordered
root→tip.

**Arguments**
- `assembly::GXBeam.Assembly`: structural assembly.
- `blade::Blade`: blade providing the aerodynamic station radii.

**Returns**
- `Vector{InterpolationPoint}` of length `length(blade.r)`.
"""
function create_interpolationpoints(assembly::GXBeam.Assembly, blade::Blade)

    ### Allocate memory
    na = length(blade.r)
    points = Vector{InterpolationPoint}(undef, na)

    ### Find the beam length 
    rgx = get_bladelength_vector(assembly)

    for i = 1:na
        pair = find_point_indices(rgx, blade.r[i])
        percent = find_interpolation_percent(rgx, pair, blade.r[i])
        points[i] = InterpolationPoint(pair, percent)
    end

    return points
end


"""
    interpolate_position(ip, assembly, state) -> SVector{3}

Position of the aerodynamic node represented by `ip` in the deformed
configuration: the linear blend of the two bracketing structural points
shifted by their current displacements.

**Arguments**
- `ip::InterpolationPoint`: cached interpolation handle for the aero node.
- `assembly::GXBeam.Assembly`: reference (undeformed) point positions.
- `state`: GXBeam state holding the per-point displacement `u`.

**Returns**
- Deformed position in the same frame as `assembly.points`.
"""
function interpolate_position(ip, assembly, state)  

    p1 = assembly.points[ip.pair[1]] + state.points[ip.pair[1]].u
    p2 = assembly.points[ip.pair[2]] + state.points[ip.pair[2]].u

    return (1-ip.percent)*p1 + ip.percent*p2
end

"""
    interpolate_deflection(ip, assembly, state) -> SVector{3}

Linearly blend the structural displacement `u` at the two bracketing
points to produce the displacement at the aerodynamic node. `assembly` is
accepted for signature symmetry with [`interpolate_position`](@ref) but
isn't used directly here.

**Arguments**
- `ip::InterpolationPoint`: cached interpolation handle.
- `assembly::GXBeam.Assembly`: structural assembly (unused).
- `state`: GXBeam state holding the per-point displacement `u`.

**Returns**
- Interpolated 3-component displacement.
"""
function interpolate_deflection(ip, assembly, state)
    p1 = state.points[ip.pair[1]].u
    p2 = state.points[ip.pair[2]].u

    return (1-ip.percent)*p1 + ip.percent*p2
end

"""
    interpolate_velocity(ip, assembly, state) -> SVector{3}

Linearly blend the structural velocity `V` at the two bracketing points to
produce the velocity at the aerodynamic node. Used inside
[`convert_velocities`](@ref) to subtract structural motion from the
inflow.

**Arguments**
- `ip::InterpolationPoint`: cached interpolation handle.
- `assembly::GXBeam.Assembly`: structural assembly (unused).
- `state`: GXBeam state holding the per-point velocity `V`.

**Returns**
- Interpolated 3-component velocity in the structural frame.
"""
function interpolate_velocity(ip, assembly, state)

    v1 = state.points[ip.pair[1]].V
    v2 = state.points[ip.pair[2]].V

    return (1-ip.percent)*v1 + ip.percent*v2
end

"""
    interpolate_angle(ip, assembly, state) -> SVector{3}

Linearly blend the Wiener–Milenkovic rotation parameters at the two
bracketing points to produce a 3-2-1 Euler-angle triple for the
aerodynamic node. The per-point rotations are first converted to Euler
angles via [`WMPtoangle`](@ref) before averaging.

**Arguments**
- `ip::InterpolationPoint`: cached interpolation handle.
- `assembly::GXBeam.Assembly`: structural assembly (unused).
- `state`: GXBeam state holding per-point WMP rotations `theta`.

**Returns**
- Interpolated 3-component Euler-angle vector (roll, pitch, yaw) in radians.
"""
function interpolate_angle(ip, assembly, state)
    
    theta1 = WMPtoangle(state.points[ip.pair[1]].theta)
    theta2 = WMPtoangle(state.points[ip.pair[2]].theta) 

    return (1-ip.percent)*theta1 + ip.percent*theta2
end

"""
    convert_velocities(blade, env, assembly, state, interpolationpoints, t, idx) -> NTuple{3}

Interpolate the structural displacement and velocity at aero station `idx`,
rotate them into the aerodynamic reference frame, and subtract the rigid
rotational contribution so what's returned is the *relative* velocity the
airfoil sees (i.e. structural-motion-induced inflow only — the freestream
is added downstream).

The axis flips (`dy = +delta[2]`, `dz = +delta[1]`, `usy = -Vs[2]`,
`usz = -Vs[1]`) account for the structural-frame → aero-frame rotation and
the sign convention that "structure moving in `+x`" is "wind moving in
`-x`" from the airfoil's perspective.

**Arguments**
- `blade::Blade`: blade providing the undeformed station coordinates `ry`, `rz`.
- `env::Environment`: environment supplying the rotor speed `env.RS(t)`.
- `assembly::GXBeam.Assembly`: structural assembly (forwarded to interpolators).
- `state`: current GXBeam state.
- `interpolationpoints::Vector{InterpolationPoint}`: precomputed handles for every aero station.
- `t`: current time (rotor speed is evaluated here).
- `idx::Int`: aero-station index.

**Returns**
- `(ux, uy, uz)`: relative velocity at the aerodynamic node in the aero frame.
"""
function convert_velocities(blade::Blade, env::Environment, assembly, state, interpolationpoints, t, idx)

    ### Interpolate the structural quantities to the aerodynamic mesh. 
    delta = interpolate_deflection(interpolationpoints[idx], assembly, state)
    Vs = interpolate_velocity(interpolationpoints[idx], assembly, state)
    

    ### Transform the quantities into the aerodynamic frame
    # dx = -delta[3]
    dy = delta[2]
    dz = delta[1]
    usx = Vs[3]
    usy = -Vs[2]
    usz = -Vs[1]
    #Note: The  velocities are negative of the transformation because as the structure moves in those directions, the wind moves in the opposite direction. 

    ### Remove the rotational velocities
    # rax = blade.rx[idx] 
    ray = blade.ry[idx] + dy
    raz = blade.rz[idx] + dz

    omega = env.RS(t)

    #Rotational velocities
    # urx = 0.0
    ury = omega*raz 
    urz = -omega*ray

    return (usx, usy-ury, usz-urz)
end

"""
    update_mesh!(blade, mesh, assembly, gxstate, env, t, na)

Refresh the per-step coupling buffers stored in `mesh` (`delta`,
`def_theta`, `aerov`) by interpolating the latest GXBeam state onto the
aerodynamic stations. Called once per time step ahead of the BEM/DS load
evaluation.

Mutates `mesh.delta`, `mesh.def_theta`, and `mesh.aerov` in place.

**Arguments**
- `blade::Blade`: blade with aerodynamic stations.
- `mesh`: simulation mesh carrying `interpolationpoints` and the buffers to be filled.
- `assembly::GXBeam.Assembly`: structural assembly (undeformed reference).
- `gxstate`: current GXBeam state.
- `env::Environment`: environment for the rotational-velocity subtraction.
- `t`: current time.
- `na::Int`: number of aerodynamic stations to update.
"""
function update_mesh!(blade::Blade, mesh, assembly::GXBeam.Assembly, gxstate, env::Environment, t, na)

    for j = 1:na
        mesh.delta[j] = interpolate_deflection(mesh.interpolationpoints[j], assembly, gxstate)

        mesh.def_theta[j] = interpolate_angle(mesh.interpolationpoints[j], assembly, gxstate)

        mesh.aerov[j] = convert_velocities(blade, env, assembly, gxstate, mesh.interpolationpoints, t, j)
    end
end



"""
    transform_BC_G(rhr_x, rhr_y, rhr_z, azimuth, precone, tilt, yaw) -> NTuple{3}

Transform a 3-vector from the **blade-coned** frame (after precone, before
azimuth/tilt/yaw) into the **global** ground frame. The rotation order is
precone → azimuth → tilt → yaw.

**Arguments**
- `rhr_x`, `rhr_y`, `rhr_z`: components in the blade-coned frame.
- `azimuth`, `precone`, `tilt`, `yaw`: rotation angles in radians.

**Returns**
- `(rg_x, rg_y, rg_z)`: components in the global frame.
"""
function transform_BC_G(rhr_x, rhr_y, rhr_z, azimuth, precone, tilt, yaw)

    rtx, rty, rtz = rotate_y(rhr_x, rhr_y, rhr_z, precone; T=false)
    rtx, rty, rtz = rotate_x(rtx, rty, rtz, azimuth; T=true)
    rtx, rty, rtz = rotate_y(rtx, rty, rtz, tilt; T=true)
    rg_x, rg_y, rg_z = rotate_z(rtx, rty, rtz, yaw; T=true)

    return rg_x, rg_y, rg_z
end

"""
    transform_BC_HR(rbc_x, rbc_y, rbc_z, precone) -> NTuple{3}

Transform a 3-vector from the **blade-coned (BC)** frame into the
**hub-rotating (HR)** frame by applying the precone rotation about the
y-axis.

**Arguments**
- `rbc_x`, `rbc_y`, `rbc_z`: components in the blade-coned frame.
- `precone`: precone angle in radians.

**Returns**
- `(rhr_x, rhr_y, rhr_z)`: components in the hub-rotating frame.
"""
function transform_BC_HR(rbc_x, rbc_y, rbc_z, precone)

    rhr_x, rhr_y, rhr_z = rotate_y(rbc_x, rbc_y, rbc_z, precone; T=false)

    return rhr_x, rhr_y, rhr_z
end

"""
    transform_HR_L(rhr_x, rhr_y, rhr_z, curve, sweep, precone) -> NTuple{3}

Transform a 3-vector from the **hub-rotating (HR)** frame into the
**local-element (L)** frame by undoing sweep, curve, and precone in turn.

**Arguments**
- `rhr_x`, `rhr_y`, `rhr_z`: components in the hub-rotating frame.
- `curve`, `sweep`, `precone`: blade-local geometry angles in radians.

**Returns**
- `(rl_x, rl_y, rl_z)`: components in the local-element frame.
"""
function transform_HR_L(rhr_x, rhr_y, rhr_z, curve, sweep, precone)

    rt_x, rt_y, rt_z = rotate_x(rhr_x, rhr_y, rhr_z, sweep; T=true)
    rt_x, rt_y, rt_z = rotate_y(rt_x, rt_y, rt_z, curve; T=true)
    rl_x, rl_y, rl_z = rotate_y(rt_x, rt_y, rt_z, precone; T=false)

    return rl_x, rl_y, rl_z
end

"""
    transform_G_L(rg_x, rg_y, rg_z, azimuth, curve, precone, sweep, tilt, yaw) -> NTuple{3}

Composite transform from the **global** ground frame all the way down to
the **local-element (L)** frame. Applies yaw → tilt → azimuth → sweep →
curve → precone in order, undoing each step against the incoming
rotation.

**Arguments**
- `rg_x`, `rg_y`, `rg_z`: components in the global frame.
- `azimuth`, `curve`, `precone`, `sweep`, `tilt`, `yaw`: rotation angles in radians.

**Returns**
- `(rl_x, rl_y, rl_z)`: components in the local-element frame.
"""
function transform_G_L(rg_x, rg_y, rg_z, azimuth, curve, precone, sweep, tilt, yaw)

    rt_x, rt_y, rt_z = rotate_z(rg_x, rg_y, rg_z, yaw)
    rt_x, rt_y, rt_z = rotate_y(rt_x, rt_y, rt_z, tilt)
    rt_x, rt_y, rt_z = rotate_x(rt_x, rt_y, rt_z, azimuth)
    rt_x, rt_y, rt_z = rotate_x(rt_x, rt_y, rt_z, sweep; T=true)
    rt_x, rt_y, rt_z = rotate_y(rt_x, rt_y, rt_z, curve; T=true)
    rl_x, rl_y, rl_z = rotate_y(rt_x, rt_y, rt_z, precone; T=true)

    return rl_x, rl_y, rl_z
end