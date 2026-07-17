#=
GPU-resident coupled aeroelastic transient solver (surrogate structural model).

Fuses the GPU-batched BEMT (bemt_gpu.jl) and dynamic-stall (dsmodel_gpu.jl)
kernels on the aero side with a user-supplied structural surrogate on the
structural side, keeping every per-step state array resident on the device
across the whole time march. The `(section, sim)` grid batches a population of
simulations (each sim may carry its own surrogate conditioning / inflow).

This file provides:
  * device-safe ports of the aero↔structure coupling math
    (`aero_velocities`, `update_mesh!`, `dimensionalize!`, `extract_ds_loads!`,
     `build_surrogate_loads`), each as a KernelAbstractions kernel over
     `(section, sim)` or `(element, sim)`;
  * an abstract batched-surrogate interface;
  * the `initialize_sim_surrogate_gpu` / `run_sim_surrogate_gpu!` orchestration.

Forward pass only (matches the aero kernels). All kernels call the same pure
scalar helpers as the CPU path, so a CPU-backend run reproduces the CPU
functions to machine precision (validated piecewise in Phase 2/3).

Adam Cardoza
=#

using KernelAbstractions
using Adapt

export AeroGeometryGPU, InterpGPU, GPUPointStates
export aero_velocities_gpu!, ds_inputs_gpu!, ds_loads_dimensionalize_gpu!,
       build_surrogate_loads_gpu!, update_feedback_gpu!

# ---------------------------------------------------------------------------
# GPUPointStates — batched decoded structural point states, all device-resident.
# Each field is (3, np, ns). This is both the input to a GPU surrogate's
# `encode_initial` (the IC seed) and the in-place output target of `decode!`.
# ---------------------------------------------------------------------------

"""
    GPUPointStates{TA<:AbstractArray}

Batched structural point states on the device, `(3, np, ns)` per field. Used as
the IC seed passed to a GPU surrogate's [`encode_initial`](@ref) and as the
in-place target of [`decode!`](@ref). Fields mirror [`SurrogatePointState`](@ref):
`u` (displacement), `theta` (WM rotation params), `V`, `Omega`, `F`, `M`.
"""
struct GPUPointStates{TA<:AbstractArray}
    u::TA
    theta::TA
    V::TA
    Omega::TA
    F::TA
    M::TA
end

"""
    GPUPointStates(np, ns; ArrayType=Array{Float64})

Allocate zeroed `(3, np, ns)` device buffers for every field.
"""
function GPUPointStates(np::Integer, ns::Integer; ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    mk() = to_backend_array(ArrayType, zeros(TF, 3, np, ns))
    return GPUPointStates(mk(), mk(), mk(), mk(), mk(), mk())
end

Adapt.adapt_structure(to, p::GPUPointStates) = GPUPointStates(
    adapt(to, p.u), adapt(to, p.theta), adapt(to, p.V),
    adapt(to, p.Omega), adapt(to, p.F), adapt(to, p.M),
)

# ---------------------------------------------------------------------------
# Device-safe WMP → Euler-angle conversion.
#
# Port of `WMPtoangle` (gxbeam.jl) + `GXBeam.wiener_milenkovic` /
# `rotation_parameter_scaling` (GXBeam/src/math.jl) as pure scalar arithmetic.
# `WMPtoangle` returns the Euler angles of R = C', where C is the WM transform
# matrix; we build the needed entries of C directly and read the 3-2-1 angles:
#   roll  = atan(C[2,3], C[3,3])
#   pitch = asin(-C[1,3])
#   yaw   = atan(C[1,2], C[1,1])
# ---------------------------------------------------------------------------
@inline function wmp_to_angle_dev(c1, c2, c3)
    TF = typeof(c1)
    # rotation_parameter_scaling: extend range past 360°.
    cnorm = sqrt(c1 * c1 + c2 * c2 + c3 * c3)
    # m = div(div(cnorm,4)+1, 2) — integer floor divisions.
    m = floor((floor(cnorm / TF(4)) + one(TF)) / TF(2))
    scaling = cnorm > TF(4) ? one(TF) - TF(8) * m / cnorm : one(TF)

    s1 = scaling * c1
    s2 = scaling * c2
    s3 = scaling * c3

    c0 = TF(2) - (s1 * s1 + s2 * s2 + s3 * s3) / TF(8)
    den = one(TF) / (TF(4) - c0)^2

    C11 = den * (c0 * c0 + s1 * s1 - s2 * s2 - s3 * s3)
    C12 = den * TF(2) * (s1 * s2 + c0 * s3)
    C13 = den * TF(2) * (s1 * s3 - c0 * s2)
    C23 = den * TF(2) * (s2 * s3 + c0 * s1)
    C33 = den * (c0 * c0 - s1 * s1 - s2 * s2 + s3 * s3)

    roll  = atan(C23, C33)
    pitch = asin(clamp(-C13, -one(TF), one(TF)))
    yaw   = atan(C12, C11)
    return roll, pitch, yaw
end

# ---------------------------------------------------------------------------
# Kernel-safe geometry / interpolation constants (uploaded once).
# ---------------------------------------------------------------------------

"""
    AeroGeometryGPU{TF, TA<:AbstractVector}

Per-section blade geometry + rotor/environment scalars the coupling kernels
need, packed into device vectors. Built from a [`Blade`](@ref), [`Rotor`](@ref),
and `SimpleEnvironment` shear exponent. Move with `Adapt.adapt`.

**Fields**
- `rx, ry, rz::TA` — undeformed aero-node coordinates (lead-lag / flap / radial).
- `thetax, thetay::TA` — sweep / curve angles (rad).
- `chord::TA` — chord at each aero node.
- `twist::TA` — twist at each aero node (rad).
- `precone, hubht, tilt, yaw, shearexp::TF` — rotor/inflow scalars.
"""
struct AeroGeometryGPU{TF, TA<:AbstractVector}
    rx::TA
    ry::TA
    rz::TA
    thetax::TA
    thetay::TA
    chord::TA
    twist::TA
    precone::TF
    hubht::TF
    tilt::TF
    yaw::TF
    shearexp::TF
end

"""
    AeroGeometryGPU(blade, rotor, env; ArrayType=Array{Float64})

Build an [`AeroGeometryGPU`](@ref). Uses `rotor.yaw` (not negated) to match
`get_aerostructural_velocities`, the coupled-march velocity routine.
"""
function AeroGeometryGPU(blade::Blade, rotor::Rotor, env::Environment; ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    v(x) = to_backend_vector(ArrayType, TF.(x))
    return AeroGeometryGPU(
        v(blade.rx), v(blade.ry), v(blade.rz),
        v(blade.thetax), v(blade.thetay),
        v(blade.c), v(blade.twist),
        TF(blade.precone), TF(rotor.hubht), TF(rotor.tilt), TF(rotor.yaw), TF(env.shearexp),
    )
end

Adapt.adapt_structure(to, g::AeroGeometryGPU) = AeroGeometryGPU(
    adapt(to, g.rx), adapt(to, g.ry), adapt(to, g.rz),
    adapt(to, g.thetax), adapt(to, g.thetay),
    adapt(to, g.chord), adapt(to, g.twist),
    g.precone, g.hubht, g.tilt, g.yaw, g.shearexp,
)

AeroGeometryGPU(rx, ry, rz, tx, ty, c, tw, pc, hh, tl, yw, se) =
    AeroGeometryGPU{typeof(pc), typeof(rx)}(rx, ry, rz, tx, ty, c, tw, pc, hh, tl, yw, se)

"""
    InterpGPU{TIV<:AbstractVector, TV<:AbstractVector}

Structural→aero linear-interpolation handles for every aero node, ported from
the cached [`InterpolationPoint`](@ref) list into three device vectors.

**Fields**
- `pair1, pair2::TIV` — 1-based bracketing structural-point indices (Int32).
- `percent::TV` — fractional position in `[0,1]` between `pair1` and `pair2`.
"""
struct InterpGPU{TIV<:AbstractVector, TV<:AbstractVector}
    pair1::TIV
    pair2::TIV
    percent::TV
end

"""
    InterpGPU(interpolationpoints; ArrayType=Array{Float64})

Pack a `Vector{InterpolationPoint}` (from [`create_interpolationpoints`](@ref))
into device vectors.
"""
function InterpGPU(ips::AbstractVector; ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    na = length(ips)
    p1 = Int32[ips[j].pair[1] for j in 1:na]
    p2 = Int32[ips[j].pair[2] for j in 1:na]
    pc = TF[ips[j].percent for j in 1:na]
    return InterpGPU(
        to_backend_vector(similar_type(ArrayType, Int32), p1),
        to_backend_vector(similar_type(ArrayType, Int32), p2),
        to_backend_vector(ArrayType, pc),
    )
end

Adapt.adapt_structure(to, ip::InterpGPU) = InterpGPU(
    adapt(to, ip.pair1), adapt(to, ip.pair2), adapt(to, ip.percent),
)

InterpGPU(p1::AbstractVector, p2::AbstractVector, pc::AbstractVector) =
    InterpGPU{typeof(p1), typeof(pc)}(p1, p2, pc)

# ---------------------------------------------------------------------------
# Kernel 1 — aero velocities (port of get_aerostructural_velocities).
#
# Per (section j, sim s): reads structural feedback delta/def_theta/aerov,
# host-evaluated per-sim inflow (Ux,Uy,Uz) and rotor speed RS, azimuth, and the
# original pitch. Writes Vx, Vy (BEMT inputs) and pitch_eff = pitch - dθ1 (the
# twist-corrected pitch the BEMT solve uses).
# ---------------------------------------------------------------------------
@kernel function aero_velocities_kernel!(Vx, Vy, pitch_eff,
        delta, def_theta, aerov,
        rx, ry, rz, thetax, thetay, precone, hubht, tilt, yaw, shearexp,
        Ux, Uy, Uz, RS, azimuth, pitch)
    j, s = @index(Global, NTuple)
    @inbounds begin
        TF = eltype(Vx)

        d1 = delta[1, j, s]; d2 = delta[2, j, s]; d3 = delta[3, j, s]
        th1 = def_theta[1, j, s]; th2 = def_theta[2, j, s]; th3 = def_theta[3, j, s]
        vs1 = aerov[1, j, s]; vs2 = aerov[2, j, s]; vs3 = aerov[3, j, s]

        # structural deflection into aero frame
        dx = -d3; dy = d2; dz = d1
        rbc_x = rx[j] + dx
        rbc_y = ry[j] + dy
        rbc_z = rz[j] + dz

        sweep = -(thetax[j] - th3)
        curve = thetay[j] + th2
        pc = precone

        az = azimuth[s]; rs = RS[s]
        yw = yaw

        # --- freestream ---
        _, _, rgz = transform_BC_G(rbc_x, rbc_y, rbc_z, az, pc, tilt, yw)
        z = rgz + hubht
        factor = (z / hubht)^shearexp
        ugx = Ux[s] * factor; ugy = Uy[s] * factor; ugz = Uz[s] * factor
        ulx_wind, uly_wind, _ = transform_G_L(ugx, ugy, ugz, az, curve, pc - th2, sweep, tilt, yw)

        # --- rigid rotational velocity ---
        # cross((rs,0,0),(rhr_x,rhr_y,rhr_z)) = (0, -rs*rhr_z, rs*rhr_y)
        _, rhr_y, rhr_z = transform_BC_HR(rbc_x, rbc_y, rbc_z, pc)
        vxr = zero(TF)
        vyr = -rs * rhr_z
        vzr = rs * rhr_y
        ulx_rot, uly_rot, _ = transform_HR_L(vxr, vyr, vzr, curve, sweep, pc - th2)

        # --- structural velocity ---
        usx, usy, _ = transform_HR_L(vs1, vs2, vs3, curve, sweep, pc - th2)

        Vxa = ulx_wind - ulx_rot + usx
        Vya = uly_wind - uly_rot + usy

        st, ct = sincos(th1)
        Vx[j, s] = Vxa * ct - Vya * st
        Vy[j, s] = -Vxa * st + Vya * ct
        pitch_eff[j, s] = pitch - th1
    end
end

"""
    aero_velocities_gpu!(Vx, Vy, pitch_eff, geo, delta, def_theta, aerov,
                         Ux, Uy, Uz, RS, azimuth, pitch)

Compute the airfoil-frame relative velocities `Vx, Vy` and the twist-corrected
pitch `pitch_eff = pitch - def_theta_x` for every `(section, sim)`, given the
structural feedback (`delta/def_theta/aerov`, each `(3, na, ns)`) and the
host-evaluated per-sim inflow (`Ux/Uy/Uz`, `RS`, `azimuth` each length `ns`).
Device port of [`get_aerostructural_velocities`](@ref).
"""
function aero_velocities_gpu!(Vx, Vy, pitch_eff, geo::AeroGeometryGPU,
                              delta, def_theta, aerov,
                              Ux, Uy, Uz, RS, azimuth, pitch)
    na, ns = size(Vx)
    backend = KernelAbstractions.get_backend(Vx)
    k = aero_velocities_kernel!(backend)
    TF = eltype(Vx)
    k(Vx, Vy, pitch_eff, delta, def_theta, aerov,
      geo.rx, geo.ry, geo.rz, geo.thetax, geo.thetay,
      geo.precone, geo.hubht, geo.tilt, geo.yaw, geo.shearexp,
      Ux, Uy, Uz, RS, azimuth, TF(pitch); ndrange=(na, ns))
    KernelAbstractions.synchronize(backend)
    return Vx, Vy, pitch_eff
end

# ---------------------------------------------------------------------------
# Kernel 2 — DS inputs (aoa, U) from BEMT outputs (port of the aoa line in
# update_ds_inputs! / dsmodel_initial_condition!). Uses the ORIGINAL pitch
# (not pitch_eff), matching the CPU path.
# ---------------------------------------------------------------------------
@kernel function ds_inputs_kernel!(U_out, aoa_out, W, phi, twist, pitch, turbine)
    j, s = @index(Global, NTuple)
    @inbounds begin
        a = (twist[j] + pitch) - phi[j, s]
        aoa_out[j, s] = turbine ? -a : a
        U_out[j, s] = W[j, s]
    end
end

"""
    ds_inputs_gpu!(U_out, aoa_out, geo, W, phi, pitch, turbine)

Build the dynamic-stall inputs `(U = W, aoa = ±((twist+pitch) - phi))` for every
`(section, sim)` from the BEMT outputs. `aoa` is negated when `turbine` is true.
Port of the aoa construction in [`update_ds_inputs!`](@ref).
"""
function ds_inputs_gpu!(U_out, aoa_out, geo::AeroGeometryGPU, W, phi, pitch, turbine::Bool)
    na, ns = size(W)
    backend = KernelAbstractions.get_backend(W)
    k = ds_inputs_kernel!(backend)
    TF = eltype(W)
    k(U_out, aoa_out, W, phi, geo.twist, TF(pitch), turbine; ndrange=(na, ns))
    KernelAbstractions.synchronize(backend)
    return U_out, aoa_out
end

# ---------------------------------------------------------------------------
# Kernel 3 — DS loads → (Cx, Cy) rotation + dimensionalize.
# Port of extract_ds_loads! (Cl/Cd → Cx/Cy via phi) followed by dimensionalize!.
# ---------------------------------------------------------------------------
@kernel function ds_loads_dimensionalize_kernel!(Fx, Fy, Mx, Cl, Cd, Cm, phi, W, chord, rho)
    j, s = @index(Global, NTuple)
    @inbounds begin
        TF = eltype(Fx)
        sphi, cphi = sincos(phi[j, s])
        cl = Cl[j, s]; cd = Cd[j, s]; cm = Cm[j, s]
        Cx = cl * cphi + cd * sphi
        Cy = -(cl * sphi - cd * cphi)
        q = TF(0.5) * rho * W[j, s]^2
        c = chord[j]
        Fx[j, s] = Cx * q * c
        Fy[j, s] = Cy * q * c
        Mx[j, s] = cm * q * c * c
    end
end

"""
    ds_loads_dimensionalize_gpu!(Fx, Fy, Mx, geo, Cl, Cd, Cm, phi, W, rho)

Rotate the DS lift/drag coefficients into the rotor frame (`Cx/Cy`, port of
[`extract_ds_loads!`](@ref)) and dimensionalize into sectional loads
(`Fx/Fy/Mx`, port of [`dimensionalize!`](@ref)) for every `(section, sim)`.
"""
function ds_loads_dimensionalize_gpu!(Fx, Fy, Mx, geo::AeroGeometryGPU, Cl, Cd, Cm, phi, W, rho)
    na, ns = size(Fx)
    backend = KernelAbstractions.get_backend(Fx)
    k = ds_loads_dimensionalize_kernel!(backend)
    TF = eltype(Fx)
    k(Fx, Fy, Mx, Cl, Cd, Cm, phi, W, geo.chord, TF(rho); ndrange=(na, ns))
    KernelAbstractions.synchronize(backend)
    return Fx, Fy, Mx
end

# ---------------------------------------------------------------------------
# Kernel 4 — surrogate load matrix (port of build_surrogate_loads).
# f[j,1] = -g*m*cos(az);  f[j,2] = fy*L/2 - g*m*sin(az);  f[j,3] = -fx*L/2
# columns 4:6 zero. Here m = elem_mass (already L·ρ_lin) per element.
# ---------------------------------------------------------------------------
@kernel function build_surrogate_loads_kernel!(f_elem, Fx, Fy, elem_L, elem_m, g, azimuth)
    e, s = @index(Global, NTuple)
    @inbounds begin
        TF = eltype(f_elem)
        L = elem_L[e]; m = elem_m[e]
        ca = cos(azimuth[s]); sa = sin(azimuth[s])
        f_elem[e, 1, s] = -g * m * ca
        f_elem[e, 2, s] = Fy[e, s] * L / TF(2) - g * m * sa
        f_elem[e, 3, s] = -Fx[e, s] * L / TF(2)
        f_elem[e, 4, s] = zero(TF)
        f_elem[e, 5, s] = zero(TF)
        f_elem[e, 6, s] = zero(TF)
    end
end

"""
    build_surrogate_loads_gpu!(f_elem, Fx, Fy, elem_L, elem_m, g, azimuth)

Assemble the `(nelem, 6, ns)` surrogate load matrix from the sectional
aero loads and per-element length/mass. Device port of
[`build_surrogate_loads`](@ref). `elem_m[e]` is the element mass
(`mass[1,1]·L`); `azimuth` is length `ns`.
"""
function build_surrogate_loads_gpu!(f_elem, Fx, Fy, elem_L, elem_m, g, azimuth)
    nelem, _, ns = size(f_elem)
    backend = KernelAbstractions.get_backend(f_elem)
    k = build_surrogate_loads_kernel!(backend)
    TF = eltype(f_elem)
    k(f_elem, Fx, Fy, elem_L, elem_m, TF(g), azimuth; ndrange=(nelem, ns))
    KernelAbstractions.synchronize(backend)
    return f_elem
end

# ---------------------------------------------------------------------------
# Kernel 5 — structural→aero feedback (port of update_mesh! →
# interpolate_deflection / interpolate_angle / convert_velocities).
#
# Reads decoded point states u/theta/V (each (3, np, ns)); writes
# delta/def_theta/aerov (each (3, na, ns)). RS is per-sim (length ns).
# ---------------------------------------------------------------------------
@kernel function update_feedback_kernel!(delta, def_theta, aerov,
        u, theta, V, pair1, pair2, percent, ry, rz, RS)
    j, s = @index(Global, NTuple)
    @inbounds begin
        i1 = pair1[j]; i2 = pair2[j]; p = percent[j]
        omp = one(p) - p

        # interpolate_deflection: linear blend of u
        du1 = omp * u[1, i1, s] + p * u[1, i2, s]
        du2 = omp * u[2, i1, s] + p * u[2, i2, s]
        du3 = omp * u[3, i1, s] + p * u[3, i2, s]
        delta[1, j, s] = du1; delta[2, j, s] = du2; delta[3, j, s] = du3

        # interpolate_angle: WMPtoangle each end, then blend
        a1r, a1p, a1y = wmp_to_angle_dev(theta[1, i1, s], theta[2, i1, s], theta[3, i1, s])
        a2r, a2p, a2y = wmp_to_angle_dev(theta[1, i2, s], theta[2, i2, s], theta[3, i2, s])
        def_theta[1, j, s] = omp * a1r + p * a2r
        def_theta[2, j, s] = omp * a1p + p * a2p
        def_theta[3, j, s] = omp * a1y + p * a2y

        # convert_velocities: interp V, transform, remove rotational
        Vs1 = omp * V[1, i1, s] + p * V[1, i2, s]
        Vs2 = omp * V[2, i1, s] + p * V[2, i2, s]
        Vs3 = omp * V[3, i1, s] + p * V[3, i2, s]
        dy = du2; dz = du1
        usx = Vs3; usy = -Vs2; usz = -Vs1
        ray = ry[j] + dy; raz = rz[j] + dz
        omega = RS[s]
        ury = omega * raz
        urz = -omega * ray
        aerov[1, j, s] = usx
        aerov[2, j, s] = usy - ury
        aerov[3, j, s] = usz - urz
    end
end

"""
    update_feedback_gpu!(delta, def_theta, aerov, geo, ip, u, theta, V, RS)

Refresh the structural→aero coupling buffers (`delta/def_theta/aerov`, each
`(3, na, ns)`) from the decoded surrogate point states (`u/theta/V`, each
`(3, np, ns)`). Device port of [`update_mesh!`](@ref) — combines
`interpolate_deflection`, `interpolate_angle` (via [`wmp_to_angle_dev`](@ref)),
and `convert_velocities`. `RS` is the per-sim rotor speed (length `ns`).
"""
function update_feedback_gpu!(delta, def_theta, aerov, geo::AeroGeometryGPU,
                              ip::InterpGPU, u, theta, V, RS)
    _, na, ns = size(delta)
    backend = KernelAbstractions.get_backend(delta)
    k = update_feedback_kernel!(backend)
    k(delta, def_theta, aerov, u, theta, V,
      ip.pair1, ip.pair2, ip.percent, geo.ry, geo.rz, RS; ndrange=(na, ns))
    KernelAbstractions.synchronize(backend)
    return delta, def_theta, aerov
end

# ===========================================================================
# Orchestration: GPU-resident coupled aeroelastic march.
# ===========================================================================

export GPUSurrogateMesh, GPUSurrogateHistory
export initialize_sim_surrogate_gpu, run_sim_surrogate_gpu!

"""
    GPUSurrogateMesh

All device-resident scratch buffers + precomputed constants for the coupled
surrogate march. Built by [`initialize_sim_surrogate_gpu`](@ref). Aero state is
`(na, ns)`; the DS ping-pong pair is `(32, na, ns)`; the surrogate load matrix
is `(nelem, 6, ns)`; the structural feedback is `(3, na, ns)`; and the decoded
point states live in a [`GPUPointStates`](@ref) `(3, np, ns)`.
"""
struct GPUSurrogateMesh{TG, TIP, TB, TD, TBE, TM2, TX, TF3, TM3, TP, TIV, THV}
    geo::TG                    # AeroGeometryGPU
    ip::TIP                    # InterpGPU
    rotor_gpu::RotorGPU
    blade_gpu::TB              # BladeGPU
    dsaf::TD                   # DSAirfoilGPU
    bemt::TBE                  # GPUBEMTOutputs (na, ns)
    Vx::TM2; Vy::TM2; pitch_eff::TM2
    U_ds::TM2; aoa_ds::TM2
    xds_a::TX; xds_b::TX       # (32, na, ns) ping-pong
    Cl::TM2; Cd::TM2; Cm::TM2
    Fx::TM2; Fy::TM2; Mx::TM2
    f_elem::TF3                # (nelem, 6, ns)
    delta::TM3; def_theta::TM3; aerov::TM3   # (3, na, ns)
    decoded::TP                # GPUPointStates (3, np, ns)
    elem_L::TIV; elem_m::TIV   # (nelem,)
    turbine::Bool
    # host inflow scratch (length ns) + device mirrors
    az_host::Vector{Float64}
    Ux::THV; Uy::THV; Uz::THV; RSv::THV; azv::THV
end

"""
    GPUSurrogateHistory

Device history arrays copied back to the host at the end of the march. Aero
loads are `(na, ns, nt)`; the structural point-state history is `(3, np, ns, nt)`.

**Fields**
- `Fx, Fy::(na,ns,nt)` — spanwise aero loads.
- `u, F, M::(3,np,ns,nt)` — decoded displacement / internal force / moment.
"""
struct GPUSurrogateHistory{TA3, TA4}
    Fx::TA3
    Fy::TA3
    u::TA4
    F::TA4
    M::TA4
end

"""
    initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, ns;
                                 ArrayType=Array{Float64}, n_alpha_bemt=721, n_alpha_ds=1441)
        -> (mesh::GPUSurrogateMesh, history::GPUSurrogateHistory)

Allocate every device buffer and precompute every constant the coupled
surrogate march needs, for a batch of `ns` sims on `ArrayType`. Mirrors
[`initialize_sim_surrogate`](@ref) but device-resident and batched.
"""
function initialize_sim_surrogate_gpu(blade::Blade, rotor::Rotor, assembly::GXBeam.Assembly,
                                      env::Environment, tvec, ns::Integer;
                                      ArrayType::Type=Array{Float64},
                                      n_alpha_bemt::Integer=721, n_alpha_ds::Integer=1441)
    TF = eltype(ArrayType)
    na = length(blade.r)
    np = length(assembly.points)
    nelem = length(assembly.elements)
    nt = length(tvec)

    geo  = AeroGeometryGPU(blade, rotor, env; ArrayType=ArrayType)
    ips  = create_interpolationpoints(assembly, blade)
    ip   = InterpGPU(ips; ArrayType=ArrayType)
    rgpu = RotorGPU(rotor)
    bgpu = BladeGPU(blade; n_alpha=n_alpha_bemt, ArrayType=ArrayType)
    dsaf = DSAirfoilGPU(blade; n_alpha=n_alpha_ds, ArrayType=ArrayType)

    bemt = GPUBEMTOutputs(na, ns; ArrayType=ArrayType)
    m2() = to_backend_matrix(ArrayType, zeros(TF, na, ns))
    m3() = to_backend_array(ArrayType, zeros(TF, 3, na, ns))
    xds() = to_backend_array(ArrayType, zeros(TF, DS_NSTATES, na, ns))
    f_elem = to_backend_array(ArrayType, zeros(TF, nelem, 6, ns))

    elem_L = to_backend_vector(ArrayType, TF[assembly.elements[j].L for j in 1:nelem])
    elem_m = to_backend_vector(ArrayType, TF[assembly.elements[j].mass[1,1]*assembly.elements[j].L for j in 1:nelem])

    decoded = GPUPointStates(np, ns; ArrayType=ArrayType)

    hv() = to_backend_vector(ArrayType, zeros(TF, ns))

    mesh = GPUSurrogateMesh(
        geo, ip, rgpu, bgpu, dsaf, bemt,
        m2(), m2(), m2(), m2(), m2(),
        xds(), xds(),
        m2(), m2(), m2(), m2(), m2(), m2(),
        f_elem, m3(), m3(), m3(), decoded,
        elem_L, elem_m, rotor.turbine,
        zeros(Float64, ns), hv(), hv(), hv(), hv(), hv())

    a3() = to_backend_array(ArrayType, zeros(TF, na, ns, nt))
    a4() = to_backend_array(ArrayType, zeros(TF, 3, np, ns, nt))
    history = GPUSurrogateHistory(a3(), a3(), a4(), a4(), a4())

    return mesh, history
end

# Evaluate the (shared) environment's time-callables onto the host inflow
# vectors, then copy to the device mirrors. Per-sim envs could fill these
# columns individually; v1 broadcasts a single env across all sims.
function _upload_inflow!(mesh::GPUSurrogateMesh, env::Environment, t)
    U = env.U(t)
    rs = env.RS(t)
    fill!(mesh.Ux, eltype(mesh.Ux)(U[1]))
    fill!(mesh.Uy, eltype(mesh.Uy)(U[2]))
    fill!(mesh.Uz, eltype(mesh.Uz)(U[3]))
    fill!(mesh.RSv, eltype(mesh.RSv)(rs))
    copyto!(mesh.azv, eltype(mesh.azv).(mesh.az_host))
    return rs
end

# One full aero substep on the device: velocity → BEMT → DS inputs →
# DS advance (prev→cur) → loads/dimensionalize. Does NOT swap the ping-pong.
function _aero_substep!(mesh::GPUSurrogateMesh, env::Environment, xds_cur, xds_prev, pitch, dt; init::Bool)
    aero_velocities_gpu!(mesh.Vx, mesh.Vy, mesh.pitch_eff, mesh.geo,
        mesh.delta, mesh.def_theta, mesh.aerov,
        mesh.Ux, mesh.Uy, mesh.Uz, mesh.RSv, mesh.azv, pitch)
    solve_BEMT_gpu!(mesh.bemt, mesh.rotor_gpu, mesh.blade_gpu, env,
        mesh.Vx, mesh.Vy, mesh.pitch_eff)
    ds_inputs_gpu!(mesh.U_ds, mesh.aoa_ds, mesh.geo, mesh.bemt.W, mesh.bemt.phi, pitch, mesh.turbine)
    if init
        ds_init_step_gpu!(xds_cur, mesh.Cl, mesh.Cd, mesh.Cm, mesh.dsaf, mesh.U_ds, mesh.aoa_ds)
    else
        ds_step_gpu!(xds_cur, xds_prev, mesh.Cl, mesh.Cd, mesh.Cm, mesh.dsaf, mesh.U_ds, mesh.aoa_ds, dt)
    end
    ds_loads_dimensionalize_gpu!(mesh.Fx, mesh.Fy, mesh.Mx, mesh.geo,
        mesh.Cl, mesh.Cd, mesh.Cm, mesh.bemt.phi, mesh.bemt.W, env.rho)
    return nothing
end

_record!(history, mesh, i) = begin
    @views history.Fx[:, :, i] .= mesh.Fx
    @views history.Fy[:, :, i] .= mesh.Fy
    @views history.u[:, :, :, i] .= mesh.decoded.u
    @views history.F[:, :, :, i] .= mesh.decoded.F
    @views history.M[:, :, :, i] .= mesh.decoded.M
    nothing
end

"""
    run_sim_surrogate_gpu!(mesh, history, env, tvec, surr;
                           u0=nothing, pitch=0.0, g=9.81, azimuth0=0.0, verbose=false)

Run the GPU-resident coupled aeroelastic march for all `ns` sims. `surr` is any
batched structural surrogate implementing `encode_initial(surr, ::GPUPointStates)`,
`step_latent(surr, z, f_elem)`, and `decode!(surr, z, ::GPUPointStates)` on
device arrays (see [`BatchedKoopman`]). `u0` is the per-sim IC seed
([`GPUPointStates`](@ref)); `nothing` uses an at-rest zero state.

Every per-step quantity stays on the device; only the scalar time / azimuth /
inflow cross to the host each step. History is written into `history` in place;
copy it back with `Array(history.Fx)` etc. after the call.
"""
function run_sim_surrogate_gpu!(mesh::GPUSurrogateMesh, history::GPUSurrogateHistory,
                                env::Environment, tvec, surr;
                                u0::Union{Nothing,GPUPointStates}=nothing,
                                pitch=0.0, g=9.81, azimuth0=0.0, verbose::Bool=false)
    nt = length(tvec)
    ns = length(mesh.az_host)

    # --- Initial condition (step 1) ---
    t0 = tvec[1]
    fill!(mesh.az_host, azimuth0)
    _upload_inflow!(mesh, env, t0)

    # feedback buffers start at zero (rigid IC), matching the CPU surrogate path.
    xds_prev = mesh.xds_a
    xds_cur  = mesh.xds_b
    _aero_substep!(mesh, env, xds_prev, xds_prev, pitch, tvec[2]-tvec[1]; init=true)

    # encode IC seed → z, decode → history[1], refresh feedback.
    u0_seed = u0 === nothing ? mesh.decoded : u0
    if u0 === nothing
        # zero the seed buffers (decoded currently holds junk); encode a rest state.
        fill!(mesh.decoded.u, 0); fill!(mesh.decoded.theta, 0); fill!(mesh.decoded.V, 0)
        fill!(mesh.decoded.Omega, 0); fill!(mesh.decoded.F, 0); fill!(mesh.decoded.M, 0)
    end
    z = encode_initial(surr, u0_seed)
    decode!(surr, z, mesh.decoded)
    update_feedback_gpu!(mesh.delta, mesh.def_theta, mesh.aerov, mesh.geo, mesh.ip,
        mesh.decoded.u, mesh.decoded.theta, mesh.decoded.V, mesh.RSv)
    _record!(history, mesh, 1)

    verbose && println("GPU surrogate IC done (ns=$ns).")

    # --- Main loop ---
    for i in 2:nt
        t = tvec[i]; dt = t - tvec[i-1]
        dt < 0 && error("run_sim_surrogate_gpu!: negative time step at i=$i")

        # advance azimuth per sim on host, then upload inflow.
        rs = env.RS(t)
        @inbounds for s in 1:ns
            mesh.az_host[s] += rs * dt
        end
        _upload_inflow!(mesh, env, t)

        _aero_substep!(mesh, env, xds_cur, xds_prev, pitch, dt; init=false)

        build_surrogate_loads_gpu!(mesh.f_elem, mesh.Fx, mesh.Fy, mesh.elem_L, mesh.elem_m, g, mesh.azv)
        z = step_latent(surr, z, mesh.f_elem)
        decode!(surr, z, mesh.decoded)
        update_feedback_gpu!(mesh.delta, mesh.def_theta, mesh.aerov, mesh.geo, mesh.ip,
            mesh.decoded.u, mesh.decoded.theta, mesh.decoded.V, mesh.RSv)

        _record!(history, mesh, i)

        xds_cur, xds_prev = xds_prev, xds_cur    # ping-pong swap

        verbose && (i % 100 == 0) && println("  step $i / $nt   t=$t")
    end

    return history
end
