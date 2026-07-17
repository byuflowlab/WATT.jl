#=
GPU-batched BEMT for WATT.jl.

Ports the core of `CCBlade.residual_and_outputs` to a device-agnostic
`KernelAbstractions` kernel that solves all (section, simulation) tuples in
parallel with a fixed-iteration Brent root-find.

Scope of this file (v1):
- Wind-turbine convention only (`turbine == true`).
- Normal operation: Vx > 0, Vy > 0 (root in q1 = (eps, pi/2)). Verified 100%
  hit rate across the operating envelope in examples/gpu_bemt_bracket_check.jl.
- Prandtl tip and tip+hub loss corrections. No Mach/Re/rotation corrections
  (constructor errors if any are set on the Rotor).
- Airfoil tables pre-sampled from the CPU-side DS.Airfoil splines onto a
  uniform alpha grid; on-device lookup is linear interpolation.
- Forward pass only. AD adjoint is planned as v2.

Adam Cardoza
=#

using KernelAbstractions
using Adapt

export BladeGPU, RotorGPU, GPUBEMTOutputs, solve_BEMT_gpu!, N_BRENT_ITERS_DEFAULT

# Default matches the value we plan to validate in examples/gpu_bemt_iters.jl.
const N_BRENT_ITERS_DEFAULT = 20

@enum TipMode::Int32 TIP_NONE = 0 TIP_PRANDTL = 1 TIP_PRANDTL_HUB = 2

# ---------------------------------------------------------------------------
# BladeGPU
# ---------------------------------------------------------------------------

"""
    BladeGPU{TF, TA<:AbstractVector, TM<:AbstractMatrix}

GPU-friendly blade representation. Holds per-section scalars and airfoil
tables pre-sampled onto a uniform alpha grid. Move between hosts and devices
with `Adapt.adapt(ArrayType, blade_gpu)`.

**Fields**
- `r::TA`         — radial station of each aero node.
- `c::TA`         — chord at each aero node.
- `twist::TA`     — twist at each aero node (rad).
- `rhub::TF`      — hub radius.
- `rtip::TF`      — tip radius.
- `alpha_min::TF` — lower bound of the uniform alpha grid (rad).
- `dalpha::TF`    — spacing of the uniform alpha grid (rad).
- `n_alpha::Int32`— number of samples in the uniform alpha grid.
- `cl_table::TM`  — `(n_alpha, n_sections)` sampled lift coefficient.
- `cd_table::TM`  — `(n_alpha, n_sections)` sampled drag coefficient.
"""
struct BladeGPU{TF, TA<:AbstractVector, TM<:AbstractMatrix}
    r::TA
    c::TA
    twist::TA
    rhub::TF
    rtip::TF
    alpha_min::TF
    dalpha::TF
    n_alpha::Int32
    cl_table::TM
    cd_table::TM
end

"""
    BladeGPU(blade::Blade; n_alpha=721, ArrayType=Array{Float64})

Build a `BladeGPU` from a `Blade`. Evaluates each section's `cl`/`cd` on a
uniform alpha grid spanning `[-π, π]` and packs the samples into arrays of
the requested `ArrayType` (e.g. `Array{Float64}` for CPU, `CuArray{Float32}`
for CUDA, `MtlArray{Float32}` for Metal).
"""
function BladeGPU(blade::Blade; n_alpha::Integer=721, ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    n = length(blade.r)

    alpha_min = TF(-pi)
    alpha_max = TF(pi)
    dalpha = (alpha_max - alpha_min) / TF(n_alpha - 1)

    cl_h = Array{TF}(undef, n_alpha, n)
    cd_h = Array{TF}(undef, n_alpha, n)
    for j in 1:n
        af = blade.airfoils[j]
        for i in 1:n_alpha
            α = alpha_min + dalpha * TF(i - 1)
            cl_h[i, j] = TF(af.cl(α))
            cd_h[i, j] = TF(af.cd(α))
        end
    end

    r_v     = to_backend_vector(ArrayType, TF.(blade.r))
    c_v     = to_backend_vector(ArrayType, TF.(blade.c))
    twist_v = to_backend_vector(ArrayType, TF.(blade.twist))
    cl_m    = to_backend_matrix(ArrayType, cl_h)
    cd_m    = to_backend_matrix(ArrayType, cd_h)

    return BladeGPU(r_v, c_v, twist_v,
                    TF(blade.rhub), TF(blade.rtip),
                    alpha_min, dalpha, Int32(n_alpha),
                    cl_m, cd_m)
end

# ArrayType-driven host→backend transfer. Users of GPU backends add methods
# for their vendor's array (e.g. `to_backend_vector(::Type{CuArray{T}}, v)`
# in a package extension); the fallback handles plain `Array` and any GPU
# array constructible from a host array.
to_backend_vector(::Type{Array{T}}, v::AbstractVector) where {T} = Vector{T}(v)
to_backend_matrix(::Type{Array{T}}, m::AbstractMatrix) where {T} = Matrix{T}(m)
to_backend_vector(AT::Type, v::AbstractVector) = AT(v)
to_backend_matrix(AT::Type, m::AbstractMatrix) = AT(m)

Adapt.adapt_structure(to, b::BladeGPU) = BladeGPU(
    adapt(to, b.r),
    adapt(to, b.c),
    adapt(to, b.twist),
    b.rhub, b.rtip,
    b.alpha_min, b.dalpha, b.n_alpha,
    adapt(to, b.cl_table),
    adapt(to, b.cd_table),
)

# Fallback inner constructor used by adapt above (all field types inferred).
BladeGPU(r, c, twist, rhub, rtip, amin, dalpha, nα, cl, cd) =
    BladeGPU{typeof(rhub), typeof(r), typeof(cl)}(
        r, c, twist, rhub, rtip, amin, dalpha, nα, cl, cd)

# ---------------------------------------------------------------------------
# RotorGPU
# ---------------------------------------------------------------------------

"""
    RotorGPU

Minimal, kernel-safe rotor description. Holds only the scalar fields the
device kernel needs. Correction models are restricted to Prandtl tip / tip+hub
in v1.

**Fields**
- `B::Int32`     — number of blades.
- `turbine::Bool`
- `tip_mode::TipMode` — `TIP_NONE`, `TIP_PRANDTL`, or `TIP_PRANDTL_HUB`.
"""
struct RotorGPU
    B::Int32
    turbine::Bool
    tip_mode::TipMode
end

"""
    RotorGPU(rotor::Rotor)

Build a `RotorGPU` from a `Rotor`. Errors if any unsupported correction
model (Mach, Re, rotation) is set. Recognizes `CCBlade.PrandtlTip` and
`CCBlade.PrandtlTipHub`; other tip corrections are also rejected.
"""
function RotorGPU(rotor::Rotor)
    !rotor.turbine && error("RotorGPU: v1 requires turbine=true (wind-turbine convention).")
    !isnothing(rotor.mach)     && error("RotorGPU: Mach correction not supported in v1.")
    !isnothing(rotor.re)       && error("RotorGPU: Re correction not supported in v1.")
    !isnothing(rotor.rotation) && error("RotorGPU: rotation correction not supported in v1.")

    tip_mode = if isnothing(rotor.tip)
        TIP_NONE
    elseif rotor.tip isa CCBlade.PrandtlTip
        TIP_PRANDTL
    elseif rotor.tip isa CCBlade.PrandtlTipHub
        TIP_PRANDTL_HUB
    else
        error("RotorGPU: unsupported tip correction $(typeof(rotor.tip)). v1 supports PrandtlTip and PrandtlTipHub only.")
    end

    return RotorGPU(Int32(rotor.B), rotor.turbine, tip_mode)
end

# ---------------------------------------------------------------------------
# Device-side helpers
# ---------------------------------------------------------------------------

# Linear interpolation on a uniform alpha grid. Clamps to [alpha_min, alpha_max].
@inline function interp_af(alpha, section, cl_table, cd_table,
                            alpha_min, dalpha, n_alpha)
    # locate bracket
    x = (alpha - alpha_min) / dalpha
    # clamp to valid range: index 0 .. n_alpha - 1 (0-based)
    x0 = clamp(x, zero(x), oftype(x, n_alpha - 1))
    i0 = unsafe_trunc(Int32, x0)             # 0-based lower index
    if i0 >= n_alpha - 1
        i0 = Int32(n_alpha - 2)
    end
    t = x0 - oftype(x0, i0)
    # 1-based indexing into the tables (row = alpha, col = section)
    r0 = i0 + Int32(1)
    r1 = r0 + Int32(1)
    cl = (one(t) - t) * cl_table[r0, section] + t * cl_table[r1, section]
    cd = (one(t) - t) * cd_table[r0, section] + t * cd_table[r1, section]
    return cl, cd
end

# Prandtl tip / tip+hub correction. Branchless once `tip_mode` is inlined.
@inline function tip_loss(tip_mode::TipMode, r, Rhub, Rtip, phi, B)
    if tip_mode == TIP_NONE
        return one(r)
    end
    asphi = abs(sin(phi))
    # Tip
    ftip = (B / oftype(r, 2)) * (Rtip / r - one(r)) / asphi
    Ftip = (oftype(r, 2) / oftype(r, pi)) * acos(exp(-ftip))
    if tip_mode == TIP_PRANDTL
        return Ftip
    end
    # Tip + hub
    fhub = (B / oftype(r, 2)) * (r / Rhub - one(r)) / asphi
    Fhub = (oftype(r, 2) / oftype(r, pi)) * acos(exp(-fhub))
    return Ftip * Fhub
end

"""
Device-side BEMT residual (port of `CCBlade.residual_and_outputs`).
Wind-turbine convention, no Mach/Re/rotation corrections.

Returns `(R, Np, Tp, alpha, W, cl_signed, cd, cn, ct, F, a, ap)` where
`cl_signed` is the sign-flipped cl consistent with the turbine convention
in CCBlade.
"""
@inline function bemt_residual_and_outputs(phi, section, Vx, Vy, pitch,
        r, chord, twist, Rhub, Rtip, rho,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        B_blades, tip_mode)

    TF = typeof(phi)
    sphi, cphi = sincos(phi)

    # angle of attack (wind-turbine convention flips sign for airfoil lookup).
    # xv[3] in the CPU code is `twist`; CCBlade's residual then does
    # alpha = (twist + pitch) - phi, so we do the same.
    alpha = (twist + pitch) - phi

    # cl, cd (turbine convention: sample at -alpha, then flip cl sign)
    cl_raw, cd = interp_af(-alpha, section, cl_table, cd_table,
                            alpha_min, dalpha, n_alpha)
    cl = -cl_raw

    # resolve into normal / tangential coefficients
    cn = cl * cphi - cd * sphi
    ct = cl * sphi + cd * cphi

    # hub/tip loss
    F = tip_loss(tip_mode, r, Rhub, Rtip, phi, TF(B_blades))

    sigma_p = TF(B_blades) * chord / (TF(2) * TF(pi) * r)
    k  = cn * sigma_p / (TF(4) * F * sphi * sphi)
    kp = ct * sigma_p / (TF(4) * F * sphi * cphi)

    # ----- q1-only path (Vx > 0, Vy > 0, phi > 0) -----
    # (Deliberately no Vx=0 / Vy=0 / phi<0 branches — those are outside our
    # supported regime; see gpu_bemt_bracket_check.jl.)

    # k regime
    a = if k >= -TF(2) / TF(3)   # momentum
        k / (one(TF) - k)
    else                          # empirical (Buhl variant)
        g1 = TF(2) * k + TF(1) / TF(9)
        g2 = -TF(2) * k - TF(1) / TF(3)
        g3 = -TF(2) * k - TF(7) / TF(9)
        (g1 + sqrt(g2)) / g3
    end
    u = a * Vx

    ap = kp / (one(TF) + kp)
    v = ap * Vy

    # residual
    R = sin(phi) / (one(TF) + a) - (Vx / Vy) * cos(phi) / (one(TF) - ap)

    # loads (dimensional per-unit-span; hub/tip losses applied via cn/ct through F above)
    W = sqrt((Vx + u) * (Vx + u) + (Vy - v) * (Vy - v))
    Np = cn * TF(0.5) * rho * W * W * chord
    Tp = ct * TF(0.5) * rho * W * W * chord

    # turbine convention on outputs (matches CCBlade branch)
    return (R,
            -Np, -Tp,
            -alpha, W,
            -cl, cd, -cn, -ct,
            F,
            -a, -ap)
end

# ---------------------------------------------------------------------------
# Fixed-iteration Brent
# ---------------------------------------------------------------------------

# Fixed sub-bracket scan over (phi_a, phi_b). Evaluates the residual at
# `NSUB + 1` uniformly-spaced phi values and returns the first sub-interval
# where the residual changes sign, or the full [phi_a, phi_b] with
# `found = false` if no sign change is seen. Every thread executes the same
# number of residual evaluations, keeping the warp in lockstep.
@inline function scan_sign_change(resfun, phi_a, phi_b, ::Val{NSUB}) where {NSUB}
    TF = typeof(phi_a)
    dphi = (phi_b - phi_a) / TF(NSUB)
    Rprev = resfun(phi_a)
    found = false
    phi_lo = phi_a
    phi_hi = phi_b
    for i in 1:NSUB
        phi_i = phi_a + dphi * TF(i)
        Rnext = resfun(phi_i)
        if !found && (Rprev * Rnext) < zero(TF)
            phi_lo = phi_i - dphi
            phi_hi = phi_i
            found = true
        end
        Rprev = Rnext
    end
    return phi_lo, phi_hi, found
end

# Brent's method as in FLOWMath but with a fixed iteration count and no
# early exit. Every thread does exactly `n_iters` iterations. `f` is a
# closure that captures the (section, sim) inputs and returns just the
# residual R (not the full outputs).
@inline function brent_fixed(f, xa, xb, ::Val{N}) where {N}
    TF = typeof(xa)
    fa = f(xa)
    fb = f(xb)

    # If the caller violates the bracket assumption, return xb (safe default).
    # Threads that pass through here still burn N iterations to keep the
    # warp lockstep.
    if fa * fb > zero(TF)
        return xb, fb
    end

    if abs(fa) < abs(fb)
        xa, xb = xb, xa
        fa, fb = fb, fa
    end

    xc = xa
    fc = fa
    mflag = true
    xd = xa   # unused on first iter; initialize to something typed
    tol = TF(2) * eps(TF) * abs(xb) + TF(1e-12)

    for _ in 1:N
        # interpolation
        s = if fa != fc && fb != fc
            # inverse quadratic
            xa * fb * fc / ((fa - fb) * (fa - fc)) +
            xb * fa * fc / ((fb - fa) * (fb - fc)) +
            xc * fa * fb / ((fc - fa) * (fc - fb))
        else
            # secant
            xb - fb * (xb - xa) / (fb - fa)
        end

        # bisection guard
        cond1 = !((s > (TF(3) * xa + xb) / TF(4) && s < xb) ||
                  (s < (TF(3) * xa + xb) / TF(4) && s > xb))
        cond2 = mflag  && abs(s - xb) >= abs(xb - xc) / TF(2)
        cond3 = !mflag && abs(s - xb) >= abs(xc - xd) / TF(2)
        cond4 = mflag  && abs(xb - xc) < tol
        cond5 = !mflag && abs(xc - xd) < tol

        if cond1 || cond2 || cond3 || cond4 || cond5
            s = (xa + xb) / TF(2)
            mflag = true
        else
            mflag = false
        end

        fs = f(s)
        xd = xc
        xc = xb
        fc = fb

        if fa * fs < zero(TF)
            xb = s
            fb = fs
        else
            xa = s
            fa = fs
        end

        if abs(fa) < abs(fb)
            xa, xb = xb, xa
            fa, fb = fb, fa
        end
    end

    return xb, fb
end

# ---------------------------------------------------------------------------
# Kernel and public entry point
# ---------------------------------------------------------------------------

# Output layout: pre-allocated matrices of shape (n_sections, n_sims).
"""
    GPUBEMTOutputs

Container of pre-allocated `(n_sections, n_sims)` device matrices written by
[`solve_BEMT_gpu!`](@ref). Constructing this once and reusing it every step
avoids per-step device allocations.
"""
struct GPUBEMTOutputs{TM<:AbstractMatrix, TB<:AbstractMatrix}
    Np::TM
    Tp::TM
    phi::TM
    alpha::TM
    W::TM
    cl::TM
    cd::TM
    cn::TM
    ct::TM
    F::TM
    a::TM
    ap::TM
    success::TB   # Bool matrix; true if Brent's bracket assumption held
end

"""
    GPUBEMTOutputs(n_sections, n_sims; ArrayType=Array{Float64})

Allocate output matrices for [`solve_BEMT_gpu!`](@ref). `ArrayType` matches
the array flavor of the `BladeGPU` you plan to solve against.
"""
function GPUBEMTOutputs(n_sections::Integer, n_sims::Integer; ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    shape = (n_sections, n_sims)
    make() = to_backend_matrix(ArrayType, zeros(TF, shape))
    return GPUBEMTOutputs(
        make(), make(), make(), make(), make(),
        make(), make(), make(), make(), make(),
        make(), make(),
        to_backend_matrix(similar_type(ArrayType, Bool), zeros(Bool, shape)),
    )
end

# Swap the element type of an array-type constructor (e.g. CuArray{Float64} →
# CuArray{Int32}) without hard-coding a package. The generic fallback rebuilds
# the type from its wrapper, so CuArray / ROCArray / MtlArray / Array all work
# out of the box; a vendor package may still add a more specific method.
similar_type(::Type{Array{T}}, ::Type{U}) where {T, U} = Array{U}
similar_type(AT::Type{<:AbstractArray}, ::Type{U}) where {U} = Base.typename(AT).wrapper{U}
similar_type(AT::Type, ::Type{U}) where {U} =
    error("similar_type: don't know how to build a $U-eltype array like $AT. " *
          "Define `WATT.similar_type(::Type{$AT}, ::Type{$U})`.")

# Per-thread BEMT solve for one (section, sim): bracket → fixed-iteration
# Brent → recompute outputs. Shared by the vector-pitch and matrix-pitch
# kernels so the physics lives in exactly one place. Returns the full output
# tuple plus the bracket-valid flag.
@inline function bemt_solve_one(Vx, Vy, pitch, section::Int32,
        r, chord, twist, rhub, rtip, rho,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        B_blades, tip_mode, ::Val{NIT}) where {NIT}

    # Sub-bracket scans on q1 = (eps, pi/2), then falling back to
    # q3 = (pi/2, pi - eps) only if q1 doesn't produce a sign-change.
    # Rationale: the aerostructural-bracket trace observed 0.25% of
    # transient solves landing in q3, all at root-cylinder sections where
    # Vy is small (phi_geom sits right at 90°). Because those sections
    # cluster into a single warp on the GPU, the divergence cost of the
    # conditional q3 is small — most warps never take the fallback.
    TF = typeof(Vx)
    eps_phi = TF(1e-6)

    # residual closure captures the inputs for this thread
    resfun = phi -> bemt_residual_and_outputs(
        phi, section, Vx, Vy, pitch,
        r, chord, twist, rhub, rtip,
        rho,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        B_blades, tip_mode)[1]

    # --- Try q1 first ---
    phi_lo, phi_hi, q1_found = scan_sign_change(resfun,
        eps_phi, TF(pi / 2) - eps_phi, Val(10))

    # --- Fallback to q3 only if q1 didn't produce a sign change ---
    q3_found = false
    if !q1_found
        phi_lo, phi_hi, q3_found = scan_sign_change(resfun,
            TF(pi / 2) + eps_phi, TF(pi) - eps_phi, Val(10))
    end
    bracket_ok = q1_found | q3_found

    phistar, _ = brent_fixed(resfun, phi_lo, phi_hi, Val(NIT))

    # Recompute all outputs at phistar. Even if bracket_ok is false, we still
    # run the residual once — cheaper than branching. If bracket was invalid,
    # the caller clobbers all outputs with zero (matches CPU CCBlade.Outputs()).
    (_, Np, Tp, alpha, W, cl, cd, cn, ct, F, a, ap) =
        bemt_residual_and_outputs(
            phistar, section, Vx, Vy, pitch,
            r, chord, twist, rhub, rtip,
            rho,
            alpha_min, dalpha, n_alpha, cl_table, cd_table,
            B_blades, tip_mode)

    return (Np, Tp, phistar, alpha, W, cl, cd, cn, ct, F, a, ap, bracket_ok)
end

# Write one thread's outputs, zeroing everything if the bracket was invalid.
@inline function _bemt_write_outputs!(out_Np, out_Tp, out_phi, out_alpha, out_W,
        out_cl, out_cd, out_cn, out_ct, out_F, out_a, out_ap, out_success,
        j, s, res)
    (Np, Tp, phistar, alpha, W, cl, cd, cn, ct, F, a, ap, bracket_ok) = res
    z = zero(typeof(Np))
    out_Np[j, s]      = ifelse(bracket_ok, Np,      z)
    out_Tp[j, s]      = ifelse(bracket_ok, Tp,      z)
    out_phi[j, s]     = ifelse(bracket_ok, phistar, z)
    out_alpha[j, s]   = ifelse(bracket_ok, alpha,   z)
    out_W[j, s]       = ifelse(bracket_ok, W,       z)
    out_cl[j, s]      = ifelse(bracket_ok, cl,      z)
    out_cd[j, s]      = ifelse(bracket_ok, cd,      z)
    out_cn[j, s]      = ifelse(bracket_ok, cn,      z)
    out_ct[j, s]      = ifelse(bracket_ok, ct,      z)
    out_F[j, s]       = ifelse(bracket_ok, F,       z)
    out_a[j, s]       = ifelse(bracket_ok, a,       z)
    out_ap[j, s]      = ifelse(bracket_ok, ap,      z)
    out_success[j, s] = bracket_ok
end

@kernel function bemt_solve_kernel!(
        out_Np, out_Tp, out_phi, out_alpha, out_W,
        out_cl, out_cd, out_cn, out_ct, out_F,
        out_a, out_ap, out_success,
        Vx_all, Vy_all, pitch_all,
        r_vec, c_vec, twist_vec, rhub, rtip,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        rho,
        B_blades, tip_mode,
        ::Val{NIT}) where {NIT}

    j, s = @index(Global, NTuple)   # j = section index, s = sim index

    res = bemt_solve_one(Vx_all[j, s], Vy_all[j, s], pitch_all[s], Int32(j),
        r_vec[j], c_vec[j], twist_vec[j], rhub, rtip, rho,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        B_blades, tip_mode, Val(NIT))

    _bemt_write_outputs!(out_Np, out_Tp, out_phi, out_alpha, out_W,
        out_cl, out_cd, out_cn, out_ct, out_F, out_a, out_ap, out_success, j, s, res)
end

# Matrix-pitch variant: pitch varies per (section, sim). Used by the coupled
# aeroelastic solver, where the twist correction `pitch - def_theta_x` differs
# per section and per sim.
@kernel function bemt_solve_kernel_pitchmat!(
        out_Np, out_Tp, out_phi, out_alpha, out_W,
        out_cl, out_cd, out_cn, out_ct, out_F,
        out_a, out_ap, out_success,
        Vx_all, Vy_all, pitch_all,
        r_vec, c_vec, twist_vec, rhub, rtip,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        rho,
        B_blades, tip_mode,
        ::Val{NIT}) where {NIT}

    j, s = @index(Global, NTuple)

    res = bemt_solve_one(Vx_all[j, s], Vy_all[j, s], pitch_all[j, s], Int32(j),
        r_vec[j], c_vec[j], twist_vec[j], rhub, rtip, rho,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        B_blades, tip_mode, Val(NIT))

    _bemt_write_outputs!(out_Np, out_Tp, out_phi, out_alpha, out_W,
        out_cl, out_cd, out_cn, out_ct, out_F, out_a, out_ap, out_success, j, s, res)
end

"""
    solve_BEMT_gpu!(outputs, rotor_gpu, blade_gpu, env, Vx, Vy, pitch;
                    n_iters=N_BRENT_ITERS_DEFAULT)

Solve the BEMT for every `(section, sim)` pair in parallel. Writes results
into `outputs` in place. `Vx`, `Vy` are `(n_sections, n_sims)` matrices;
`pitch` is a length-`n_sims` vector. All array arguments must live on the
same backend (as determined by `KernelAbstractions.get_backend`).

`env` is any `Environment` — only `rho`, `mu`, `a` are read (as scalars).
"""
function solve_BEMT_gpu!(outputs::GPUBEMTOutputs,
                          rotor_gpu::RotorGPU,
                          blade_gpu::BladeGPU,
                          env::Environment,
                          Vx::AbstractMatrix, Vy::AbstractMatrix,
                          pitch::AbstractVector;
                          n_iters::Integer=N_BRENT_ITERS_DEFAULT)
    _, n_sims = size(Vx)
    length(pitch) == n_sims ||
        error("solve_BEMT_gpu!: length(pitch) $(length(pitch)) != n_sims $n_sims")
    backend = KernelAbstractions.get_backend(Vx)
    _solve_BEMT_gpu!(bemt_solve_kernel!(backend), outputs, rotor_gpu, blade_gpu,
                     env, Vx, Vy, pitch; n_iters)
    return outputs
end

"""
    solve_BEMT_gpu!(outputs, rotor_gpu, blade_gpu, env, Vx, Vy, pitch::AbstractMatrix; n_iters)

Matrix-pitch overload: `pitch` is a `(n_sections, n_sims)` matrix so the blade
pitch can vary per section and per sim. Used by the coupled aeroelastic solver,
where the effective pitch `pitch - def_theta_x` carries the structural twist
correction at each station. Otherwise identical to the vector-pitch method.
"""
function solve_BEMT_gpu!(outputs::GPUBEMTOutputs,
                          rotor_gpu::RotorGPU,
                          blade_gpu::BladeGPU,
                          env::Environment,
                          Vx::AbstractMatrix, Vy::AbstractMatrix,
                          pitch::AbstractMatrix;
                          n_iters::Integer=N_BRENT_ITERS_DEFAULT)
    n_sections, n_sims = size(Vx)
    size(pitch) == (n_sections, n_sims) ||
        error("solve_BEMT_gpu!: size(pitch) $(size(pitch)) != $((n_sections, n_sims))")
    backend = KernelAbstractions.get_backend(Vx)
    _solve_BEMT_gpu!(bemt_solve_kernel_pitchmat!(backend), outputs, rotor_gpu, blade_gpu,
                     env, Vx, Vy, pitch; n_iters)
    return outputs
end

# Shared launch glue for both pitch flavors. `kernel` is the already-selected
# compiled kernel; `pitch` is a vector or matrix matching it.
function _solve_BEMT_gpu!(kernel, outputs::GPUBEMTOutputs,
                          rotor_gpu::RotorGPU, blade_gpu::BladeGPU, env::Environment,
                          Vx::AbstractMatrix, Vy::AbstractMatrix, pitch;
                          n_iters::Integer=N_BRENT_ITERS_DEFAULT)
    n_sections, n_sims = size(Vx)
    size(Vy) == (n_sections, n_sims) ||
        error("solve_BEMT_gpu!: size(Vy) $(size(Vy)) != size(Vx) $((n_sections, n_sims))")
    size(outputs.Np) == (n_sections, n_sims) ||
        error("solve_BEMT_gpu!: outputs sized $(size(outputs.Np)) != $((n_sections, n_sims))")

    backend = KernelAbstractions.get_backend(Vx)
    TF = eltype(blade_gpu.r)

    kernel(
        outputs.Np, outputs.Tp, outputs.phi, outputs.alpha, outputs.W,
        outputs.cl, outputs.cd, outputs.cn, outputs.ct, outputs.F,
        outputs.a, outputs.ap, outputs.success,
        Vx, Vy, pitch,
        blade_gpu.r, blade_gpu.c, blade_gpu.twist, blade_gpu.rhub, blade_gpu.rtip,
        blade_gpu.alpha_min, blade_gpu.dalpha, blade_gpu.n_alpha,
        blade_gpu.cl_table, blade_gpu.cd_table,
        TF(env.rho),
        rotor_gpu.B, rotor_gpu.tip_mode,
        Val(Int(n_iters));
        ndrange=(n_sections, n_sims),
    )
    KernelAbstractions.synchronize(backend)
    return outputs
end
