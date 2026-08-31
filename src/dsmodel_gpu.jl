#=
GPU-batched dynamic-stall model for WATT.jl.

Ports the Beddoes-Leishman v3 (AeroDyn + Gonzalez) indicial model — the only
DS variant WATT builds via `make_dsairfoil` — to device-agnostic
`KernelAbstractions` kernels that march every (section, simulation) pair in
parallel. The model is a pure explicit 32-state recurrence, so the kernel is a
stateful elementwise map over (section, sim), launched once per timestep.

Scope (v1):
- BeddoesLeishman `Discrete`, version 3, `ADGSP` separation point only. The
  `DSAirfoilGPU` constructor errors on any other configuration.
- Forward pass only (no AD adjoint yet).
- Airfoil cl/cd/cm pre-sampled from the CPU Akima splines onto a uniform alpha
  grid; on-device lookup is linear interpolation (same approach as bemt_gpu.jl).
- Full state history retained on device: `(32, n_sections, n_sims, nt)`.

The `@inline` functions `ds_init_adg`, `ds_step_adg`, `ds_loads_adg` are the
single source of truth: the KernelAbstractions kernels and any CPU reference
loop both call them, so a CPU-vs-GPU comparison is a pure device-correctness
check. Fidelity to the original Akima model is a separate comparison against a
golden trace.

Ports of:
- `update_states_ADG!`  -> `ds_step_adg`
- `BLADG_coefficients`  -> `ds_loads_adg`
- `initialize_ADG`      -> `ds_init_adg`
- `separationpoint(::ADGSP,…)` / `chordwiseseparationpoint(::ADGSP,…)` inlined.
- `blend_cosine` (utils.jl) inlined.

Adam Cardoza
=#

# Exports are centralized in WATT.jl — see the `# GPU backend` block there.

"""
    DS_NSTATES

Number of Beddoes-Leishman dynamic-stall states carried per blade section on the
GPU, currently `32`.

The GPU state array is a dense `(DS_NSTATES, nsections, nblades, ntime)` block
rather than the per-model-variant layout the CPU path uses, so this constant is
what sizes the leading dimension.

**Notes**
Internal — not exported. Reference it as `WATT.DS_NSTATES`, which is how
`examples/gpu_dsmodel_*.jl` already spells it.
"""
const DS_NSTATES = 32
const DS_FCLIMIT = (1.0 + 0.2)^2   # matches DynamicStallModels.fclimit

# ---------------------------------------------------------------------------
# Per-section constant layout. Packed into a single (DS_NCONST, n_sections)
# matrix so the whole constant set moves to the device in one array and the
# kernel reads consts[k, j].
# ---------------------------------------------------------------------------
const C_DCNDALPHA = 1
const C_ALPHA0    = 2
const C_CHORD     = 3
const C_XCP       = 4
const C_ALPHACUT2 = 5
const C_CUTRAD    = 6
const C_A1        = 7
const C_A2        = 8
const C_A5        = 9
const C_B1        = 10
const C_B2        = 11
const C_B5        = 12
const C_TP        = 13
const C_TF0       = 14
const C_TV0       = 15
const C_TVL       = 16
const C_TSH       = 17
const C_ZETA      = 18
const C_CD0       = 19
const C_CM0       = 20
const C_CN1       = 21
const C_ETA       = 22
const C_ASOUND    = 23
const DS_NCONST   = 23

# ---------------------------------------------------------------------------
# DSAirfoilGPU
# ---------------------------------------------------------------------------

"""
    DSAirfoilGPU{TF, TM<:AbstractMatrix}

Device-friendly per-section dynamic-stall airfoil data. Holds the packed
constant matrix and cl/cd/cm sampled on a uniform alpha grid. Move between host
and device with `Adapt.adapt(ArrayType, dsafgpu)`.

**Fields**
- `consts::TM`    — `(DS_NCONST, n_sections)` per-section scalar constants.
- `cl_table::TM`  — `(n_alpha, n_sections)` sampled lift coefficient.
- `cd_table::TM`  — `(n_alpha, n_sections)` sampled drag coefficient.
- `cm_table::TM`  — `(n_alpha, n_sections)` sampled moment coefficient.
- `alpha_min::TF` — lower bound of the uniform alpha grid (rad).
- `dalpha::TF`    — spacing of the uniform alpha grid (rad).
- `n_alpha::Int32`— number of samples in the uniform alpha grid.
"""
struct DSAirfoilGPU{TF, TM<:AbstractMatrix}
    consts::TM
    cl_table::TM
    cd_table::TM
    cm_table::TM
    alpha_min::TF
    dalpha::TF
    n_alpha::Int32
end

"""
    DSAirfoilGPU(blade::Blade; n_alpha=1441, ArrayType=Array{Float64})

Build a `DSAirfoilGPU` from a `Blade`. Validates that every section carries a
`BeddoesLeishman` `Discrete` version-3 model with an `ADGSP` separation point,
packs the per-section constants, and samples cl/cd/cm onto a uniform grid over
`[-π, π]`.
"""
function DSAirfoilGPU(blade::Blade; n_alpha::Integer=1441, ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    n  = length(blade.r)

    consts_h = Array{TF}(undef, DS_NCONST, n)
    for j in 1:n
        af = blade.airfoils[j]
        m  = af.model
        (m isa DS.BeddoesLeishman) ||
            error("DSAirfoilGPU: section $j model is $(typeof(m)); v1 requires BeddoesLeishman.")
        (m.detype isa DS.Discrete) ||
            error("DSAirfoilGPU: section $j is not Discrete; v1 requires the discrete model.")
        (m.version == 3) ||
            error("DSAirfoilGPU: section $j is version $(m.version); v1 supports version 3 (Gonzalez) only.")
        (af.sfun isa DS.ADGSP) ||
            error("DSAirfoilGPU: section $j separation point is $(typeof(af.sfun)); v1 requires ADGSP.")

        consts_h[C_DCNDALPHA, j] = af.dcndalpha
        consts_h[C_ALPHA0, j]    = af.alpha0
        consts_h[C_CHORD, j]     = blade.c[j]
        consts_h[C_XCP, j]       = blade.xcp[j]
        consts_h[C_ALPHACUT2, j] = af.alphacut[2]
        consts_h[C_CUTRAD, j]    = af.cutrad
        consts_h[C_A1, j]        = m.A[1]
        consts_h[C_A2, j]        = m.A[2]
        consts_h[C_A5, j]        = m.A[3]
        consts_h[C_B1, j]        = m.b[1]
        consts_h[C_B2, j]        = m.b[2]
        consts_h[C_B5, j]        = m.b[3]
        consts_h[C_TP, j]        = m.T[1]
        consts_h[C_TF0, j]       = m.T[2]
        consts_h[C_TV0, j]       = m.T[3]
        consts_h[C_TVL, j]       = m.T[4]
        consts_h[C_TSH, j]       = m.T[5]
        consts_h[C_ZETA, j]      = m.zeta
        consts_h[C_CD0, j]       = m.Cd0
        consts_h[C_CM0, j]       = m.Cm0
        consts_h[C_CN1, j]       = m.Cn1
        consts_h[C_ETA, j]       = m.eta
        consts_h[C_ASOUND, j]    = m.a
    end

    alpha_min = TF(-pi)
    alpha_max = TF(pi)
    dalpha = (alpha_max - alpha_min) / TF(n_alpha - 1)

    cl_h = Array{TF}(undef, n_alpha, n)
    cd_h = Array{TF}(undef, n_alpha, n)
    cm_h = Array{TF}(undef, n_alpha, n)
    for j in 1:n
        af = blade.airfoils[j]
        for i in 1:n_alpha
            α = alpha_min + dalpha * TF(i - 1)
            cl_h[i, j] = TF(af.cl(α))
            cd_h[i, j] = TF(af.cd(α))
            cm_h[i, j] = TF(af.cm(α))
        end
    end

    return DSAirfoilGPU(
        to_backend_matrix(ArrayType, consts_h),
        to_backend_matrix(ArrayType, cl_h),
        to_backend_matrix(ArrayType, cd_h),
        to_backend_matrix(ArrayType, cm_h),
        alpha_min, dalpha, Int32(n_alpha),
    )
end

Adapt.adapt_structure(to, a::DSAirfoilGPU) = DSAirfoilGPU(
    adapt(to, a.consts),
    adapt(to, a.cl_table),
    adapt(to, a.cd_table),
    adapt(to, a.cm_table),
    a.alpha_min, a.dalpha, a.n_alpha,
)

# Fallback inner constructor used by adapt (all field types inferred).
DSAirfoilGPU(consts, cl, cd, cm, amin, dalpha, nα) =
    DSAirfoilGPU{typeof(amin), typeof(consts)}(consts, cl, cd, cm, amin, dalpha, nα)

# ---------------------------------------------------------------------------
# DSHistory — preallocated (32, n_sections, n_sims, nt) state history + loads.
# ---------------------------------------------------------------------------

"""
    DSHistory{TA<:AbstractArray}

Preallocated device arrays holding the full DS march.

**Fields**
- `xds::TA`  — `(32, n_sections, n_sims, nt)` state history.
- `Cl::TA`   — `(n_sections, n_sims, nt)` lift coefficient.
- `Cd::TA`   — `(n_sections, n_sims, nt)` drag coefficient.
- `Cm::TA`   — `(n_sections, n_sims, nt)` moment coefficient.
"""
struct DSHistory{T4<:AbstractArray, T3<:AbstractArray}
    xds::T4
    Cl::T3
    Cd::T3
    Cm::T3
end

"""
    DSHistory(n_sections, n_sims, nt; ArrayType=Array{Float64})

Allocate the DS history/loads arrays. `ArrayType` matches the `DSAirfoilGPU`
you plan to march against.
"""
function DSHistory(n_sections::Integer, n_sims::Integer, nt::Integer; ArrayType::Type=Array{Float64})
    TF = eltype(ArrayType)
    xds = to_backend_array(ArrayType, zeros(TF, DS_NSTATES, n_sections, n_sims, nt))
    Cl  = to_backend_array(ArrayType, zeros(TF, n_sections, n_sims, nt))
    Cd  = to_backend_array(ArrayType, zeros(TF, n_sections, n_sims, nt))
    Cm  = to_backend_array(ArrayType, zeros(TF, n_sections, n_sims, nt))
    return DSHistory(xds, Cl, Cd, Cm)
end

# N-dim host->backend transfer (companion to to_backend_vector/matrix in bemt_gpu.jl).
to_backend_array(::Type{Array{T}}, a::AbstractArray) where {T} = Array{T}(a)
to_backend_array(AT::Type, a::AbstractArray) = AT(a)

# ---------------------------------------------------------------------------
# Device-side helpers
# ---------------------------------------------------------------------------

# Linear interpolation of cl/cd/cm on the uniform alpha grid. Shares the bracket
# computation across the three tables. Clamps to [alpha_min, alpha_max].
@inline function ds_af3(alpha, section, cl_table, cd_table, cm_table,
                        alpha_min, dalpha, n_alpha)
    x  = (alpha - alpha_min) / dalpha
    x0 = clamp(x, zero(x), oftype(x, n_alpha - 1))
    i0 = unsafe_trunc(Int32, x0)
    # Two-sided index clamp. Besides the normal upper edge, this catches a
    # non-finite `alpha` (from a diverging solve): `clamp(NaN,…)` is NaN and
    # `unsafe_trunc(Int32, NaN)` is typemin(Int32), which would index far out of
    # bounds — a caught BoundsError on CPU but an illegal memory access under
    # the kernel's @inbounds on GPU. `t` stays NaN so the NaN still propagates
    # into the output (divergence remains visible) without any OOB read.
    i0 = clamp(i0, Int32(0), Int32(n_alpha - 2))
    t  = x0 - oftype(x0, i0)
    r0 = i0 + Int32(1)
    r1 = r0 + Int32(1)
    omt = one(t) - t
    cl = omt * cl_table[r0, section] + t * cl_table[r1, section]
    cd = omt * cd_table[r0, section] + t * cd_table[r1, section]
    cm = omt * cm_table[r0, section] + t * cm_table[r1, section]
    return cl, cd, cm
end

# blend_cosine (utils.jl:168), branch form. Bounds assumed in order (lb <= ub).
@inline function ds_blend_cosine(x, lb, ub)
    TF = typeof(x)
    if x >= ub
        return one(TF)
    elseif x <= lb
        return zero(TF)
    else
        cterm = cos((x - lb) * TF(pi) / (ub - lb))
        return (one(TF) - cterm) / TF(2)
    end
end

# ADGSP separation point (airfoils.jl:692). Reconstructs Cn from cl/cd tables.
@inline function ds_sep_adgsp(cl, cd, alpha, alpha0, Cd0, dcndalpha_circ)
    TF = typeof(alpha)
    delalpha = alpha - alpha0
    Cn = cl * cos(alpha) + (cd - Cd0) * sin(alpha)
    tr = if dcndalpha_circ == zero(TF) || alpha == alpha0 || Cn == zero(TF)
        zero(TF)
    else
        trr = Cn / (dcndalpha_circ * delalpha)
        trr < zero(TF) ? zero(TF) : trr
    end
    fst = ((TF(3) * sqrt(tr) - one(TF)) / TF(2))^2
    return fst > one(TF) ? one(TF) : fst
end

# ADGSP chordwise separation point (airfoils.jl:751).
@inline function ds_sepc_adgsp(cl, cd, alpha, alpha0, Cd0, eta, dcndalpha_circ)
    TF = typeof(alpha)
    Cc = cl * sin(alpha) - (cd - Cd0) * cos(alpha)
    D  = eta * dcndalpha_circ * (alpha - alpha0) * alpha
    fc = (Cc / D + TF(0.2))^2
    return min(fc, TF(DS_FCLIMIT))
end

# ---------------------------------------------------------------------------
# Ported ADG state update: sold (32) + (U, aoa, dt) -> snew (32).
# Direct port of update_states_ADG!. `cst` is the (DS_NCONST, n_sections)
# constant matrix; `j` is the section column.
# ---------------------------------------------------------------------------
@inline function ds_step_adg(sold, U, aoa, dt, j, cst,
                             cl_table, cd_table, cm_table, alpha_min, dalpha, n_alpha)
    TF = typeof(U)

    # --- constants ---
    dcndalpha = cst[C_DCNDALPHA, j]
    alpha0    = cst[C_ALPHA0, j]
    c         = cst[C_CHORD, j]
    a         = cst[C_ASOUND, j]
    A1 = cst[C_A1, j]; A2 = cst[C_A2, j]; A5 = cst[C_A5, j]
    b1 = cst[C_B1, j]; b2 = cst[C_B2, j]; b5 = cst[C_B5, j]
    Tp = cst[C_TP, j]; Tf0 = cst[C_TF0, j]; Tv0 = cst[C_TV0, j]
    Tvl = cst[C_TVL, j]; Tsh = cst[C_TSH, j]
    zeta = cst[C_ZETA, j]; Cd0 = cst[C_CD0, j]; Cm0 = cst[C_CM0, j]; Cn1 = cst[C_CN1, j]

    # --- old states (sold[1] unused; reset to aoa below) ---
    alpha_m = sold[2];  q_m = sold[3];  qf_m = sold[4]
    Ka_m = sold[5];     Kq_m = sold[6]; Kpa_m = sold[7]; Kpq_m = sold[8]
    X1_m = sold[9];     X2_m = sold[10]; X3_m = sold[11]; X4_m = sold[12]
    Npot_m = sold[13];  Kppq_m = sold[14]; Kpppq_m = sold[15]
    Dp_m = sold[16];    fp_m = sold[17]; fpc_m = sold[18]; fpm_m = sold[19]
    Df_m = sold[20];    Dfc_m = sold[21]; Dfm_m = sold[22]
    fpp_m = sold[23]   # sold[24] (fppc_m) and sold[25] (fppm_m) are not fed back by ADG
    Cv_m = sold[26];    Nv_m = sold[27]; tauv = sold[28]
    LESF_m = sold[29];  TESF_m = sold[30]; VRTX_m = sold[31]; firstpass_m = sold[32]

    TI     = c / a
    deltas = TF(2) * U * dt / c

    s1 = aoa   # states[1]

    # low-pass filtering
    lp_cutoff = max(one(TF), U) * zeta / (TF(pi) * c)
    Clp = exp(-TF(2) * TF(pi) * dt * lp_cutoff)
    plC = one(TF) - Clp

    alpha = Clp * alpha_m + plC * aoa
    s2 = alpha

    Ka = (alpha - alpha_m) / dt
    q  = Ka * c / U
    s3 = q
    if firstpass_m == one(TF)
        q_m = q
    end
    qf = Clp * qf_m + plC * q
    s4 = qf

    Ka = qf * U / c
    s5 = Ka

    Kq = (q - q_m) / dt
    if firstpass_m == one(TF)
        Kq_m = zero(TF)
    end
    Kq = Clp * Kq_m + plC * Kq
    s6 = Kq

    # sigma1 (Tf) / sigma3 (Tv) modifications
    Delta_alpha0 = alpha - alpha0
    sigma1 = one(TF); sigma3 = one(TF)
    sigma1c = one(TF); sigma1m = one(TF)

    if TESF_m == one(TF)
        if Ka_m * Delta_alpha0 < zero(TF)
            sigma1 = TF(2)
        else
            if LESF_m == zero(TF)
                sigma1 = one(TF)
            else
                sigma1 = fpp_m <= TF(0.7) ? TF(2) : TF(1.75)
            end
        end
    else
        if LESF_m == zero(TF)
            sigma1 = TF(0.5)
        elseif VRTX_m == one(TF) && zero(TF) <= tauv <= Tvl
            sigma1 = TF(0.25)
        elseif Ka_m * Delta_alpha0 > zero(TF)
            sigma1 = TF(0.75)
        end
    end

    if Tvl <= tauv <= TF(2) * Tvl
        sigma3 = TF(3)
        if TESF_m == zero(TF)
            sigma3 = TF(4)
            if VRTX_m == one(TF) && zero(TF) <= tauv <= Tvl
                sigma3 = Ka_m * Delta_alpha0 < zero(TF) ? TF(2) : one(TF)
            end
        end
    else
        if Ka_m * Delta_alpha0 < zero(TF)
            sigma3 = TF(4)
        end
    end
    if TESF_m == zero(TF) && Kq_m * Delta_alpha0 < zero(TF)
        sigma3 = one(TF)
    end

    M = U / a
    beta = sqrt(max(one(TF) - M * M, zero(TF)))   # guarded (M>1 path discarded at end)

    bot = (one(TF) - M) + dcndalpha * (M * M) * beta * (A1 * b1 + A2 * b2) / TF(2)
    k_alpha = one(TF) / bot
    bot = (one(TF) - M) + dcndalpha * M * M * beta * (A1 * b1 + A2 * b2)
    k_q = one(TF) / bot

    Talpha = TF(3) * k_alpha * TI / TF(4)
    Tq     = TF(3) * k_q * TI / TF(4)
    Tf  = Tf0 / sigma1
    Tv  = Tv0 / sigma3
    Tfc = Tf0 / sigma1c
    Tfm = Tf0 / sigma1m

    # non-circulatory normal force
    Kpa = Kpa_m * exp(-dt / Talpha) + (Ka - Ka_m) * exp(-dt / (TF(2) * Talpha))
    s7 = Kpa
    Kpq = Kpq_m * exp(-dt / Tq) + (Kq - Kq_m) * exp(-dt / (TF(2) * Tq))
    s8 = Kpq
    Nnoncirc_a  = TF(4) * Talpha * (Ka - Kpa) / M
    Nnoncirc_q  = -Tq * (Kq - Kpq) / M
    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q

    beta2 = beta^2
    delta_alpha = alpha - alpha_m
    deltaq = qf - qf_m

    X1 = X1_m * exp(-b1 * beta2 * deltas) + A1 * exp(-b1 * beta2 * deltas / TF(2)) * delta_alpha
    s9 = X1
    X2 = X2_m * exp(-b2 * beta2 * deltas) + A2 * exp(-b2 * beta2 * deltas / TF(2)) * delta_alpha
    s10 = X2
    X3 = X3_m * exp(-b1 * beta2 * deltas) + A1 * exp(-b1 * beta2 * deltas / TF(2)) * deltaq
    s11 = X3
    X4 = X4_m * exp(-b2 * beta2 * deltas) + A2 * exp(-b2 * beta2 * deltas / TF(2)) * deltaq
    s12 = X4

    alphae = (alpha - alpha0) - X1 - X2
    dcndalpha_circ = dcndalpha / beta
    Ncirc_aq = dcndalpha_circ * alphae
    Npot = Ncirc_aq + Nnoncirc_aq
    s13 = Npot

    bot = TF(15) * (one(TF) - M) + TF(3) * dcndalpha * A5 * b5 * beta * M * M / TF(2)
    kmq = TF(7) / bot
    Kppq = Kppq_m * exp(-dt / (kmq * kmq * TI)) + (Kq - Kq_m) * exp(-dt / (TF(2) * kmq * kmq * TI))
    s14 = Kppq
    Kpppq = Kpppq_m * exp(-b5 * beta2 * deltas) + A5 * deltaq * exp(-b5 * beta2 * deltas / TF(2))
    s15 = Kpppq

    # boundary layer response
    if firstpass_m == one(TF)
        Npot_m = Npot
    end
    Dp = Dp_m * exp(-deltas / Tp) + (Npot - Npot_m) * exp(-deltas / (Tp * TF(2)))
    s16 = Dp
    Cpn = Npot - Dp
    alphaf = Cpn / dcndalpha_circ + alpha0

    # separation points at alphaf (table lookups)
    clf, cdf, cmf = ds_af3(alphaf, j, cl_table, cd_table, cm_table, alpha_min, dalpha, n_alpha)
    eta_c = cst[C_ETA, j]
    fp  = ds_sep_adgsp(clf, cdf, alphaf, alpha0, Cd0, dcndalpha_circ)
    s17 = fp
    fpc = ds_sepc_adgsp(clf, cdf, alphaf, alpha0, Cd0, eta_c, dcndalpha_circ)
    s18 = fpc

    Cntemp = clf * cos(alphaf) + (cdf - Cd0) * sin(alphaf)
    fpm = abs(Cntemp) < TF(0.01) ? zero(TF) : (cmf - Cm0) / Cntemp
    s19 = fpm

    Df  = firstpass_m == one(TF) ? zero(TF) :
          Df_m  * exp(-deltas / Tf)  + (fp  - fp_m)  * exp(-deltas / (TF(2) * Tf))
    s20 = Df
    Dfc = firstpass_m == one(TF) ? zero(TF) :
          Dfc_m * exp(-deltas / Tfc) + (fpc - fpc_m) * exp(-deltas / (TF(2) * Tfc))
    s21 = Dfc
    Dfm = firstpass_m == one(TF) ? zero(TF) :
          Dfm_m * exp(-deltas / Tfm) + (fpm - fpm_m) * exp(-deltas / (TF(2) * Tfm))
    s22 = Dfm

    fpp  = fp  - Df;   s23 = fpp
    fppc = fpc - Dfc;  s24 = fppc
    fppm = fpm - Dfm;  s25 = fppm

    fterm = (one(TF) + TF(2) * sqrt(max(fpp, zero(TF)))) / TF(3)
    Cv = dcndalpha_circ * alphae * (one(TF) - fterm)^2
    s26 = Cv

    Nv = if firstpass_m == one(TF)
        zero(TF)
    elseif (tauv > Tvl) && (Ka * Delta_alpha0 > zero(TF))
        Nv_m * exp(-TF(2) * deltas / Tv)
    else
        Nv_m * exp(-deltas / Tv) + (Cv - Cv_m) * exp(-deltas / (TF(2) * Tv))
    end
    if Nv < zero(TF)
        Nv = zero(TF)
    end
    s27 = Nv

    # other states
    LESF = Cpn > Cn1 ? one(TF) : zero(TF)
    TESF = fpp < fpp_m ? one(TF) : zero(TF)
    VRTX = (zero(TF) < tauv <= TF(2) * Tvl) ? one(TF) : zero(TF)

    if (tauv >= one(TF) + Tsh / Tvl) && (LESF == one(TF))
        tauv = zero(TF)
    end
    if LESF == one(TF)
        tauv += dt * TF(2) * U / c
    end
    s28 = tauv
    s29 = LESF
    s30 = TESF
    s31 = VRTX
    s32 = zero(TF)   # firstpass

    snew = (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16,
            s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31, s32)

    # M>1 path: leave states unchanged (matches CPU `states .= oldstates; return`).
    return M > one(TF) ? sold : snew
end

# ---------------------------------------------------------------------------
# Ported load calculation: states (32) + (U) -> (Cl, Cd, Cm).
# Direct port of BLADG_coefficients.
# ---------------------------------------------------------------------------
@inline function ds_loads_adg(s, U, j, cst,
                              cl_table, cd_table, cm_table, alpha_min, dalpha, n_alpha)
    TF = typeof(U)

    dcndalpha = cst[C_DCNDALPHA, j]
    alpha0    = cst[C_ALPHA0, j]
    c         = cst[C_CHORD, j]
    xcp       = cst[C_XCP, j]
    a         = cst[C_ASOUND, j]
    alphacut2 = cst[C_ALPHACUT2, j]
    cutrad    = cst[C_CUTRAD, j]
    A1 = cst[C_A1, j]; A2 = cst[C_A2, j]; A5 = cst[C_A5, j]
    b1 = cst[C_B1, j]; b2 = cst[C_B2, j]; b5 = cst[C_B5, j]
    Tvl = cst[C_TVL, j]; eta = cst[C_ETA, j]
    Cd0 = cst[C_CD0, j]; Cm0 = cst[C_CM0, j]

    aoa = s[1]; alpha = s[2]
    qf = s[4]; Ka = s[5]; Kq = s[6]; Kpa = s[7]; Kpq = s[8]
    X1 = s[9]; X2 = s[10]; X3 = s[11]; X4 = s[12]
    Kppq = s[14]; Kpppq = s[15]
    fpp = s[23]; fppc = s[24]; fppm = s[25]
    Cvn = s[27]; tauv = s[28]; firstpass = s[32]

    M = U / a
    beta = sqrt(max(one(TF) - M * M, zero(TF)))
    TI = c / a

    bot = (one(TF) - M) + dcndalpha * M * M * beta * (A1 * b1 + A2 * b2) / TF(2)
    k_alpha = one(TF) / bot
    Talpha = TF(3) * k_alpha * TI / TF(4)
    Nnoncirc_a = TF(4) * Talpha * (Ka - Kpa) / M

    bot = (one(TF) - M) + dcndalpha * M * M * beta * (A1 * b1 + A2 * b2)
    k_q = one(TF) / bot
    Tq = TF(3) * k_q * TI / TF(4)
    Nnoncirc_q = -Tq * (Kq - Kpq) / M
    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q

    alphae = (alpha - alpha0) - X1 - X2
    dcndalpha_circ = dcndalpha / beta
    Ncirc_q = dcndalpha_circ * qf / TF(2) - X3 - X4
    Ncirc_aq = dcndalpha_circ * alphae

    fterm = (one(TF) + TF(2) * sqrt(max(fpp, zero(TF)))) / TF(3)
    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq * (fterm^2)

    Cn = tauv > zero(TF) ? Cfsn + Cvn : Cfsn

    Cpot = dcndalpha_circ * alphae * aoa
    Cfsc = Cpot * eta * (sqrt(max(fppc, zero(TF))) - TF(0.2))
    Cc = Cfsc

    Cl = Cn * cos(aoa) + Cc * sin(aoa)
    Cd = Cn * sin(aoa) - Cc * cos(aoa) + Cd0

    cla, cda, cma = ds_af3(aoa, j, cl_table, cd_table, cm_table, alpha_min, dalpha, n_alpha)
    if firstpass == one(TF)
        Cl = cla
        Cd = cda
    end

    Mcirc_q = -dcndalpha * (qf - Kpppq) * c / (TF(16) * beta * U)
    Mnoncirc_alpha = -Nnoncirc_a / TF(4)
    bot = TF(15) * (one(TF) - M) + TF(3) * dcndalpha * A5 * b5 * beta * M * M / TF(2)
    kmq = TF(7) / bot
    Mnoncirc_q = -TF(7) * TI * kmq * kmq * (Kq - Kppq) / (TF(12) * M)

    Cm_common = Mcirc_q + Mnoncirc_alpha + Mnoncirc_q
    Mfs = Cm0 + Cfsn * fppm + Cm_common
    Mv = xcp * (one(TF) - cos(TF(pi) * tauv / Tvl)) * Cvn
    Cm = tauv <= zero(TF) ? Mfs : Mfs + Mv

    W1 = one(TF) - ds_blend_cosine(aoa, alphacut2 - cutrad, alphacut2)
    W2 = ds_blend_cosine(U, zero(TF), one(TF))
    W = W1 * W2
    if W < one(TF)
        Cl = Cl * W + (one(TF) - W) * cla
        Cd = Cd * W + (one(TF) - W) * cda
        Cm = Cm * W + (one(TF) - W) * cma
    end

    return Cl, Cd, Cm
end

# ---------------------------------------------------------------------------
# Ported initial condition: (U, aoa) -> snew (32), (Cl, Cd, Cm).
# Direct port of initialize_ADG.
# ---------------------------------------------------------------------------
@inline function ds_init_adg(U, aoa, j, cst,
                             cl_table, cd_table, cm_table, alpha_min, dalpha, n_alpha)
    TF = typeof(U)
    Cd0 = cst[C_CD0, j]
    z = zero(TF)

    alpha = aoa
    cl0, cd0, cm0 = ds_af3(aoa, j, cl_table, cd_table, cm_table, alpha_min, dalpha, n_alpha)
    Cpotn = cl0 + (cd0 - Cd0) * sin(alpha)   # state 13 (Npot); note: verbatim from initialize_ADG

    fclim = TF(DS_FCLIMIT)
    snew = (aoa,          # 1  aoa
            alpha,        # 2  alpha
            z,            # 3  q
            z,            # 4  qf
            z,            # 5  Ka
            z,            # 6  Kq
            z,            # 7  Kpa
            z,            # 8  Kpq
            z,            # 9  X1
            z,            # 10 X2
            z,            # 11 X3
            z,            # 12 X4
            Cpotn,        # 13 Npot
            z,            # 14 Kppq
            z,            # 15 Kpppq
            z,            # 16 Dp
            one(TF),      # 17 fp
            fclim,        # 18 fpc
            z,            # 19 fpm
            z,            # 20 Df
            z,            # 21 Dfc
            z,            # 22 Dfm
            one(TF),      # 23 fpp
            fclim,        # 24 fppc
            z,            # 25 fppm
            z,            # 26 Cv
            z,            # 27 Nv
            z,            # 28 tauv
            z,            # 29 LESF
            z,            # 30 TESF
            z,            # 31 VRTX
            one(TF))      # 32 firstpass
    return snew, cl0, cd0, cm0
end

# ---------------------------------------------------------------------------
# Kernels and host driver
# ---------------------------------------------------------------------------

# Fill step 1 (quasi-steady initial condition) for every (section, sim).
@kernel function dsmodel_init_kernel!(xds, Cl, Cd, Cm, U, aoa,
        consts, cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
    j, s = @index(Global, NTuple)
    @inbounds begin
        Uv = U[j, s, 1]
        av = aoa[j, s, 1]
        snew, cl0, cd0, cm0 = ds_init_adg(Uv, av, Int32(j), consts,
                                          cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
        Base.Cartesian.@nexprs 32 k -> (xds[k, j, s, 1] = snew[k])
        Cl[j, s, 1] = cl0; Cd[j, s, 1] = cd0; Cm[j, s, 1] = cm0
    end
end

# Advance every (section, sim) from step i-1 to step i.
@kernel function dsmodel_step_kernel!(xds, Cl, Cd, Cm, U, aoa, dt,
        consts, cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha, i)
    j, s = @index(Global, NTuple)
    @inbounds begin
        im1  = i - 1
        sold = ntuple(k -> xds[k, j, s, im1], Val(DS_NSTATES))
        Uv = U[j, s, i]; av = aoa[j, s, i]
        snew = ds_step_adg(sold, Uv, av, dt, Int32(j), consts,
                           cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
        Base.Cartesian.@nexprs 32 k -> (xds[k, j, s, i] = snew[k])
        cl, cd, cm = ds_loads_adg(snew, Uv, Int32(j), consts,
                                  cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
        Cl[j, s, i] = cl; Cd[j, s, i] = cd; Cm[j, s, i] = cm
    end
end

"""
    march_ds_gpu!(hist::DSHistory, dsaf::DSAirfoilGPU, U, aoa, dt)

March the batched Beddoes-Leishman v3 model over every `(section, sim)` pair
for all `nt` steps, writing the full state history and `(Cl, Cd, Cm)` into
`hist`. `U`, `aoa` are `(n_sections, n_sims, nt)` device arrays on the same
backend as `hist`/`dsaf`; `dt` is the (fixed) step size.

Step 1 is the quasi-steady initial condition; steps `2:nt` are advanced by a
serial launch loop (step `i` reads step `i-1`, so launches queue in order on
the backend's stream and a single final synchronize suffices).
"""
function march_ds_gpu!(hist::DSHistory, dsaf::DSAirfoilGPU,
                       U::AbstractArray, aoa::AbstractArray, dt)
    n_sections, n_sims, nt = size(U)
    size(aoa) == (n_sections, n_sims, nt) ||
        error("march_ds_gpu!: size(aoa) $(size(aoa)) != size(U) $((n_sections, n_sims, nt))")
    size(hist.xds) == (DS_NSTATES, n_sections, n_sims, nt) ||
        error("march_ds_gpu!: hist.xds is $(size(hist.xds)), expected $((DS_NSTATES, n_sections, n_sims, nt)).")

    TF = eltype(hist.xds)
    dtv = TF(dt)
    backend = KernelAbstractions.get_backend(hist.xds)
    initk = dsmodel_init_kernel!(backend)
    stepk = dsmodel_step_kernel!(backend)

    initk(hist.xds, hist.Cl, hist.Cd, hist.Cm, U, aoa,
          dsaf.consts, dsaf.cl_table, dsaf.cd_table, dsaf.cm_table,
          dsaf.alpha_min, dsaf.dalpha, dsaf.n_alpha; ndrange=(n_sections, n_sims))
    for i in 2:nt
        stepk(hist.xds, hist.Cl, hist.Cd, hist.Cm, U, aoa, dtv,
              dsaf.consts, dsaf.cl_table, dsaf.cd_table, dsaf.cm_table,
              dsaf.alpha_min, dsaf.dalpha, dsaf.n_alpha, i; ndrange=(n_sections, n_sims))
    end
    KernelAbstractions.synchronize(backend)
    return hist
end

# ---------------------------------------------------------------------------
# Per-step (ping-pong) kernels for the coupled aeroelastic solver.
#
# The march above precomputes the whole (U, aoa) history and stores the full
# 4D state trace. In a coupled sim the inputs at step i depend on the
# structural feedback from step i-1, so U/aoa aren't known ahead of time and
# the march must advance one step at a time. These kernels take 3D state
# buffers `(DS_NSTATES, n_sections, n_sims)` — a current/previous ping-pong
# pair — and 2D inputs/outputs `(n_sections, n_sims)` for the current step.
#
# They call the identical `ds_init_adg`/`ds_step_adg`/`ds_loads_adg` functions
# as the march kernels, so a per-step drive over a prescribed history is
# bit-identical to `march_ds_gpu!` (verified in the Phase-1 unit test).
# ---------------------------------------------------------------------------

# Step-1 quasi-steady IC into a 3D state buffer.
@kernel function dsmodel_init_step_kernel!(xds_cur, Cl, Cd, Cm, U, aoa,
        consts, cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
    j, s = @index(Global, NTuple)
    @inbounds begin
        Uv = U[j, s]
        av = aoa[j, s]
        snew, cl0, cd0, cm0 = ds_init_adg(Uv, av, Int32(j), consts,
                                          cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
        Base.Cartesian.@nexprs 32 k -> (xds_cur[k, j, s] = snew[k])
        Cl[j, s] = cl0; Cd[j, s] = cd0; Cm[j, s] = cm0
    end
end

# Advance every (section, sim) one step: read `xds_prev`, write `xds_cur`.
@kernel function dsmodel_step_step_kernel!(xds_cur, xds_prev, Cl, Cd, Cm, U, aoa, dt,
        consts, cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
    j, s = @index(Global, NTuple)
    @inbounds begin
        sold = ntuple(k -> xds_prev[k, j, s], Val(DS_NSTATES))
        Uv = U[j, s]; av = aoa[j, s]
        snew = ds_step_adg(sold, Uv, av, dt, Int32(j), consts,
                           cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
        Base.Cartesian.@nexprs 32 k -> (xds_cur[k, j, s] = snew[k])
        cl, cd, cm = ds_loads_adg(snew, Uv, Int32(j), consts,
                                  cl_t, cd_t, cm_t, alpha_min, dalpha, n_alpha)
        Cl[j, s] = cl; Cd[j, s] = cd; Cm[j, s] = cm
    end
end

"""
    ds_init_step_gpu!(xds_cur, Cl, Cd, Cm, dsaf, U, aoa)

Fill the quasi-steady dynamic-stall initial condition for every `(section, sim)`
into the 3D state buffer `xds_cur` `(DS_NSTATES, n_sections, n_sims)` and write
the step-1 load coefficients into the 2D buffers `Cl/Cd/Cm` `(n_sections, n_sims)`.
`U`, `aoa` are the current-step `(n_sections, n_sims)` inputs. All arrays must
share the backend of `xds_cur`. This is the per-step (ping-pong) companion to
[`march_ds_gpu!`](@ref) used by the coupled aeroelastic solver.
"""
function ds_init_step_gpu!(xds_cur::AbstractArray, Cl::AbstractMatrix, Cd::AbstractMatrix,
                           Cm::AbstractMatrix, dsaf::DSAirfoilGPU,
                           U::AbstractMatrix, aoa::AbstractMatrix)
    n_sections, n_sims = size(U)
    size(xds_cur) == (DS_NSTATES, n_sections, n_sims) ||
        error("ds_init_step_gpu!: xds_cur is $(size(xds_cur)), expected $((DS_NSTATES, n_sections, n_sims)).")
    backend = KernelAbstractions.get_backend(xds_cur)
    initk = dsmodel_init_step_kernel!(backend)
    initk(xds_cur, Cl, Cd, Cm, U, aoa,
          dsaf.consts, dsaf.cl_table, dsaf.cd_table, dsaf.cm_table,
          dsaf.alpha_min, dsaf.dalpha, dsaf.n_alpha; ndrange=(n_sections, n_sims))
    KernelAbstractions.synchronize(backend)
    return xds_cur
end

"""
    ds_step_gpu!(xds_cur, xds_prev, Cl, Cd, Cm, dsaf, U, aoa, dt)

Advance the batched Beddoes-Leishman v3 model one step for every `(section, sim)`:
read the previous state `xds_prev`, write the new state into `xds_cur`, and
write the new load coefficients into `Cl/Cd/Cm`. All state buffers are 3D
`(DS_NSTATES, n_sections, n_sims)`; inputs/outputs are 2D `(n_sections, n_sims)`.

Caller owns the ping-pong: pass distinct `xds_cur`/`xds_prev` buffers and swap
them between steps. This advances identically to one iteration of
[`march_ds_gpu!`](@ref)'s inner loop.
"""
function ds_step_gpu!(xds_cur::AbstractArray, xds_prev::AbstractArray,
                      Cl::AbstractMatrix, Cd::AbstractMatrix, Cm::AbstractMatrix,
                      dsaf::DSAirfoilGPU, U::AbstractMatrix, aoa::AbstractMatrix, dt)
    n_sections, n_sims = size(U)
    size(xds_cur) == (DS_NSTATES, n_sections, n_sims) ||
        error("ds_step_gpu!: xds_cur is $(size(xds_cur)), expected $((DS_NSTATES, n_sections, n_sims)).")
    size(xds_prev) == size(xds_cur) ||
        error("ds_step_gpu!: xds_prev $(size(xds_prev)) != xds_cur $(size(xds_cur)).")
    TF = eltype(xds_cur)
    backend = KernelAbstractions.get_backend(xds_cur)
    stepk = dsmodel_step_step_kernel!(backend)
    stepk(xds_cur, xds_prev, Cl, Cd, Cm, U, aoa, TF(dt),
          dsaf.consts, dsaf.cl_table, dsaf.cd_table, dsaf.cm_table,
          dsaf.alpha_min, dsaf.dalpha, dsaf.n_alpha; ndrange=(n_sections, n_sims))
    KernelAbstractions.synchronize(backend)
    return xds_cur
end
