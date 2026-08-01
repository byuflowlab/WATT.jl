#=
Correctness check for GPU-BEMT.

Runs `WATT.solve_BEMT_gpu!` against the CPU `WATT.solve_BEMT` (which wraps
CCBlade) across a sweep of operating conditions. Reports max abs / max rel
differences on the quantities that matter to the coupled sim (phi, Np, Tp,
W, alpha).

Two independent toggles at the top:
- BACKEND  — which KernelAbstractions backend to run the GPU kernel on
- DEVFLOAT — device eltype (Float32 or Float64). Decoupled from BACKEND so
             you can pin down whether a discrepancy comes from precision
             (Float32 vs Float64) or from the GPU itself (CPU-KA vs CUDA
             vs Metal, holding DEVFLOAT constant).

BACKEND options:
- :cpu    — KernelAbstractions CPU backend (Array{DEVFLOAT}). Local dry-run.
- :cuda   — CUDA.jl / CuArray{DEVFLOAT}. Supercomputer NVIDIA.
- :amdgpu — AMDGPU.jl / ROCArray{DEVFLOAT}.
- :metal  — Metal.jl / MtlArray{DEVFLOAT}.

The CPU reference always runs Float64. When DEVFLOAT is Float32, part of
the residual difference is Float32 roundoff (~1e-5 relative floor is
typical); part is the linear-vs-Akima airfoil interp on the pre-sampled
table (controlled by n_alpha).

Run:
    julia --project examples/gpu_bemt_cpu_backend_check.jl
=#

# ---------------------------------------------------------------------------
# Top-level knobs
# ---------------------------------------------------------------------------
const BACKEND  = :cuda      # :cpu | :cuda | :amdgpu | :metal
const DEVFLOAT = Float32    # Float32 or Float64
const DIAGNOSE = true       # true = also dump residual snapshots at known-bad (j,s)

using WATT, OpenFASTTools, DynamicStallModels
using StaticArrays
using Printf

const _of = OpenFASTTools
const _DS = DynamicStallModels
const PrandtlTipHub = WATT.CCBlade.PrandtlTipHub

# ---------------------------------------------------------------------------
# Backend glue (mirror of gpu_bemt_benchmark.jl). Loads the vendor package,
# picks a device array type of the requested eltype, and adds three tiny
# `WATT.*` methods so `BladeGPU` / `GPUBEMTOutputs` can pack their host arrays
# onto that device. Also defines `device_sync()` (called after each kernel to
# ensure timing/comparison sees the finished result) and `to_host(...)`
# (moves a device array back to host for the reference comparison).
# ---------------------------------------------------------------------------
if BACKEND === :cuda
    using CUDA
    const DevFloat = DEVFLOAT
    const DevArrayType = CuArray{DevFloat}
    WATT.to_backend_vector(::Type{CuArray{T}}, v::AbstractVector) where {T} = CuArray{T}(v)
    WATT.to_backend_matrix(::Type{CuArray{T}}, m::AbstractMatrix) where {T} = CuArray{T}(m)
    WATT.similar_type(::Type{CuArray{T}}, ::Type{U}) where {T, U} = CuArray{U}
    device_sync() = CUDA.synchronize()
    to_host(a::AbstractArray) = Array(a)
elseif BACKEND === :amdgpu
    using AMDGPU
    const DevFloat = DEVFLOAT
    const DevArrayType = ROCArray{DevFloat}
    WATT.to_backend_vector(::Type{ROCArray{T}}, v::AbstractVector) where {T} = ROCArray{T}(v)
    WATT.to_backend_matrix(::Type{ROCArray{T}}, m::AbstractMatrix) where {T} = ROCArray{T}(m)
    WATT.similar_type(::Type{ROCArray{T}}, ::Type{U}) where {T, U} = ROCArray{U}
    device_sync() = AMDGPU.synchronize()
    to_host(a::AbstractArray) = Array(a)
elseif BACKEND === :metal
    using Metal
    const DevFloat = DEVFLOAT
    const DevArrayType = MtlArray{DevFloat}
    WATT.to_backend_vector(::Type{MtlArray{T}}, v::AbstractVector) where {T} = MtlArray{T}(v)
    WATT.to_backend_matrix(::Type{MtlArray{T}}, m::AbstractMatrix) where {T} = MtlArray{T}(m)
    WATT.similar_type(::Type{MtlArray{T}}, ::Type{U}) where {T, U} = MtlArray{U}
    device_sync() = Metal.synchronize()
    to_host(a::AbstractArray) = Array(a)
else
    const DevFloat = DEVFLOAT
    const DevArrayType = Array{DevFloat}
    device_sync() = nothing
    to_host(a::AbstractArray) = a
end

@info "Backend selected" BACKEND DevArrayType DevFloat

# ---------------------------------------------------------------------------
# Build the NREL 5MW rotor from the OpenFAST inputs in data/openfast.
# This is the same setup used by gpu_bemt_bracket_check.jl so the two scripts
# see the same geometry (25 aero nodes, cylinder root through DU/NACA
# airfoils out to the tip).
# ---------------------------------------------------------------------------
ofpath = joinpath(@__DIR__, "..", "data", "openfast")
adblade = _of.read_adblade("sn5_adblade.dat", ofpath)
edfile  = _of.read_edfile("sn5_EDfile.dat", ofpath)

aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat", "DU25_A17.dat", "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [_of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]
af_idx = Int.(adblade["BlAFID"])
afs    = aftypes[af_idx]

chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"] .* (pi/180)     # OpenFAST twist is in degrees
rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub               # aero station absolute radii
rtip = rvec[end]
n    = length(rvec)

airfoils = Vector{_DS.Airfoil}(undef, n)
xcp = Vector{Float64}(undef, n)
for i = 1:n
    airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
end

blade = WATT.Blade(rvec, chordvec, twistvec, xcp, airfoils; rhub, rtip)

# Exercise the two supported tip-correction paths. In the "no tip loss" case
# the residual has no Prandtl singularity, so we expect every (j,s) to solve.
# In the tip+hub case, the aero nodes AT r=Rhub and r=Rtip make F → 0, which
# makes the residual undefined there — those sections should hit
# bracket-failed (about n_sims × 2 = 250 misses expected).
rotor_notip = WATT.Rotor(3, 80.0, true)                    # tip=nothing
rotor_tip   = WATT.Rotor(3, 80.0, true; tip=PrandtlTipHub())

# Only rho/mu/a are used from the environment; wind/omega args here are
# placeholders because get_aero_velocities isn't called in this script.
env = WATT.environment(1.225, 1.464e-5, 343.0, 10.0, 1.0, 0.0)

# --- Build the GPU-BEMT side ------------------------------------------------
rotorgpu_notip = WATT.RotorGPU(rotor_notip)
rotorgpu_tip   = WATT.RotorGPU(rotor_tip)

# ---------------------------------------------------------------------------
# Build a batch of (Vx, Vy, pitch) inputs.
#
# One "sim" per (wind, TSR, pitch) tuple. Cartesian product of the three
# vectors below gives 5×5×5 = 125 sims. Every sim shares the same 25 aero
# sections, so the (n_sections, n_sims) = (25, 125) matrices are what the
# kernel operates on.
#
# Vx = freestream velocity at the section (constant along the blade)
# Vy = Omega * r_j        (rotational velocity in the plane of rotation)
# pitch is per-sim (constant along the blade)
# ---------------------------------------------------------------------------
wind_speeds = [5.0, 8.0, 11.4, 15.0, 20.0]
tsrs        = [4.0, 6.0, 7.5, 9.0, 11.0]
pitches_deg = [-1.0, 0.0, 2.0, 5.0, 10.0]

Us      = Float64[]
Omegas  = Float64[]
pitches = Float64[]
for U in wind_speeds, tsr in tsrs, p in pitches_deg
    push!(Us, U)
    push!(Omegas, tsr * U / rtip)
    push!(pitches, p * pi/180)
end
n_sims = length(Us)

# CPU reference inputs (Float64, always host).
Vx_h = zeros(Float64, n, n_sims)
Vy_h = zeros(Float64, n, n_sims)
for s in 1:n_sims, j in 1:n
    Vx_h[j, s] = Us[s]
    Vy_h[j, s] = Omegas[s] * blade.r[j]
end
pitch_h = collect(pitches)

# Device inputs at DevFloat precision on the selected backend.
Vx    = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vx_h))
Vy    = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vy_h))
pitch = WATT.to_backend_vector(DevArrayType, DevFloat.(pitch_h))

# ---------------------------------------------------------------------------
# `compare(...)` — run the GPU BEMT once and compare against the CPU BEMT.
#
# Pipeline:
#   1. Allocate a fresh `GPUBEMTOutputs` on the selected backend.
#   2. Launch the kernel (`solve_BEMT_gpu!`) and synchronize.
#   3. Copy GPU outputs to host in Float64 so the comparison is not skewed
#      by DevFloat precision on the difference itself.
#   4. Run the reference CPU BEMT for every (j, s), collecting phi, alpha, W,
#      Np, Tp.
#   5. Report:
#        - bracket-ok count      : how many (j, s) the kernel's sub-bracket
#                                  scan found a residual sign change in q1.
#                                  Failed sections write zeros — matches CPU
#                                  CCBlade.Outputs() behavior for singular
#                                  sections — but we mask them out here so
#                                  the report shows differences only where
#                                  both paths committed to a value.
#        - max abs / max rel     : location + values of the worst-case
#                                  section-sim pair for each quantity. rel
#                                  denominator has a 1e-8 floor to avoid
#                                  divide-by-tiny at zero-load sections;
#                                  see the report for the CPU value there
#                                  to spot false-alarm "big" relative errors.
# ---------------------------------------------------------------------------
function compare(label, rotor_cpu, rotor_gpu, blade_gpu, n_iters)
    outputs = WATT.GPUBEMTOutputs(n, n_sims; ArrayType=DevArrayType)
    WATT.solve_BEMT_gpu!(outputs, rotor_gpu, blade_gpu, env, Vx, Vy, pitch; n_iters=n_iters)
    device_sync()

    # Pull to host in Float64 to compare on equal footing
    phi_g   = Float64.(to_host(outputs.phi))
    alpha_g = Float64.(to_host(outputs.alpha))
    W_g     = Float64.(to_host(outputs.W))
    Np_g    = Float64.(to_host(outputs.Np))
    Tp_g    = Float64.(to_host(outputs.Tp))
    ok      = to_host(outputs.success)

    # CPU reference (Float64)
    phi_cpu = zeros(n, n_sims); Np_cpu = similar(phi_cpu); Tp_cpu = similar(phi_cpu)
    W_cpu   = similar(phi_cpu); alpha_cpu = similar(phi_cpu)
    xv = zeros(11)   # scratch vector layout expected by CCBlade.residual_and_outputs
    for s in 1:n_sims, j in 1:n
        ccout = WATT.solve_BEMT(rotor_cpu, blade, env, j, Vx_h[j, s], Vy_h[j, s], pitch_h[s], xv)
        phi_cpu[j, s]   = ccout.phi
        Np_cpu[j, s]    = ccout.Np
        Tp_cpu[j, s]    = ccout.Tp
        W_cpu[j, s]     = ccout.W
        alpha_cpu[j, s] = ccout.alpha
    end

    println("\n=== $label  (n_iters=$n_iters, tip=$(rotor_cpu.tip)) ===")
    n_ok = count(ok); n_bad = length(ok) - n_ok
    println("bracket-ok: $n_ok / $(length(ok))   (bracket-failed: $n_bad)")

    # Zero out differences at sections the kernel gave up on — they would
    # otherwise dominate the max as "gpu=0 vs cpu=<solved value>" comparisons,
    # which is not what we want to report here.
    mask = ok
    for (name, gpu, cpu) in (
            ("phi",   phi_g,   phi_cpu),
            ("alpha", alpha_g, alpha_cpu),
            ("W",     W_g,     W_cpu),
            ("Np",    Np_g,    Np_cpu),
            ("Tp",    Tp_g,    Tp_cpu),
        )
        # max abs = worst absolute difference across all (j,s)
        d = abs.(gpu .- cpu) .* mask
        rmax, imax = findmax(d)
        # max rel = worst relative-to-|cpu| difference. The 1e-8 floor keeps
        # zero-load sections from producing spurious infinities; look at the
        # cpu column in the report to see whether a big rel number actually
        # matters (~10 N/m difference at a 6000 N/m section is 0.2% real error;
        # ~1 N/m difference at 4 N/m is 25% "rel" but physically meaningless).
        denom = max.(abs.(cpu), 1e-8)
        rrel_arr = (d ./ denom) .* mask
        rrel, irrel = findmax(rrel_arr)
        j_abs, s_abs = Tuple(imax)
        j_rel, s_rel = Tuple(irrel)
        @printf("  %-6s  max abs %.3e @ (j=%d,s=%d, cpu=%+.3e gpu=%+.3e)  max rel %.3e @ (j=%d,s=%d, cpu=%+.3e gpu=%+.3e)\n",
                name, rmax, j_abs, s_abs, cpu[j_abs, s_abs], gpu[j_abs, s_abs],
                rrel, j_rel, s_rel, cpu[j_rel, s_rel], gpu[j_rel, s_rel])
    end
    return outputs, phi_cpu, alpha_cpu, W_cpu, Np_cpu, Tp_cpu
end
#=
Legend:
    max abs    = max absolute difference between GPU and CPU results
    max rel    = max relative-to-|cpu| difference (with a 1e-8 floor on denom)
    bracket-ok = count of (j, s) where the kernel's sub-bracket scan found a
                 residual sign change in q1 and Brent's fixed iteration
                 committed to a phi; failed sections write zero and are
                 masked out of the max-abs/max-rel report.

Two experiments below:
  n_alpha sweep    — isolates linear-vs-Akima airfoil-table interp error
  n_iters sweep    — isolates Brent root-find convergence rate
=#

# ---------------------------------------------------------------------------
# Residual snapshot diagnostic (enabled when DIAGNOSE=true).
#
# For each entry in `problem_points`, we evaluate the BEMT residual at
# `n_snap` phi values across q1 = (eps, pi/2 - eps), on both the GPU (in
# DevFloat) and the CPU (Float64 reference). Printed side-by-side so we can
# see (a) where the residual crosses zero, (b) whether it's noisy/nonfinite
# at some phi values on the GPU that the CPU sees as smooth. That tells us
# whether Float32 precision is the culprit or whether the sub-bracket scan
# is missing a sign change.
#
# We reuse the batch inputs (Vx, Vy, pitch, blade_gpu) so the inputs match
# the ones that produced the anomalies in the main comparison above.
# ---------------------------------------------------------------------------

# One thread per phi sample, one (j, s) per launch. Writes residuals to a
# 1D device vector `out`.
using KernelAbstractions: @kernel, @index, get_backend, synchronize
@kernel function snapshot_kernel!(out, phis,
        Vx_all, Vy_all, pitch_all,
        r_vec, c_vec, twist_vec, rhub, rtip,
        alpha_min, dalpha, n_alpha, cl_table, cd_table,
        rho, B_blades, tip_mode, j_target::Int32, s_target::Int32)
    i = @index(Global)
    Vx = Vx_all[j_target, s_target]
    Vy = Vy_all[j_target, s_target]
    pt = pitch_all[s_target]
    r     = r_vec[j_target]
    chord = c_vec[j_target]
    twist = twist_vec[j_target]
    phi = phis[i]
    R = WATT.bemt_residual_and_outputs(
            phi, j_target, Vx, Vy, pt,
            r, chord, twist, rhub, rtip,
            rho,
            alpha_min, dalpha, n_alpha, cl_table, cd_table,
            B_blades, tip_mode)[1]
    out[i] = R
end

# Evaluate the CPU-side residual (Float64) at the same phis, using CCBlade's
# residual_and_outputs directly so we don't repeat WATT.solve_BEMT's quadrant
# selection logic.
function cpu_residuals(rotor_cpu, j, s, phis)
    airfoil = blade.airfoils[j]
    r_j    = blade.r[j]
    chord  = blade.c[j]
    twist  = blade.twist[j]
    xv = zeros(11)
    xv[1]=r_j; xv[2]=chord; xv[3]=twist; xv[4]=blade.rhub; xv[5]=blade.rtip
    xv[6]=Vx_h[j,s]; xv[7]=Vy_h[j,s]; xv[8]=env.rho; xv[9]=pitch_h[s]
    xv[10]=env.mu; xv[11]=env.a
    pv = (airfoil, rotor_cpu.B, rotor_cpu.turbine, rotor_cpu.re,
          rotor_cpu.mach, rotor_cpu.rotation, rotor_cpu.tip)
    return [WATT.CCBlade.residual_and_outputs(phi, xv, pv)[1] for phi in phis]
end

# Given a (j, s) and a blade_gpu, run the on-device snapshot kernel and
# pull the residual samples back to host.
function gpu_residuals(rotor_gpu, blade_gpu, j, s, phis)
    n_snap = length(phis)
    phis_dev = WATT.to_backend_vector(DevArrayType, DevFloat.(phis))
    out_dev  = WATT.to_backend_vector(DevArrayType, zeros(DevFloat, n_snap))
    backend  = get_backend(Vx)
    snapshot_kernel!(backend)(
        out_dev, phis_dev,
        Vx, Vy, pitch,
        blade_gpu.r, blade_gpu.c, blade_gpu.twist, blade_gpu.rhub, blade_gpu.rtip,
        blade_gpu.alpha_min, blade_gpu.dalpha, blade_gpu.n_alpha,
        blade_gpu.cl_table, blade_gpu.cd_table,
        DevFloat(env.rho),
        rotor_gpu.B, rotor_gpu.tip_mode,
        Int32(j), Int32(s);
        ndrange = n_snap,
    )
    synchronize(backend)
    return Float64.(to_host(out_dev))
end

# Sample q1 uniformly and print (phi, cpu_R, gpu_R). Because the kernel's
# residual is called with DevFloat, we also print `gpu_R - cpu_R` — a
# floor of ~1e-6 * |cpu_R| would just be Float32 noise, but wild
# divergences (sign flip, order-of-magnitude difference) mean Float32 is
# actually shifting where the root lives.
function residual_snapshot(label, rotor_cpu, rotor_gpu, blade_gpu, j, s;
                            n_snap::Int=21)
    phis = collect(range(1e-6, pi/2 - 1e-6, length=n_snap))
    Rcpu = cpu_residuals(rotor_cpu, j, s, phis)
    Rgpu = gpu_residuals(rotor_gpu, blade_gpu, j, s, phis)
    println("\n--- residual snapshot [$label]  j=$j  s=$s  ",
            "(U=$(round(Us[s], digits=2)) m/s, ",
            "TSR=$(round(Omegas[s]*rtip/Us[s], digits=2)), ",
            "pitch=$(round(pitch_h[s]*180/pi, digits=2))°, ",
            "r=$(round(blade.r[j], digits=2))m, ",
            "r/R=$(round(blade.r[j]/rtip, digits=3))) ---")
    @printf("  %-3s  %-10s  %-14s  %-14s  %s\n", "i", "phi (rad)", "cpu R (F64)", "gpu R ($(DevFloat))", "diff")
    for i in 1:n_snap
        # print with sign so sign changes are obvious
        d = Rgpu[i] - Rcpu[i]
        @printf("  %-3d  %+.4e  %+.6e  %+.6e  %+.3e\n", i, phis[i], Rcpu[i], Rgpu[i], d)
    end
    # Report where each series changes sign, so we can see if the sub-bracket
    # scan (which uses NSUB=10 subdivisions internally) would have latched
    # onto the same sub-interval on both.
    sign_changes(v) = findall(i -> v[i] * v[i+1] < 0, 1:length(v)-1)
    println("  sign-change indices — cpu: ", sign_changes(Rcpu), "   gpu: ", sign_changes(Rgpu))
end

# ---------------------------------------------------------------------------
# n_alpha sweep — hold n_iters at 20 (Brent's converged for this rotor),
# vary the airfoil table resolution. As n_alpha grows the linear interp
# converges to the Akima spline, so the Np error narrows toward the
# root-finding + Float32 floor.
# ---------------------------------------------------------------------------
println("\n### n_alpha sweep (no tip loss, n_iters=20)")
for nα in (361, 721, 1441, 2881, 5761)
    blade_gpu_nα = WATT.BladeGPU(blade; n_alpha=nα, ArrayType=DevArrayType)
    compare("n_alpha=$nα", rotor_notip, rotorgpu_notip, blade_gpu_nα, 20)
end

# ---------------------------------------------------------------------------
# n_iters sweep — hold n_alpha at 721, vary the Brent iteration count.
# If numbers stop changing above some n_iters, Brent has converged for
# this precision. Ran both with and without tip loss.
# ---------------------------------------------------------------------------
println("\n### Iteration-count sanity (n_alpha=721)")
blade_gpu = WATT.BladeGPU(blade; n_alpha=721, ArrayType=DevArrayType)
for nit in (10, 15, 20, 25, 30)
    compare("no tip loss", rotor_notip, rotorgpu_notip, blade_gpu, nit)
end
println()
for nit in (10, 15, 20, 25, 30)
    compare("Prandtl tip+hub", rotor_tip, rotorgpu_tip, blade_gpu, nit)
end

# ---------------------------------------------------------------------------
# Residual snapshot diagnostic (only if DIAGNOSE=true).
#
# These (j, s) points were the worst-case discrepancy locations in the
# earlier GPU-Float32 run:
#   (1,  1)   — innermost root cylinder; CPU reports phi=1.19 rad, GPU
#               drifted to ~7e-5 (essentially phi_lo).
#   (8, 24)   — inner DU40 section at high TSR / low wind; CPU alpha
#               near zero, GPU alpha = -0.27 rad (-15°). Big divergence
#               at a section that is nominally well-behaved.
#   (20, 108) — outer NACA section at high load; CPU Np = +16 kN/m,
#               GPU returned ~0. Bracket_ok was true, so the kernel
#               thought it found a solution.
# ---------------------------------------------------------------------------
if DIAGNOSE
    println("\n\n### Residual snapshots at problem (j, s) locations")
    problem_points = ((1, 1), (8, 24), (20, 108))
    for (j, s) in problem_points
        residual_snapshot("no tip loss", rotor_notip, rotorgpu_notip, blade_gpu, j, s)
    end
    for (j, s) in problem_points
        residual_snapshot("Prandtl tip+hub", rotor_tip, rotorgpu_tip, blade_gpu, j, s)
    end
end
