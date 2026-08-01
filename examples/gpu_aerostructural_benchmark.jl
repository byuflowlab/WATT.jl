#=
Coupled aeroelastic (surrogate) benchmark ladder for the supercomputer.

Compares five configurations of the WATT coupled transient solver so you can
see where each layer of GPU offloading pays off. Aero = BEMT + dynamic stall.

  A  original aero (CPU) + GXBeam structure          run_sim!                (full physics, no surrogate)
  B  original aero (CPU) + surrogate (CPU)            run_sim_surrogate!      (ConditionedKoopman)
  C  original aero (CPU) + surrogate (GPU)            run_sim_surrogate!      (HostBridgeKoopman → device)
  D  GPU-ported aero + surrogate, CPU backend         run_sim_surrogate_gpu!  (ArrayType = Array)
  E  GPU-ported aero + surrogate, GPU                 run_sim_surrogate_gpu!  (ArrayType = CuArray/…)

A–C are inherently single-sim. D/E batch `ns` designs at once on the device.
Section 1 times A–E at ns=1 (the apples-to-apples single-turbine ladder);
Section 2 sweeps `ns` for D (CPU backend) and E (device) and reports speedup
vs the CPU surrogate baseline B run `ns` times.

Knobs:
- BACKEND  — :cpu | :cuda | :amdgpu | :metal   (use :cuda on the cluster)
- DEVFLOAT — Float32 | Float64  (Float64 matches the surrogate's precision)
- NT       — steps per sim
- N_SIMS_SWEEP — batch sizes for section 2. Trim the top end on device OOM
  (history footprint printed as the GB column).
- CPU_BACKEND_MAX_NS — cap for the D (GPU-code-on-CPU) column; it does not
  parallelize, so it is slow at large ns.

Run:
    julia -t auto --project=examples examples/gpu_aerostructural_benchmark.jl
    (set BACKEND = :cuda on the cluster)

NOTE: NREL-5MW-5seg root cylinders are given a real airfoil (the GPU dynamic-
stall model has no NoModel/cylinder path yet). Throughput benchmark, not a
validated structural prediction.

Adam Cardoza
=#

const BACKEND  = :cpu      # :cpu | :cuda | :amdgpu | :metal
const DEVFLOAT = Float64    # Float32 | Float64
const NT       = 201        # steps per sim (≈10 s at dt=0.05)
const DT       = 0.05
const N_SIMS_SWEEP = (1, 8, 32, 128) #, 512, 1024)
const CPU_BACKEND_MAX_NS = 8   # cap for the D column (GPU code on CPU backend)

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, JLD2, LinearAlgebra, FLOWMath, Printf
include(joinpath(@__DIR__, "koopman_surrogate.jl"))
const of = OpenFASTTools
const DS = DynamicStallModels

# ---------------------------------------------------------------------------
# Backend glue. CpuArrayType is always Array (track D); DevArrayType is the
# selected device (tracks C, E). similar_type / to_backend_* are generic in src.
# ---------------------------------------------------------------------------
if BACKEND === :cuda
    using CUDA
    const DevArrayType = CuArray{DEVFLOAT}
    device_sync() = CUDA.synchronize()
elseif BACKEND === :amdgpu
    using AMDGPU
    const DevArrayType = ROCArray{DEVFLOAT}
    device_sync() = AMDGPU.synchronize()
elseif BACKEND === :metal
    using Metal
    const DevArrayType = MtlArray{DEVFLOAT}
    device_sync() = Metal.synchronize()
else
    const DevArrayType = Array{DEVFLOAT}
    device_sync() = nothing
end
const CpuArrayType = Array{DEVFLOAT}
@info "Coupled aeroelastic benchmark ladder" BACKEND DevArrayType DEVFLOAT NT

# ---------------------------------------------------------------------------
# Build blade / rotor / assembly / env (NREL 5MW 5-seg; cylinders → real airfoil).
# ---------------------------------------------------------------------------
datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)
aftype_names = ("Cylinder1.dat","Cylinder2.dat","DU40_A17.dat","DU35_A17.dat",
                "DU30_A17.dat","DU25_A17.dat","DU21_A17.dat","NACA64_A17.dat")
aftypes = [of.read_airfoilinput(joinpath(ofpath,"airfoils",name)) for name in aftype_names]
af_idx  = of.integerfit(raf, afidx, rvec)
af_idx_real = [af_idx[i] <= 2 ? 6 : af_idx[i] for i in 1:nr]
afs = aftypes[af_idx_real]
dsairfoils = Vector{DS.Airfoil}(undef, nr); xcp = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
end
blade = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)
rotor = WATT.Rotor(3, 90.0, true; tilt=0.0, yaw=0.0)
assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)

na = length(blade.r); np = length(assembly.points); nelem = length(assembly.elements)
tvec = collect(range(0.0, step=DT, length=NT))

# Fully turbulent inflow (TurbSim time series + power-law shear), matching
# aerostructural_nrel5mw5seg_surrogate.jl. env.U(t) is host-evaluated per step
# and uploaded to the device, so the GPU march needs no kernel changes. The
# tiling constructor fills the full NT·DT window (tiles the file if shorter).
rho = 1.225; mu_air = 1.837e-5; a_snd = 335.0; shearexp = 0.2
Vrated = 11.4; tsr = 7.55; omega = Vrated * tsr / (Rtip * cos(precone))
turbfile = joinpath(ofpath, "TurbSim.dat")
env = WATT.environment(turbfile, tvec[end], rho, mu_air, a_snd, omega, shearexp)
surr_path = joinpath(datadir, "model_new.jld2")
x_norm = vec(JLD2.load(joinpath(datadir, "ae_ic_diagnostic_new.jld2"))["x_b"])
nx = length(x_norm)

hist_bytes(ns) = (2 * na + 3 * 3 * np) * ns * NT * sizeof(DEVFLOAT)

# Mean wall time over `nsamp` timed runs after one warm-up, device-synced.
function timeit(f; nsamp=3)
    f(); device_sync()
    total = 0.0
    for _ in 1:nsamp
        t0 = time_ns(); f(); device_sync()
        total += (time_ns() - t0) / 1e9
    end
    return total / nsamp
end

# ===========================================================================
# Section 1 — single-sim ladder (ns = 1)
# ===========================================================================
println("\n=== Section 1: single-sim ladder (ns=1, NT=$NT) ===")

# A — GXBeam full aeroelastic (original aero, no surrogate)
aero_A, gxhist_A, mesh_A = WATT.initialize_sim(blade, assembly, tvec)
t_A = timeit(() -> WATT.run_sim!(rotor, blade, mesh_A, env, tvec, aero_A, gxhist_A); nsamp=2)

# B — original aero + surrogate (CPU)
ck = build_conditioned_koopman(surr_path, x_norm, assembly)
aero_B, hist_B, mesh_B = WATT.initialize_sim_surrogate(blade, assembly, tvec)
t_B = timeit(() -> WATT.run_sim_surrogate!(rotor, blade, mesh_B, env, tvec, aero_B, hist_B, ck); nsamp=3)

# C — original aero (CPU) + surrogate on device
hb = build_hostbridge_koopman(surr_path, x_norm, assembly; ArrayType=DevArrayType)
aero_C, hist_C, mesh_C = WATT.initialize_sim_surrogate(blade, assembly, tvec)
t_C = timeit(() -> WATT.run_sim_surrogate!(rotor, blade, mesh_C, env, tvec, aero_C, hist_C, hb); nsamp=3)

# D — GPU-ported aero + surrogate, CPU backend (Array)
bk_D = build_batched_koopman(surr_path, reshape(x_norm, :, 1), assembly; ArrayType=CpuArrayType)
gm_D, gh_D = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, 1; ArrayType=CpuArrayType)
t_D1 = timeit(() -> WATT.run_sim_surrogate_gpu!(gm_D, gh_D, env, tvec, bk_D); nsamp=3)

# E — GPU-ported aero + surrogate, device
bk_E = build_batched_koopman(surr_path, reshape(x_norm, :, 1), assembly; ArrayType=DevArrayType)
gm_E, gh_E = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, 1; ArrayType=DevArrayType)
t_E1 = timeit(() -> WATT.run_sim_surrogate_gpu!(gm_E, gh_E, env, tvec, bk_E); nsamp=3)

@printf("  %-4s %-52s %10s\n", "", "configuration", "time (s)")
@printf("  %-4s %-52s %10.4f\n", "A", "orig aero (CPU) + GXBeam structure",           t_A)
@printf("  %-4s %-52s %10.4f\n", "B", "orig aero (CPU) + surrogate (CPU)",            t_B)
@printf("  %-4s %-52s %10.4f\n", "C", "orig aero (CPU) + surrogate ($BACKEND)",       t_C)
@printf("  %-4s %-52s %10.4f\n", "D", "GPU-ported aero+surrogate, CPU backend",       t_D1)
@printf("  %-4s %-52s %10.4f\n", "E", "GPU-ported aero+surrogate, $BACKEND",          t_E1)
@printf("\n  surrogate speedup over GXBeam (B vs A): %.1fx\n", t_A / t_B)
@printf("  structural-offload speedup (C vs B):    %.2fx  (>1 means GPU helps at ns=1)\n", t_B / t_C)

# ===========================================================================
# Section 2 — batched GPU-stack throughput sweep (tracks D, E)
# ===========================================================================
println("\n=== Section 2: batched GPU-stack throughput ($BACKEND, $DEVFLOAT) ===")
@printf("na=%d np=%d nelem=%d nt=%d nx=%d   (baseline B: %.4f s/sim)\n", na, np, nelem, NT, nx, t_B)
@printf("%-7s | %-11s | %-12s | %-13s | %-12s | %-9s | %s\n",
        "ns", "E-GPU (s)", "E s/sim", "E Msteps/s", "E vs B×ns", "hist(GB)", "D-CPUbk (s)")
println(repeat('-', 100))
for ns in N_SIMS_SWEEP
    # Replicate the validated (in-distribution) conditioning across all sims.
    # The per-step compute is data-independent (fixed Brent iterations, fixed DS
    # recurrence, dense Koopman matmul), so throughput is identical to distinct
    # designs — but random off-manifold perturbations can make the surrogate
    # diverge to NaN over the horizon, which is a stability artifact, not a
    # timing signal. Swap in real design vectors here if you have a stable set.
    xb = repeat(reshape(x_norm, :, 1), 1, ns)

    bk = build_batched_koopman(surr_path, xb, assembly; ArrayType=DevArrayType)
    gm, gh = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, ns; ArrayType=DevArrayType)
    t_gpu = timeit(() -> WATT.run_sim_surrogate_gpu!(gm, gh, env, tvec, bk); nsamp=3)

    d_str = "—"
    if ns <= CPU_BACKEND_MAX_NS
        bkd = build_batched_koopman(surr_path, xb, assembly; ArrayType=CpuArrayType)
        gmd, ghd = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, ns; ArrayType=CpuArrayType)
        t_d = timeit(() -> WATT.run_sim_surrogate_gpu!(gmd, ghd, env, tvec, bkd); nsamp=1)
        d_str = @sprintf("%.3f", t_d)
    end

    work = na * ns * (NT - 1)
    @printf("%-7d | %-11.4f | %-12.5f | %-13.2f | %-12.1f | %-9.3f | %s\n",
            ns, t_gpu, t_gpu/ns, work/t_gpu/1e6, ns*t_B/t_gpu, hist_bytes(ns)/2^30, d_str)
    GC.gc()
end

println("\nJulia threads available: ", Threads.nthreads())
println("(Run with `julia -t auto --project=examples examples/gpu_aerostructural_benchmark.jl`.)")
println("Trim N_SIMS_SWEEP if the history footprint (GB) exceeds device memory.")
