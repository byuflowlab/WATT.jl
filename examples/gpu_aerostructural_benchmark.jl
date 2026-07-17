#=
GPU coupled aeroelastic (surrogate) throughput benchmark for the supercomputer.

Times the full GPU-resident coupled march `WATT.run_sim_surrogate_gpu!` — GPU
BEMT + dynamic stall + batched NeuralKoopman structural surrogate, all state on
device — sweeping the batch size `ns` (number of simultaneous designs) at a
fixed step count. Reports wall time, per-sim time, and section·sim·step
throughput, plus the device-memory footprint of the history buffers.

Baseline: the CPU surrogate solver `WATT.run_sim_surrogate!` (single sim,
ConditionedKoopman) timed once. Since it is inherently single-sim, the honest
"what the CPU would run for a batch" figure is `ns × t_cpu1` (sequential solves);
the reported speedup is `ns·t_cpu1 / t_gpu`.

Knobs:
- BACKEND  — :cpu | :cuda | :amdgpu | :metal   (use :cuda on the cluster)
- DEVFLOAT — Float32 | Float64  (Float64 matches the surrogate's precision)
- NT       — steps per sim (cost scales ~linearly)
- N_SIMS_SWEEP — batch sizes to sweep. Trim the top end if you hit device OOM:
  history is (na + 3·np·(3 forces))·ns·nt·8 bytes ≈ see the printed footprint.

Run:
    julia -t auto --project=examples examples/gpu_aerostructural_benchmark.jl
    (set BACKEND = :cuda on the cluster)

NOTE: the NREL-5MW-5seg root cylinders are given a real airfoil here (the GPU
dynamic-stall model has no NoModel/cylinder path yet). This is a throughput
benchmark, not a validated structural prediction.

Adam Cardoza
=#

const BACKEND  = :cuda      # :cpu | :cuda | :amdgpu | :metal
const DEVFLOAT = Float64    # Float32 | Float64
const NT       = 201        # steps per sim (≈10 s at dt=0.05)
const DT       = 0.05
const N_SIMS_SWEEP = (1, 8, 32, 128, 512, 1024)
const CPU_BASELINE_MAX_NS = 32   # also run the CPU solver serially for ns ≤ this

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, StructArrays, JLD2, LinearAlgebra, FLOWMath, Printf
include(joinpath(@__DIR__, "koopman_surrogate.jl"))
const of = OpenFASTTools
const DS = DynamicStallModels

# ---------------------------------------------------------------------------
# Backend glue (host->device transfer is covered by the generic to_backend_*).
# ---------------------------------------------------------------------------
if BACKEND === :cuda
    using CUDA
    const ArrayType = CuArray{DEVFLOAT}
    device_sync() = CUDA.synchronize()
elseif BACKEND === :amdgpu
    using AMDGPU
    const ArrayType = ROCArray{DEVFLOAT}
    device_sync() = AMDGPU.synchronize()
elseif BACKEND === :metal
    using Metal
    const ArrayType = MtlArray{DEVFLOAT}
    device_sync() = Metal.synchronize()
else
    const ArrayType = Array{DEVFLOAT}
    device_sync() = nothing
end
@info "GPU coupled surrogate benchmark" BACKEND ArrayType DEVFLOAT NT

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
dsairfoils = StructArray{DS.Airfoil}(undef, nr); xcp = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
end
blade = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)
rotor = WATT.Rotor(3, 90.0, true; tilt=0.0, yaw=0.0)
assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)
env = WATT.environment(1.225, 1.837e-5, 335.0, 11.4, 1.5, 0.2)

na = length(blade.r); np = length(assembly.points); nelem = length(assembly.elements)
tvec = collect(range(0.0, step=DT, length=NT))

surr_path = joinpath(datadir, "model_new.jld2")
x_norm = vec(JLD2.load(joinpath(datadir, "ae_ic_diagnostic_new.jld2"))["x_b"])
nx = length(x_norm)

# History device-memory footprint (bytes): Fx,Fy (na,ns,nt) + u,F,M (3,np,ns,nt).
hist_bytes(ns) = (2 * na + 3 * 3 * np) * ns * NT * sizeof(DEVFLOAT)

# Min wall time over `nsamp` timed runs (after one warm-up), device-synced.
function timeit(f; nsamp=3)
    f(); device_sync()                    # warm-up / compile
    best = Inf
    for _ in 1:nsamp
        t0 = time_ns(); f(); device_sync()
        best = min(best, (time_ns() - t0) / 1e9)
    end
    return best
end

# ---------------------------------------------------------------------------
# CPU single-sim baseline (ConditionedKoopman + run_sim_surrogate!).
# ---------------------------------------------------------------------------
println("\n=== CPU baseline: run_sim_surrogate! (1 sim) ===")
ck = build_conditioned_koopman(surr_path, x_norm, assembly)
aero_cpu, surr_hist_cpu, mesh_cpu = WATT.initialize_sim_surrogate(blade, assembly, tvec)
t_cpu1 = timeit(() -> WATT.run_sim_surrogate!(rotor, blade, mesh_cpu, env, tvec,
                                              aero_cpu, surr_hist_cpu, ck); nsamp=3)
@printf("  t_cpu(1 sim) = %.4f s  (%.1f sim-steps/s)\n", t_cpu1, NT / t_cpu1)

# ---------------------------------------------------------------------------
# GPU sweep.
# ---------------------------------------------------------------------------
println("\n=== GPU coupled surrogate throughput ($BACKEND, $DEVFLOAT) ===")
@printf("na=%d  np=%d  nelem=%d  nt=%d  nx=%d\n", na, np, nelem, NT, nx)
@printf("%-7s | %-11s | %-12s | %-13s | %-11s | %-9s | %s\n",
        "ns", "GPU (s)", "s / sim", "Msteps/s", "vs CPU×ns", "hist (GB)", "CPU-serial (s)")
println(repeat('-', 104))

for ns in N_SIMS_SWEEP
    # per-design conditioning: perturb x so each sim is a distinct design
    xb = reshape(x_norm, :, 1) .+ (ns == 1 ? 0.0 : 0.01) .* randn(nx, ns)
    bk = build_batched_koopman(surr_path, xb, assembly; ArrayType=ArrayType)
    gm, gh = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, ns; ArrayType=ArrayType)

    t_gpu = timeit(() -> WATT.run_sim_surrogate_gpu!(gm, gh, env, tvec, bk); nsamp=3)

    work = na * ns * (NT - 1)
    persim = t_gpu / ns
    speedup = ns * t_cpu1 / t_gpu
    gb = hist_bytes(ns) / 2^30

    # Optional actual CPU-serial batch (ns sequential solves) for small ns.
    cpu_ser_str = "—"
    if ns <= CPU_BASELINE_MAX_NS
        t_cpu_ser = timeit(() -> (for _ in 1:ns
            WATT.run_sim_surrogate!(rotor, blade, mesh_cpu, env, tvec, aero_cpu, surr_hist_cpu, ck)
        end); nsamp=1)
        cpu_ser_str = @sprintf("%.3f", t_cpu_ser)
    end

    @printf("%-7d | %-11.4f | %-12.5f | %-13.2f | %-11.1f | %-9.3f | %s\n",
            ns, t_gpu, persim, work / t_gpu / 1e6, speedup, gb, cpu_ser_str)
    GC.gc()
end

println("\nJulia threads available: ", Threads.nthreads())
println("(Run with `julia -t auto --project=examples examples/gpu_aerostructural_benchmark.jl`.)")
println("Trim N_SIMS_SWEEP if the history footprint (GB column) exceeds device memory.")
