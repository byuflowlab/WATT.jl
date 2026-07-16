#=
GPU dynamic-stall throughput benchmark for the supercomputer.

Companion to gpu_bemt_benchmark.jl, but for the batched Beddoes-Leishman v3
(ADG) march. Times `WATT.march_ds_gpu!` on the selected backend against a
native serial CPU march and a native `@threads` CPU march, sweeping n_sims at a
fixed step count. Uses the same NREL 5MW blade and prescribed (U, aoa) batch as
the validation/warp scripts (see gpu_dsmodel_common.jl).

Two baselines, mirroring the BEMT benchmark:
- CPU-serial / CPU-@threads use the *native* DynamicStallModels.jl routines
  (`DS.initialize`, `DS.update_states_ADG!`, `DS.get_loads`) marched per
  (section, sim) — the same code paths the reference generator uses and the
  coupled sim's `take_aero_step!` uses. This is the honest "what the CPU would
  actually run" baseline, independent of the GPU kernel.
- The GPU path runs the KernelAbstractions kernels through `march_ds_gpu!`.

Both paths write the full (32, n_sections, n_sims, nt) state history plus loads,
so the comparison is apples-to-apples on memory traffic.

Two top-level knobs, same as gpu_bemt_benchmark.jl:
- BACKEND  — :cpu | :cuda | :amdgpu | :metal
- DEVFLOAT — Float32 or Float64. Default Float64 to match the coupled sim's
             structural surrogate precision.

Run:
    julia -t auto --project=examples examples/gpu_dsmodel_benchmark.jl        # :cpu
    (set BACKEND = :cuda on the cluster)

Adam Cardoza
=#

const BACKEND  = :cuda      # :cpu | :cuda | :amdgpu | :metal
const DEVFLOAT = Float64    # Float32 or Float64
const N_ALPHA  = 721        # airfoil-table resolution. Table lookup is O(1)
                            # (index by alpha_min/dalpha), so this affects
                            # memory, not per-step cost — kept modest.
const NT       = 250        # steps per sim. March cost scales linearly with NT.
const DT       = 0.01

include(joinpath(@__DIR__, "gpu_dsmodel_common.jl"))
using BenchmarkTools, Printf

# ---------------------------------------------------------------------------
# Backend glue. The DS path only needs host->device array transfer, which the
# generic `to_backend_array(AT, a) = AT(a)` in dsmodel_gpu.jl already covers.
# So this block just picks the device array type and a synchronize hook.
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

@info "Backend selected" BACKEND DevArrayType DEVFLOAT N_ALPHA NT

# ---------------------------------------------------------------------------
# Build blade + DS airfoil tables (device side once; batch rebuilt per n_sims).
# ---------------------------------------------------------------------------
blade, _ = build_nrel5mw_blade()
n_sections = length(blade.r)
dsaf_dev = WATT.DSAirfoilGPU(blade; n_alpha=N_ALPHA, ArrayType=DevArrayType)

# Cycle the reference scenarios to reach an arbitrary n_sims (each gets a
# unique name so make_ds_batch keeps them distinct).
function make_scenarios(n_sims)
    base = DS_SCENARIOS
    scs = NamedTuple[]
    k = 0
    while length(scs) < n_sims
        sc = base[(k % length(base)) + 1]
        push!(scs, (name = "$(sc.name)_$k", aoa_mean = sc.aoa_mean,
                    aoa_amp = sc.aoa_amp, freq = sc.freq, Uscale = sc.Uscale))
        k += 1
    end
    return Tuple(scs)
end

# ---------------------------------------------------------------------------
# Native CPU march (independent of the KA kernel). Marches one 2D section
# through its prescribed (U, aoa) history, writing the (32 × nt) state history
# and (Cl, Cd, Cm) into the supplied views. Mirrors gpu_dsmodel_reference.jl.
# ---------------------------------------------------------------------------
function march_section!(xds_col, Cl_row, Cd_row, Cm_row, airfoil, c, xcp, Uts, aoats, dt)
    nt = length(Uts)
    p  = [c, xcp]
    y  = [Uts[1], 0.0, aoats[1], 0.0]

    x0, l0 = _DS.initialize(airfoil, [0.0], y, p)
    xds_col[:, 1] .= x0
    Cl_row[1] = l0[1]; Cd_row[1] = l0[2]; Cm_row[1] = l0[3]

    xold = collect(x0)
    xnew = zero(xold)
    for i in 2:nt
        y[1] = Uts[i]; y[2] = 0.0; y[3] = aoats[i]; y[4] = 0.0
        _DS.update_states_ADG!(airfoil, xold, xnew, y, p, dt)
        Cl, Cd, Cm = _DS.get_loads(airfoil.model, airfoil, xnew, y, p)
        xds_col[:, i] .= xnew
        Cl_row[i] = Cl; Cd_row[i] = Cd; Cm_row[i] = Cm
        xold, xnew = xnew, xold
    end
    return nothing
end

function march_batch_cpu_serial!(xds, Cl, Cd, Cm, blade, U, aoa, dt)
    n_sections, n_sims, _ = size(U)
    for s in 1:n_sims, j in 1:n_sections
        march_section!(view(xds, :, j, s, :), view(Cl, j, s, :), view(Cd, j, s, :),
                       view(Cm, j, s, :), blade.airfoils[j], blade.c[j], blade.xcp[j],
                       view(U, j, s, :), view(aoa, j, s, :), dt)
    end
    return nothing
end

function march_batch_cpu_threads!(xds, Cl, Cd, Cm, blade, U, aoa, dt)
    n_sections, n_sims, _ = size(U)
    # Each (section, sim) march is fully independent and writes disjoint slices,
    # and march_section! keeps all its scratch local, so no per-thread scratch
    # is needed. :static keeps threadid stable but we don't rely on it here.
    Threads.@threads :static for s in 1:n_sims
        for j in 1:n_sections
            march_section!(view(xds, :, j, s, :), view(Cl, j, s, :), view(Cd, j, s, :),
                           view(Cm, j, s, :), blade.airfoils[j], blade.c[j], blade.xcp[j],
                           view(U, j, s, :), view(aoa, j, s, :), dt)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Warm-up: JIT both paths and warm the GPU.
# ---------------------------------------------------------------------------
let
    U, aoa, _, dt = make_ds_batch(blade; scenarios = make_scenarios(4), nt = NT, dt = DT)
    Ud   = WATT.to_backend_array(DevArrayType, DEVFLOAT.(U))
    aoad = WATT.to_backend_array(DevArrayType, DEVFLOAT.(aoa))
    hist = WATT.DSHistory(n_sections, 4, NT; ArrayType = DevArrayType)
    WATT.march_ds_gpu!(hist, dsaf_dev, Ud, aoad, dt)
    device_sync()

    xds = Array{Float64}(undef, WATT.DS_NSTATES, n_sections, 4, NT)
    Cl = zeros(n_sections, 4, NT); Cd = similar(Cl); Cm = similar(Cl)
    march_batch_cpu_serial!(xds, Cl, Cd, Cm, blade, U, aoa, dt)
    march_batch_cpu_threads!(xds, Cl, Cd, Cm, blade, U, aoa, dt)
end

# ---------------------------------------------------------------------------
# Benchmark sweep over n_sims. Each GPU step launches n_sections × n_sims
# threads; the march queues NT such launches serially on the stream (step i
# reads step i-1). NREL 5MW has n_sections=25, so n_sims=100 → 2500 threads
# per launch. Values past the coupled-sim use case are included to see where
# GPU throughput flattens.
# ---------------------------------------------------------------------------
const N_SIMS_SWEEP = (1, 5, 10, 20, 50, 100)

function run_row(n_sims)
    U, aoa, _, dt = make_ds_batch(blade; scenarios = make_scenarios(n_sims), nt = NT, dt = DT)

    # ---- GPU (or CPU-backend KA) path ----
    Ud   = WATT.to_backend_array(DevArrayType, DEVFLOAT.(U))
    aoad = WATT.to_backend_array(DevArrayType, DEVFLOAT.(aoa))
    hist = WATT.DSHistory(n_sections, n_sims, NT; ArrayType = DevArrayType)
    t_gpu = @belapsed begin
        WATT.march_ds_gpu!($hist, $dsaf_dev, $Ud, $aoad, $dt)
        device_sync()
    end samples=20 evals=1

    # ---- native CPU serial ----
    xds = Array{Float64}(undef, WATT.DS_NSTATES, n_sections, n_sims, NT)
    Cl = zeros(n_sections, n_sims, NT); Cd = similar(Cl); Cm = similar(Cl)
    t_cpu_ser = @belapsed march_batch_cpu_serial!($xds, $Cl, $Cd, $Cm, $blade, $U, $aoa, $dt) samples=20 evals=1

    # ---- native CPU @threads across sims ----
    t_cpu_par = @belapsed march_batch_cpu_threads!($xds, $Cl, $Cd, $Cm, $blade, $U, $aoa, $dt) samples=20 evals=1

    # section-sim-steps processed, for a throughput figure.
    work = n_sections * n_sims * (NT - 1)
    @printf("%-7d | %-14.3f | %-16.3f | %-16.3f | %-7.1fx | %.2f\n",
            n_sims, 1e3*t_gpu, 1e3*t_cpu_ser, 1e3*t_cpu_par,
            t_cpu_ser / t_gpu, work / t_gpu / 1e6)
end

println("\n=== GPU dynamic-stall throughput benchmark ===")
@printf("n_sections = %d,  nt = %d,  backend = %s,  eltype = %s\n",
        n_sections, NT, BACKEND, DEVFLOAT)
@printf("%-7s | %-14s | %-16s | %-16s | %-8s | %s\n",
        "n_sims", "GPU (ms)", "CPU-serial (ms)", "CPU-@threads (ms)", "speedup", "GPU Msteps/s")
println(repeat('-', 100))
for n_sims in N_SIMS_SWEEP
    run_row(n_sims)
end

println("\nJulia threads available: ", Threads.nthreads())
println("(Use `julia -t auto --project=examples examples/gpu_dsmodel_benchmark.jl` on the cluster.)")
