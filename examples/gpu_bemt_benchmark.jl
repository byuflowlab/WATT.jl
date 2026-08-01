#=
GPU-BEMT throughput benchmark for the supercomputer.

Times `WATT.solve_BEMT_gpu!` against the CPU `WATT.solve_BEMT` loop across
a sweep of n_sims at n_iters = 20 (the tested default). Uses the NREL
5MW-style rotor from data/openfast with Prandtl tip+hub loss.

Same two top-level knobs as gpu_bemt_cpu_backend_check.jl:
- BACKEND  — :cpu | :cuda | :amdgpu | :metal
- DEVFLOAT — Float32 or Float64. Default Float64 to match the coupled
             sim's structural surrogate precision. On H200 / A100-class
             hardware the FP64/FP32 throughput ratio is ~1/2 so the
             perf cost is modest; on consumer GPUs it can be 32×.

The CUDA branch adds three glue methods that would normally live in a
package extension. Promote to `ext/WATTCUDAExt.jl` once the interface is
frozen.

Run:
    julia -t auto --project examples/gpu_bemt_benchmark.jl
=#

const BACKEND  = :cuda      # :cpu | :cuda | :amdgpu | :metal
const DEVFLOAT = Float64    # Float32 or Float64

using WATT, OpenFASTTools, DynamicStallModels
using StaticArrays
using BenchmarkTools
using Printf

const _of = OpenFASTTools
const _DS = DynamicStallModels

# ---------------------------------------------------------------------------
# Backend glue
# ---------------------------------------------------------------------------
if BACKEND === :cuda
    using CUDA
    const DevFloat = DEVFLOAT
    const DevArrayType = CuArray{DevFloat}
    # Glue methods (would live in ext/WATTCUDAExt.jl)
    WATT.to_backend_vector(::Type{CuArray{T}}, v::AbstractVector) where {T} = CuArray{T}(v)
    WATT.to_backend_matrix(::Type{CuArray{T}}, m::AbstractMatrix) where {T} = CuArray{T}(m)
    WATT.similar_type(::Type{CuArray{T}}, ::Type{U}) where {T, U} = CuArray{U}
    device_sync() = CUDA.synchronize()
elseif BACKEND === :amdgpu
    using AMDGPU
    const DevFloat = DEVFLOAT
    const DevArrayType = ROCArray{DevFloat}
    WATT.to_backend_vector(::Type{ROCArray{T}}, v::AbstractVector) where {T} = ROCArray{T}(v)
    WATT.to_backend_matrix(::Type{ROCArray{T}}, m::AbstractMatrix) where {T} = ROCArray{T}(m)
    WATT.similar_type(::Type{ROCArray{T}}, ::Type{U}) where {T, U} = ROCArray{U}
    device_sync() = AMDGPU.synchronize()
elseif BACKEND === :metal
    using Metal
    const DevFloat = DEVFLOAT
    const DevArrayType = MtlArray{DevFloat}
    WATT.to_backend_vector(::Type{MtlArray{T}}, v::AbstractVector) where {T} = MtlArray{T}(v)
    WATT.to_backend_matrix(::Type{MtlArray{T}}, m::AbstractMatrix) where {T} = MtlArray{T}(m)
    WATT.similar_type(::Type{MtlArray{T}}, ::Type{U}) where {T, U} = MtlArray{U}
    device_sync() = Metal.synchronize()
else
    const DevFloat = DEVFLOAT
    const DevArrayType = Array{DevFloat}
    device_sync() = nothing
end

@info "Backend selected" BACKEND DevArrayType DevFloat

# ---------------------------------------------------------------------------
# Build NREL 5MW rotor
# ---------------------------------------------------------------------------
ofpath = "/home/cardoza2/wattplay/WATT.jl/data/openfast"
adblade = _of.read_adblade("sn5_adblade.dat", ofpath)
edfile  = _of.read_edfile("sn5_EDfile.dat", ofpath)

aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat", "DU25_A17.dat", "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [_of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]
af_idx = Int.(adblade["BlAFID"])
afs    = aftypes[af_idx]

chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"] .* (pi/180)
rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
n_sections = length(rvec)

airfoils = Vector{_DS.Airfoil}(undef, n_sections)
xcp = Vector{Float64}(undef, n_sections)
for i = 1:n_sections
    airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
end

blade = WATT.Blade(rvec, chordvec, twistvec, xcp, airfoils; rhub, rtip)
rotor = WATT.Rotor(3, 80.0, true; tip=WATT.CCBlade.PrandtlTipHub())
env   = WATT.environment(1.225, 1.464e-5, 343.0, 10.0, 1.0, 0.0)

# ---------------------------------------------------------------------------
# Build a batch of representative operating points
# ---------------------------------------------------------------------------
function make_batch(n_sims)
    # cycle through a representative set of (wind, TSR, pitch) tuples
    Us     = Float64[]; Omegas = Float64[]; pitches = Float64[]
    wind_speeds = [5.0, 8.0, 11.4, 15.0, 20.0]
    tsrs        = [4.0, 6.0, 7.5, 9.0, 11.0]
    pitches_deg = [-1.0, 0.0, 2.0, 5.0, 10.0]
    k = 0
    while length(Us) < n_sims
        for U in wind_speeds, tsr in tsrs, p in pitches_deg
            length(Us) == n_sims && break
            push!(Us, U)
            push!(Omegas, tsr * U / rtip)
            push!(pitches, p * pi/180)
            k += 1
        end
    end

    Vx_h = zeros(Float64, n_sections, n_sims)
    Vy_h = zeros(Float64, n_sections, n_sims)
    for s in 1:n_sims, j in 1:n_sections
        Vx_h[j, s] = Us[s]
        Vy_h[j, s] = Omegas[s] * blade.r[j]
    end
    return Vx_h, Vy_h, collect(pitches)
end

# ---------------------------------------------------------------------------
# Build a batch that triggers q3 at root sections.
#
# The aerostructural bracket trace found that q3 fires when the root
# cylinder's effective Vy goes small or negative — which happens under
# structural feedback (blade flap-back during a gust reduces or reverses
# the in-plane velocity). The static sweep alone can't reproduce that
# because the geometric phi = atan(Vx/Vy) with Vx>0 and Vy>0 stays inside
# q1 no matter how small Vy gets.
#
# To force the q3 code path without dragging in a full coupled sim, this
# batch takes the nominal (Vx, Vy) matrix and overrides `Vy` on the six
# innermost aero nodes (the sections that actually hit q3 in the trace)
# to a small negative value. This mimics the structural-velocity kickback
# and pushes the root-section root into q3.
# ---------------------------------------------------------------------------
function make_batch_root_negVy(n_sims; n_root_sections=6, root_Vy=-0.5)
    Vx_h, Vy_h, pitch_h = make_batch(n_sims)
    for s in 1:n_sims, j in 1:min(n_root_sections, n_sections)
        Vy_h[j, s] = root_Vy
    end
    return Vx_h, Vy_h, pitch_h
end

# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------
blade_gpu  = WATT.BladeGPU(blade; n_alpha=721, ArrayType=DevArrayType)
rotor_gpu  = WATT.RotorGPU(rotor)

# Pre-populate JIT and warm the GPU
let
    Vx_h, Vy_h, pitch_h = make_batch(4)
    Vx = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vx_h))
    Vy = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vy_h))
    pitch = WATT.to_backend_vector(DevArrayType, DevFloat.(pitch_h))
    outputs = WATT.GPUBEMTOutputs(n_sections, 4; ArrayType=DevArrayType)
    WATT.solve_BEMT_gpu!(outputs, rotor_gpu, blade_gpu, env, Vx, Vy, pitch; n_iters=20)
    device_sync()
end

println("\n=== GPU-BEMT throughput benchmark ===")
@printf("n_sections = %d,  n_iters = 20,  backend = %s,  eltype = %s\n",
        n_sections, BACKEND, DevFloat)

# n_sims sweep. Each entry launches n_sections × n_sims threads on the GPU.
# NREL 5MW has 25 aero sections, so n_sims=200 → 5000 threads (H200 has ~16k
# CUDA cores; a full "wave" is ~5–10k threads depending on occupancy). Values
# larger than the coupled-sim use case are included to see where GPU perf
# actually flattens.
const N_SIMS_SWEEP = (1, 5, 10, 20, 50, 100) #, 200, 500, 1000) # I don't need to run past this since I won't ever be running that many sims.

# Runs one row of the benchmark for a given batch. `q3_count` = number of
# (section, sim) outputs where the converged phi landed in q3 = (pi/2, pi).
# On the GPU path this tells us whether the conditional q3 fallback fired.
function run_row(n_sims, batch_maker)
    Vx_h, Vy_h, pitch_h = batch_maker(n_sims)

    # ---- GPU path ----
    Vx     = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vx_h))
    Vy     = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vy_h))
    pitch  = WATT.to_backend_vector(DevArrayType, DevFloat.(pitch_h))
    outputs = WATT.GPUBEMTOutputs(n_sections, n_sims; ArrayType=DevArrayType)
    t_gpu = @belapsed begin
        WATT.solve_BEMT_gpu!($outputs, $rotor_gpu, $blade_gpu, $env, $Vx, $Vy, $pitch; n_iters=20)
        device_sync()
    end samples=50 evals=1

    # Count q3 landings in the GPU output for this batch.
    phi_host = Array(outputs.phi)
    q3_count = count(p -> p > pi/2, phi_host)

    # ---- CPU serial ----
    xv = zeros(11)
    t_cpu_ser = @belapsed begin
        for s in 1:$n_sims, j in 1:$n_sections
            WATT.solve_BEMT($rotor, $blade, $env, j, $Vx_h[j, s], $Vy_h[j, s], $pitch_h[s], $xv)
        end
    end samples=50 evals=1

    # ---- CPU @threads across sims (each thread uses its own scratch xv) ----
    # Size by `maxthreadid()` — since Julia 1.9 interactive/GC threads can
    # produce threadids above `nthreads(:default)`, so `[_ for _ in
    # 1:nthreads()]` under-allocates and BoundsErrors inside @threads.
    # `:static` scheduling forbids task migration, which keeps `threadid()`
    # stable inside a chunk.
    xvs = [zeros(11) for _ in 1:Threads.maxthreadid()]
    t_cpu_par = @belapsed begin
        Threads.@threads :static for s in 1:$n_sims
            for j in 1:$n_sections
                WATT.solve_BEMT($rotor, $blade, $env, j, $Vx_h[j, s], $Vy_h[j, s], $pitch_h[s], $xvs[Threads.threadid()])
            end
        end
    end samples=50 evals=1

    @printf("%-7d | %-14.3f | %-14.3f | %-14.3f | %.1fx | %d\n",
            n_sims, 1e3*t_gpu, 1e3*t_cpu_ser, 1e3*t_cpu_par, t_cpu_ser / t_gpu, q3_count)
end

function run_sweep(label, batch_maker)
    println("\n--- $label ---")
    @printf("%-7s | %-14s | %-14s | %-14s | %-6s | %s\n",
            "n_sims", "GPU (ms)", "CPU-serial (ms)", "CPU-@threads (ms)", "speedup", "q3 landings")
    println(repeat('-', 105))
    for n_sims in N_SIMS_SWEEP
        run_row(n_sims, batch_maker)
    end
end

# 1) Nominal operation — should stay 100% in q1, q3-landings column = 0.
#    Measures the baseline cost of the (now-conditional) q1+q3 kernel when
#    q3 never fires (i.e. same warps as before; the fallback branch is
#    masked out).
run_sweep("Nominal (all q1)", make_batch)

# 2) Root-heavy operation — forces q3 at the six innermost aero nodes.
#    Those sections cluster into a single warp, so the divergence cost of
#    the q3 fallback is confined to one warp per timestep. Comparing the
#    GPU column to the nominal sweep at the same n_sims measures that cost.
run_sweep("Root-heavy (forces q3 at innermost 6 sections)", make_batch_root_negVy)

println("\nJulia threads available: ", Threads.nthreads())
println("(Use `julia -t auto --project examples/gpu_bemt_benchmark.jl` on the cluster.)")
