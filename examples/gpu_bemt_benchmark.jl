#=
GPU-BEMT throughput benchmark for the supercomputer.

Times `WATT.solve_BEMT_gpu!` against the CPU `WATT.solve_BEMT` loop for
n_sims ∈ {1, 5, 10, 20, 50}, at n_iters = 20 (the tested default). Uses
the NREL 5MW-style rotor loaded from data/openfast, with Prandtl tip+hub
loss.

Backend selection (top of file):
- BACKEND = :cuda     — use CUDA.jl / CuArray
- BACKEND = :cpu      — KernelAbstractions CPU backend (sanity-run locally)
- BACKEND = :amdgpu   — AMDGPU.jl / ROCArray  (glue below is stubbed)
- BACKEND = :metal    — Metal.jl / MtlArray  (glue below is stubbed)

The CUDA branch also adds three tiny glue methods that would normally live
in a package extension (WATT has no explicit ext yet). Once we're ready to
freeze the interface, promote them into `ext/WATTCUDAExt.jl`.

Run:
    julia --project examples/gpu_bemt_benchmark.jl
(Set BACKEND below to :cpu locally, :cuda on the cluster.)
=#

const BACKEND = :cuda   # <-- change to :cpu for local dry-run

using WATT, OpenFASTTools, DynamicStallModels
using StaticArrays, StructArrays
using BenchmarkTools
using Printf

const _of = OpenFASTTools
const _DS = DynamicStallModels

# ---------------------------------------------------------------------------
# Backend glue
# ---------------------------------------------------------------------------
if BACKEND === :cuda
    using CUDA
    const DevFloat = Float32
    const DevArrayType = CuArray{DevFloat}
    const DevBoolArrayType = CuArray{Bool}
    # Glue methods (would live in ext/WATTCUDAExt.jl)
    WATT.to_backend_vector(::Type{CuArray{T}}, v::AbstractVector) where {T} = CuArray{T}(v)
    WATT.to_backend_matrix(::Type{CuArray{T}}, m::AbstractMatrix) where {T} = CuArray{T}(m)
    WATT.similar_type(::Type{CuArray{T}}, ::Type{U}) where {T, U} = CuArray{U}
    device_sync() = CUDA.synchronize()
elseif BACKEND === :amdgpu
    using AMDGPU
    const DevFloat = Float32
    const DevArrayType = ROCArray{DevFloat}
    const DevBoolArrayType = ROCArray{Bool}
    WATT.to_backend_vector(::Type{ROCArray{T}}, v::AbstractVector) where {T} = ROCArray{T}(v)
    WATT.to_backend_matrix(::Type{ROCArray{T}}, m::AbstractMatrix) where {T} = ROCArray{T}(m)
    WATT.similar_type(::Type{ROCArray{T}}, ::Type{U}) where {T, U} = ROCArray{U}
    device_sync() = AMDGPU.synchronize()
elseif BACKEND === :metal
    using Metal
    const DevFloat = Float32
    const DevArrayType = MtlArray{DevFloat}
    const DevBoolArrayType = MtlArray{Bool}
    WATT.to_backend_vector(::Type{MtlArray{T}}, v::AbstractVector) where {T} = MtlArray{T}(v)
    WATT.to_backend_matrix(::Type{MtlArray{T}}, m::AbstractMatrix) where {T} = MtlArray{T}(m)
    WATT.similar_type(::Type{MtlArray{T}}, ::Type{U}) where {T, U} = MtlArray{U}
    device_sync() = Metal.synchronize()
else
    const DevFloat = Float64
    const DevArrayType = Array{DevFloat}
    const DevBoolArrayType = Array{Bool}
    device_sync() = nothing
end

@info "Backend selected" BACKEND DevArrayType DevFloat

# ---------------------------------------------------------------------------
# Build NREL 5MW rotor
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
twistvec = adblade["BlTwist"] .* (pi/180)
rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
n_sections = length(rvec)

airfoils = StructArray{_DS.Airfoil}(undef, n_sections)
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
@printf("n_sections = %d,  n_iters = 20,  backend = %s,  eltype = %s\n\n",
        n_sections, BACKEND, DevFloat)

@printf("%-7s | %-14s | %-14s | %-14s | %s\n",
        "n_sims", "GPU (ms)", "CPU-serial (ms)", "CPU-@threads (ms)", "speedup GPU vs serial")
println(repeat('-', 100))

for n_sims in (1, 5, 10, 20, 50)
    Vx_h, Vy_h, pitch_h = make_batch(n_sims)

    # ---- GPU path ----
    Vx     = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vx_h))
    Vy     = WATT.to_backend_matrix(DevArrayType, DevFloat.(Vy_h))
    pitch  = WATT.to_backend_vector(DevArrayType, DevFloat.(pitch_h))
    outputs = WATT.GPUBEMTOutputs(n_sections, n_sims; ArrayType=DevArrayType)
    t_gpu = @belapsed begin
        WATT.solve_BEMT_gpu!($outputs, $rotor_gpu, $blade_gpu, $env, $Vx, $Vy, $pitch; n_iters=20)
        device_sync()
    end samples=50 evals=1

    # ---- CPU serial ----
    xv = zeros(11)
    t_cpu_ser = @belapsed begin
        for s in 1:$n_sims, j in 1:$n_sections
            WATT.solve_BEMT($rotor, $blade, $env, j, $Vx_h[j, s], $Vy_h[j, s], $pitch_h[s], $xv)
        end
    end samples=50 evals=1

    # ---- CPU @threads across sims (each thread uses its own scratch xv) ----
    xvs = [zeros(11) for _ in 1:Threads.nthreads()]
    t_cpu_par = @belapsed begin
        Threads.@threads for s in 1:$n_sims
            for j in 1:$n_sections
                WATT.solve_BEMT($rotor, $blade, $env, j, $Vx_h[j, s], $Vy_h[j, s], $pitch_h[s], $xvs[Threads.threadid()])
            end
        end
    end samples=50 evals=1

    @printf("%-7d | %-14.3f | %-14.3f | %-14.3f | %.1fx\n",
            n_sims, 1e3*t_gpu, 1e3*t_cpu_ser, 1e3*t_cpu_par, t_cpu_ser / t_gpu)
end

println("\nJulia threads available: ", Threads.nthreads())
println("(Use `julia -t auto --project examples/gpu_bemt_benchmark.jl` on the cluster.)")
