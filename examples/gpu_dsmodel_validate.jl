#=
Validation for the GPU dynamic-stall port.

Loads the committed golden trace (examples/data/gpu_dsmodel_reference.jld2,
regenerate with gpu_dsmodel_reference.jl) and runs two gates:

  GATE 1 — device correctness. Marches the batch on the selected backend and
           again on the plain CPU (`Array`) backend. Both paths call the
           *identical* `ds_*_adg` functions, so any difference is pure
           backend / floating-point behavior. On :cpu this is exactly 0; on
           :cuda it isolates GPU transcendental/FP differences.

  GATE 2 — logic fidelity. Compares the CPU march against the Akima golden
           trace. The residual is the linear-table-vs-Akima interpolation
           error (largest near stall where the polar is sharply curved), not a
           port bug. Boolean states (LESF/TESF/VRTX) may flip at threshold
           ties; those are counted separately, not treated as continuous error.

Run:
    julia -t auto --project=examples examples/gpu_dsmodel_validate.jl        # :cpu
    (set BACKEND = :cuda on the cluster)

Adam Cardoza
=#

const BACKEND  = :cpu       # :cpu | :cuda | :amdgpu | :metal
const DEVFLOAT = Float64
const N_ALPHA  = 4001       # airfoil-table resolution (host-side; fidelity knob)

include(joinpath(@__DIR__, "gpu_dsmodel_common.jl"))
using JLD2, Printf

# ---------------------------------------------------------------------------
# Backend glue (would live in ext/WATTCUDAExt.jl once frozen). The DS structs
# only need host->device array transfer; the generic fallbacks in bemt_gpu.jl /
# dsmodel_gpu.jl already cover `AT(x)`, so this just picks the array type.
# ---------------------------------------------------------------------------
if BACKEND === :cuda
    using CUDA
    const DevArrayType = CuArray{DEVFLOAT}
elseif BACKEND === :amdgpu
    using AMDGPU
    const DevArrayType = ROCArray{DEVFLOAT}
elseif BACKEND === :metal
    using Metal
    const DevArrayType = MtlArray{DEVFLOAT}
else
    const DevArrayType = Array{DEVFLOAT}
end
const CpuArrayType = Array{DEVFLOAT}

@info "DS validation" BACKEND DevArrayType N_ALPHA

# ---------------------------------------------------------------------------
# Load reference + build blade
# ---------------------------------------------------------------------------
isfile(REF_FILE) || error("Reference $REF_FILE not found — run gpu_dsmodel_reference.jl first.")
@load REF_FILE U aoa dt xds_hist Cl Cd Cm n_sections n_sims nt scenario_names
@printf("Loaded reference: %d sections × %d sims × %d steps\n", n_sections, n_sims, nt)

blade, _ = build_nrel5mw_blade()

# March on a given ArrayType and return host copies.
function march(ArrayType)
    dsaf = WATT.DSAirfoilGPU(blade; n_alpha=N_ALPHA, ArrayType=ArrayType)
    Ud   = WATT.to_backend_array(ArrayType, DEVFLOAT.(U))
    aoad = WATT.to_backend_array(ArrayType, DEVFLOAT.(aoa))
    hist = WATT.DSHistory(n_sections, n_sims, nt; ArrayType=ArrayType)
    WATT.march_ds_gpu!(hist, dsaf, Ud, aoad, dt)
    return Array(hist.xds), Array(hist.Cl), Array(hist.Cd), Array(hist.Cm)
end

# ---------------------------------------------------------------------------
# GATE 1 — device vs CPU (same code)
# ---------------------------------------------------------------------------
xds_dev, Cl_dev, Cd_dev, Cm_dev = march(DevArrayType)
xds_cpu, Cl_cpu, Cd_cpu, Cm_cpu = (BACKEND === :cpu) ?
    (xds_dev, Cl_dev, Cd_dev, Cm_dev) : march(CpuArrayType)

g1_cont = maximum(abs.(xds_dev[1:28, :, :, :] .- xds_cpu[1:28, :, :, :]))
g1_bool = count(abs.(xds_dev[29:32, :, :, :] .- xds_cpu[29:32, :, :, :]) .> 0.5)
g1_load = max(maximum(abs.(Cl_dev .- Cl_cpu)), maximum(abs.(Cd_dev .- Cd_cpu)),
              maximum(abs.(Cm_dev .- Cm_cpu)))

println("\n=== GATE 1: device ($(BACKEND)) vs CPU, identical ds_*_adg code ===")
@printf("continuous states max|Δ| = %.3e\n", g1_cont)
@printf("boolean-state flips       = %d / %d\n", g1_bool, length(xds_dev[29:32, :, :, :]))
@printf("loads max|Δ|              = %.3e\n", g1_load)
g1_pass = g1_cont < 1e-9 && g1_load < 1e-9
println(g1_pass ? "GATE 1 PASS (device matches CPU to ~machine precision)" :
                  (BACKEND === :cpu ? "GATE 1: (cpu backend — expected exact)" :
                                      "GATE 1 WARN: differences exceed 1e-9 — inspect"))

# ---------------------------------------------------------------------------
# GATE 2 — CPU march vs Akima golden trace (interpolation fidelity)
# ---------------------------------------------------------------------------
g2_cont = maximum(abs.(xds_cpu[1:28, :, :, :] .- xds_hist[1:28, :, :, :]))
g2_bool = count(abs.(xds_cpu[29:32, :, :, :] .- xds_hist[29:32, :, :, :]) .> 0.5)
g2_cl = maximum(abs.(Cl_cpu .- Cl)); g2_cd = maximum(abs.(Cd_cpu .- Cd)); g2_cm = maximum(abs.(Cm_cpu .- Cm))

println("\n=== GATE 2: CPU march vs Akima reference (n_alpha=$N_ALPHA) ===")
@printf("continuous states max|Δ| = %.3e\n", g2_cont)
@printf("boolean-state flips       = %d / %d  (%.4f%%)\n",
        g2_bool, length(xds_cpu[29:32, :, :, :]), 100 * g2_bool / length(xds_cpu[29:32, :, :, :]))
@printf("loads max|Δ|  Cl=%.3e  Cd=%.3e  Cm=%.3e\n", g2_cl, g2_cd, g2_cm)
println("  per-sim max|ΔCl|:")
for s in 1:n_sims
    @printf("    %-14s %.3e\n", scenario_names[s], maximum(abs.(Cl_cpu[:, s, :] .- Cl[:, s, :])))
end
g2_load_tol = 5e-3
println((g2_cl < g2_load_tol && g2_cd < g2_load_tol && g2_cm < g2_load_tol) ?
        "GATE 2 PASS (loads within $g2_load_tol — interpolation-limited, port faithful)" :
        "GATE 2 WARN: loads exceed $g2_load_tol — inspect port logic")
