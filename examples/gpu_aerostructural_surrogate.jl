#=
GPU-resident coupled aeroelastic simulation with a surrogate structural model.

Demonstrates `run_sim_surrogate_gpu!`: the GPU-batched BEMT + dynamic-stall aero
kernels coupled to a batched NeuralKoopman structural surrogate, with all state
resident on the device across the time march.

Two demonstrations:

  (1) THREE-WAY OVERLAY (single design). Overlays, on the same axes:
        * the original full aeroelastic model  — GXBeam `run_sim!`
        * the CPU surrogate                    — `run_sim_surrogate!`  (ConditionedKoopman)
        * the GPU surrogate                    — `run_sim_surrogate_gpu!` (BatchedKoopman)
      for tip out-of-plane deflection, final spanwise aero loads, root internal
      moment, and tip internal force.

  (2) BATCHED MULTI-DESIGN RUN. Steps `ns` designs at once on the device (each
      with its own FiLM conditioning column), and reports throughput vs `ns`.

Backend switch at the top: `:cpu` runs the KernelAbstractions CPU backend
(validation / laptop); `:cuda` runs on an NVIDIA GPU (cluster). The surrogate
math and coupling kernels are identical across backends.

NOTE: the NREL-5MW-5seg root cylinders are given a real airfoil here — the GPU
dynamic-stall model (v1) has no NoModel/cylinder path yet, and forcing
BeddoesLeishman on a true cylinder polar makes the CPU Akima model blow up.
Treat trajectories as an integration demonstration, not a validated structural
prediction (the surrogate's design vector is the training-set conditioning, not
a mapping from this specific blade).

Adam Cardoza
=#

const BACKEND  = :cpu       # :cpu | :cuda
const DEVFLOAT = Float64

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, JLD2, LinearAlgebra, FLOWMath, Printf
using Plots

include(joinpath(@__DIR__, "koopman_surrogate.jl"))
const of = OpenFASTTools
const DS = DynamicStallModels

if BACKEND === :cuda
    using CUDA
    const ArrayType = CuArray{DEVFLOAT}
else
    const ArrayType = Array{DEVFLOAT}
end
@info "GPU coupled surrogate demo" BACKEND ArrayType

datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")

# ---------------------------------------------------------------------------
# Blade / rotor / assembly / env  (NREL 5MW 5-seg)
# ---------------------------------------------------------------------------
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)
aftype_names = ("Cylinder1.dat","Cylinder2.dat","DU40_A17.dat","DU35_A17.dat",
                "DU30_A17.dat","DU25_A17.dat","DU21_A17.dat","NACA64_A17.dat")
aftypes = [of.read_airfoilinput(joinpath(ofpath,"airfoils",name)) for name in aftype_names]
af_idx  = of.integerfit(raf, afidx, rvec)
af_idx_real = [af_idx[i] <= 2 ? 6 : af_idx[i] for i in 1:nr]   # cylinders → DU25 (see header note)
afs = aftypes[af_idx_real]

dsairfoils = Vector{DS.Airfoil}(undef, nr); xcp = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
end
blade = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)
rotor = WATT.Rotor(3, 90.0, true; tilt=0.0, yaw=0.0)
assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)
env = WATT.environment(1.225, 1.837e-5, 335.0, 11.4, 1.5, 0.2)

na = length(blade.r); np = length(assembly.points); nelem = length(assembly.elements)

surr_path = joinpath(datadir, "model_new.jld2")
diag = JLD2.load(joinpath(datadir, "ae_ic_diagnostic_new.jld2"))
x_b = diag["x_b"]; x_norm = vec(x_b)

tvec = collect(0:0.05:10.0)
nt = length(tvec)

# ---------------------------------------------------------------------------
# (1) Three-way overlay — single design
# ---------------------------------------------------------------------------
println("\n=== GXBeam full aeroelastic model ===")
aero_gx, gxhistory, mesh_gx = WATT.initialize_sim(blade, assembly, tvec)
WATT.run_sim!(rotor, blade, mesh_gx, env, tvec, aero_gx, gxhistory)

println("=== CPU surrogate (ConditionedKoopman) ===")
ck = build_conditioned_koopman(surr_path, x_norm, assembly)
aero_cpu, surr_hist_cpu, mesh_cpu = WATT.initialize_sim_surrogate(blade, assembly, tvec)
WATT.run_sim_surrogate!(rotor, blade, mesh_cpu, env, tvec, aero_cpu, surr_hist_cpu, ck)

println("=== GPU surrogate (BatchedKoopman, ns=1, $BACKEND) ===")
bk1 = build_batched_koopman(surr_path, reshape(x_norm, :, 1), assembly; ArrayType=ArrayType)
gmesh, ghist = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, 1; ArrayType=ArrayType)
WATT.run_sim_surrogate_gpu!(gmesh, ghist, env, tvec, bk1)

# copy device history back to host
gFx = Array(ghist.Fx); gFy = Array(ghist.Fy)
gu  = Array(ghist.u);  gF  = Array(ghist.F); gM = Array(ghist.M)

# comparable series
tip_uz_gx   = [gxhistory[i].points[end].u[3]     for i in 1:nt]
tip_uz_cpu  = [surr_hist_cpu[i].points[end].u[3] for i in 1:nt]
tip_uz_gpu  = [gu[3, end, 1, i]                  for i in 1:nt]

rootMx_gx   = [-gxhistory[i].points[1].M[1]      for i in 1:nt]
rootMx_cpu  = [surr_hist_cpu[i].points[1].M[1]   for i in 1:nt]
rootMx_gpu  = [gM[1, 1, 1, i]                     for i in 1:nt]

tipFy_gx    = [gxhistory[i].points[end].F[2]      for i in 1:nt]
tipFy_cpu   = [surr_hist_cpu[i].points[end].F[2]  for i in 1:nt]
tipFy_gpu   = [gF[2, end, 1, i]                    for i in 1:nt]

# plots
p1 = plot(tvec, tip_uz_gx,  label="GXBeam",       lw=2, color=:black,
          xlabel="Time (s)", ylabel="Tip def. z (m)", title="Out-of-plane tip deflection")
plot!(p1, tvec, tip_uz_cpu, label="CPU surrogate", lw=2, color=:royalblue, ls=:dash)
plot!(p1, tvec, tip_uz_gpu, label="GPU surrogate", lw=2, color=:crimson,   ls=:dot)

p2 = plot(blade.r, aero_gx.Fx[end, :],  label="Fx GXBeam",  lw=2, color=:black,
          xlabel="Span (m)", ylabel="Force (N/m)", title="Final spanwise aero loads")
plot!(p2, blade.r, aero_gx.Fy[end, :],  label="Fy GXBeam",  lw=2, color=:gray)
plot!(p2, blade.r, gFx[:, 1, end],      label="Fx GPU",     lw=2, color=:crimson, ls=:dot)
plot!(p2, blade.r, gFy[:, 1, end],      label="Fy GPU",     lw=2, color=:orange,  ls=:dot)

p3 = plot(tvec, rootMx_gx,  label="GXBeam",        lw=2, color=:black,
          xlabel="Time (s)", ylabel="Root M_x (N·m)", title="Root internal moment (x)")
plot!(p3, tvec, rootMx_cpu, label="CPU surrogate", lw=2, color=:royalblue, ls=:dash)
plot!(p3, tvec, rootMx_gpu, label="GPU surrogate", lw=2, color=:crimson,   ls=:dot)

p4 = plot(tvec, tipFy_gx,   label="GXBeam",        lw=2, color=:black,
          xlabel="Time (s)", ylabel="Tip F_y (N)", title="Tip internal force (y)")
plot!(p4, tvec, tipFy_cpu,  label="CPU surrogate", lw=2, color=:royalblue, ls=:dash)
plot!(p4, tvec, tipFy_gpu,  label="GPU surrogate", lw=2, color=:crimson,   ls=:dot)

plt = plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 800),
           plot_title="GXBeam vs CPU-surrogate vs GPU-surrogate")
savefig(plt, joinpath(@__DIR__, "gpu_surrogate_threeway.png"))
println("Saved overlay to examples/gpu_surrogate_threeway.png")

# quick agreement summary (GPU vs CPU surrogate)
reld(a, b) = maximum(abs.(a .- b)) / (maximum(abs.(a)) + eps())
@printf("\nGPU-vs-CPU-surrogate:  tip u_z rel=%.2e   root M_x rel=%.2e   tip F_y rel=%.2e\n",
        reld(tip_uz_cpu, tip_uz_gpu), reld(rootMx_cpu, rootMx_gpu), reld(tipFy_cpu, tipFy_gpu))

# ---------------------------------------------------------------------------
# (2) Batched multi-design run + throughput
# ---------------------------------------------------------------------------
println("\n=== Batched multi-design throughput ($BACKEND) ===")
for ns in (1, 8, 64, 256)
    # each design perturbs the conditioning vector slightly
    xb = reshape(x_norm, :, 1) .+ 0.01 .* randn(length(x_norm), ns)
    bk = build_batched_koopman(surr_path, xb, assembly; ArrayType=ArrayType)
    gm, gh = WATT.initialize_sim_surrogate_gpu(blade, rotor, assembly, env, tvec, ns; ArrayType=ArrayType)
    WATT.run_sim_surrogate_gpu!(gm, gh, env, tvec, bk)          # warm-up (compile)
    t0 = time_ns()
    WATT.run_sim_surrogate_gpu!(gm, gh, env, tvec, bk)
    dt = (time_ns() - t0) / 1e9
    @printf("  ns=%4d :  %.3f s   (%.1f sim-steps/s)\n", ns, dt, ns * nt / dt)
end

println("\nDone.")
