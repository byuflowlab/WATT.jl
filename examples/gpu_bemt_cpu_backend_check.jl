#=
CPU-backend correctness check for GPU-BEMT (Step 1 of the plan's verification).

Runs `WATT.solve_BEMT_gpu!` with a KernelAbstractions CPU backend and
compares its outputs against the CPU `WATT.solve_BEMT` (which wraps
CCBlade) across a sweep of operating conditions. Reports max/mean
absolute and relative differences on the quantities that matter to the
coupled sim (phi, Np, Tp, W, alpha).

Run with:
    julia --project examples/gpu_bemt_cpu_backend_check.jl
=#

using WATT, OpenFASTTools, DynamicStallModels
using StaticArrays, StructArrays
using Printf

const _of = OpenFASTTools
const _DS = DynamicStallModels
const PrandtlTipHub = WATT.CCBlade.PrandtlTipHub

# --- Build the NREL 5MW rotor (same setup as bracket-check script) ---------
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
n    = length(rvec)

airfoils = StructArray{_DS.Airfoil}(undef, n)
xcp = Vector{Float64}(undef, n)
for i = 1:n
    airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
end

blade = WATT.Blade(rvec, chordvec, twistvec, xcp, airfoils; rhub, rtip)

# Compare with-and-without tip loss to exercise both paths.
rotor_notip = WATT.Rotor(3, 80.0, true)                    # tip=nothing
rotor_tip   = WATT.Rotor(3, 80.0, true; tip=PrandtlTipHub())

env = WATT.environment(1.225, 1.464e-5, 343.0, 10.0, 1.0, 0.0)  # last two args: U, Omega — only rho/mu/a used

# --- Build the GPU-BEMT side (CPU backend) ----------------------------------
rotorgpu_notip = WATT.RotorGPU(rotor_notip)
rotorgpu_tip   = WATT.RotorGPU(rotor_tip)

# --- Build a batch of operating points --------------------------------------
wind_speeds = [5.0, 8.0, 11.4, 15.0, 20.0]      # m/s
tsrs        = [4.0, 6.0, 7.5, 9.0, 11.0]
pitches_deg = [-1.0, 0.0, 2.0, 5.0, 10.0]

# Cartesian product -> one sim per combination
Us      = Float64[]
Omegas  = Float64[]
pitches = Float64[]
for U in wind_speeds, tsr in tsrs, p in pitches_deg
    push!(Us, U)
    push!(Omegas, tsr * U / rtip)
    push!(pitches, p * pi/180)
end
n_sims = length(Us)

Vx = zeros(Float64, n, n_sims)
Vy = zeros(Float64, n, n_sims)
for s in 1:n_sims, j in 1:n
    Vx[j, s] = Us[s]
    Vy[j, s] = Omegas[s] * blade.r[j]
end
pitch_vec = collect(pitches)

# --- Run each rotor variant, compare -----------------------------------------
function compare(label, rotor_cpu, rotor_gpu, blade_gpu, n_iters)
    outputs = WATT.GPUBEMTOutputs(n, n_sims; ArrayType=Array{Float64})
    WATT.solve_BEMT_gpu!(outputs, rotor_gpu, blade_gpu, env, Vx, Vy, pitch_vec; n_iters=n_iters)

    # CPU reference
    phi_cpu = zeros(n, n_sims); Np_cpu = similar(phi_cpu); Tp_cpu = similar(phi_cpu)
    W_cpu   = similar(phi_cpu); alpha_cpu = similar(phi_cpu)
    xv = zeros(11)
    for s in 1:n_sims, j in 1:n
        ccout = WATT.solve_BEMT(rotor_cpu, blade, env, j, Vx[j, s], Vy[j, s], pitch_vec[s], xv)
        phi_cpu[j, s]   = ccout.phi
        Np_cpu[j, s]    = ccout.Np
        Tp_cpu[j, s]    = ccout.Tp
        W_cpu[j, s]     = ccout.W
        alpha_cpu[j, s] = ccout.alpha
    end

    println("\n=== $label  (n_iters=$n_iters, tip=$(rotor_cpu.tip)) ===")
    n_ok = count(outputs.success)
    n_bad = length(outputs.success) - n_ok
    println("bracket-ok: $n_ok / $(length(outputs.success))   (bracket-failed: $n_bad)")

    # Restrict comparison to sections where the GPU actually found a bracket;
    # failed sections write zeros (matches CPU CCBlade.Outputs()) but we don't
    # want NaN/zero blowups poisoning the max on the tip section.
    mask = outputs.success

    for (name, gpu, cpu) in (
            ("phi",   outputs.phi,   phi_cpu),
            ("alpha", outputs.alpha, alpha_cpu),
            ("W",     outputs.W,     W_cpu),
            ("Np",    outputs.Np,    Np_cpu),
            ("Tp",    outputs.Tp,    Tp_cpu),
        )
        d    = abs.(gpu .- cpu) .* mask
        rmax, imax = findmax(d)
        denom = max.(abs.(cpu), 1e-8)
        rrel_arr = (d ./ denom) .* mask
        rrel, irrel = findmax(rrel_arr)
        j_abs, s_abs = Tuple(imax)
        j_rel, s_rel = Tuple(irrel)
        @printf("  %-6s  max abs %.3e @ (j=%d,s=%d, cpu=%+.3e gpu=%+.3e)  max rel %.3e @ (j=%d,s=%d, cpu=%+.3e gpu=%+.3e)\n",
                name, rmax, j_abs, s_abs, cpu[j_abs, s_abs], gpu[j_abs, s_abs],
                rrel, j_rel, s_rel, cpu[j_rel, s_rel], gpu[j_rel, s_rel])
    end
end

println("\n### n_alpha sweep (no tip loss, n_iters=20) — check whether the ~1% Np gap")
println("### is dominated by the linear-vs-Akima airfoil interpolation:")
for nα in (361, 721, 1441, 2881, 5761)
    blade_gpu_nα = WATT.BladeGPU(blade; n_alpha=nα, ArrayType=Array{Float64})
    compare("n_alpha=$nα", rotor_notip, rotorgpu_notip, blade_gpu_nα, 20)
end

println("\n### Iteration-count sanity (n_alpha=721):")
blade_gpu = WATT.BladeGPU(blade; n_alpha=721, ArrayType=Array{Float64})
for nit in (10, 15, 20, 25, 30)
    compare("no tip loss", rotor_notip, rotorgpu_notip, blade_gpu, nit)
end
println()
for nit in (10, 15, 20, 25, 30)
    compare("Prandtl tip+hub", rotor_tip, rotorgpu_tip, blade_gpu, nit)
end
