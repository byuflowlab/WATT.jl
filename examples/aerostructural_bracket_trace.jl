#=
Aerostructural quadrant trace — Step 0 in a "real transient sim" flavor.

`gpu_bemt_bracket_check.jl` proved q1 is the winning bracket 100% of the time
across a static Cartesian sweep of (wind, TSR, pitch, section). This script
does the same check but from inside an *unsteady, two-way-coupled*
simulation, which is closer to what the GPU BEMT will actually see when it
plugs into the surrogate coupling. The residual sees non-trivial section
velocities (Vy shaped by structural deformation and rotation, Vx varied by
turbulent inflow) that the static sweep doesn't cover.

Mechanism: `src/bemt.jl` has a `BEMT_QUAD_TRACE` `Ref{Union{Nothing, Dict}}`
that solve_BEMT increments on each successful solve. It's `nothing` by
default (zero overhead); this script sets it to a fresh `Dict{Symbol,Int}`
before running and clears it afterwards, and prints the final tally.

Rotor + inflow setup mirrors examples/aerostructural_nrel5mw5seg.jl so the
sim is a real coupled unsteady case.

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, JLD2
using Printf

const of = OpenFASTTools
const DS = DynamicStallModels

datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")

# --- Load PreComp blade fixture (identical to aerostructural_nrel5mw5seg.jl) ---
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)

# --- Build dsairfoils ---
aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat",  "DU25_A17.dat",  "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]
af_idx  = of.integerfit(raf, afidx, rvec)
afs     = aftypes[af_idx]

dsairfoils = Vector{DS.Airfoil}(undef, nr)
xcp        = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    if cylinder_mask[i]
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; dsmodel=DS.NoModel(), polar=polars[i])
    else
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
    end
end

# --- Blade, Rotor, Assembly ---
blade = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)
B, hubHt = 3, 90.0
rotor = WATT.Rotor(B, hubHt, true; tilt=0.0, yaw=0.0)
assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)

# --- Turbulent inflow environment ---
rho, mu, a  = 1.225, 1.837e-5, 335.0
shearexp    = 0.2
Vrated, tsr = 11.4, 7.55
rotorR      = Rtip * cos(precone)
omega       = Vrated * tsr / rotorR
turbfile    = joinpath(ofpath, "TurbSim.dat")
env         = environment(turbfile, rho, mu, a, omega, shearexp)

# --- Simulate ---
# Longer window than the initial 5s sanity run so we sample more of the
# turbulent inflow record. 30s at 50 Hz = 601 steps × 36 nodes ≈ 22k solves,
# still cheap on CPU.
tvec = collect(0:0.05:30.0)
n_steps = length(tvec)

println("\n=== Aerostructural bracket trace ===")
@printf("Rotor: NREL 5MW dissertation blade,  %d aero nodes\n", nr)
@printf("Sim:   two-way coupled, tvec = [0, %.2f]s  (%d time steps)\n", tvec[end], n_steps)

aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=false)

# --- Enable the trace, run the sim, disable ---
# `initialize_sim` also calls solve_BEMT to establish the initial condition,
# so we enable the trace *before* it. That means the tally includes both
# initial-condition solves and every per-time-step solve.
WATT.BEMT_QUAD_TRACE[] = Dict(:q1=>0, :q2=>0, :q3=>0, :q4=>0, :other=>0)

WATT.run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; verbose=false)

tally = WATT.BEMT_QUAD_TRACE[]
WATT.BEMT_QUAD_TRACE[] = nothing   # restore no-overhead default

# --- Report ---
total = sum(values(tally))
println("\nTotal successful solve_BEMT calls: $total")
println("Expected order of magnitude:       ~n_steps * n_aero_nodes = $(n_steps * nr)")
println("(May be smaller — bracket-fail sections and hub/tip cylinders are excluded)\n")

println("Quadrant tally:")
for q in (:q1, :q2, :q3, :q4, :other)
    pct = 100 * tally[q] / max(total, 1)
    @printf("  %-6s  %10d   (%6.3f%%)\n", string(q), tally[q], pct)
end

q1_rate = tally[:q1] / max(total, 1)
@printf("\nq1 hit rate under unsteady coupled loading: %.4f%%\n", 100 * q1_rate)
if q1_rate >= 0.999
    println("→ q1-only GPU kernel remains a safe assumption in the coupled path.")
else
    println("→ q1-only kernel would miss $(round((1-q1_rate)*100, digits=3))% of transient solves.")
    println("  Consider extending the GPU kernel with the observed fallback quadrants.")
end

# ---------------------------------------------------------------------------
# Where did the non-q1 hits happen?
#
# The source-level trace only counts by quadrant; it doesn't know the
# timestep. We recover (timestep, section) by inspecting aerostates.phi:
# each solved section stores its converged inflow angle, and we can
# classify each (i, j) by which quadrant that phi falls into.
#
# What we're looking for:
#   - Which sections?       Root cylinders, mid-span, or tip?
#   - Which timesteps?      Concentrated in a few gusts, or scattered?
#   - How do they cluster into GPU warps?  Warps of 32 threads = 32
#     (section, sim) pairs. In this CPU sim we have only 1 sim, so a
#     "warp" = ceil(n_aero / 32) = 2 warps per timestep. On GPU with
#     n_sims=20 we'd have ceil(25*20 / 32) ≈ 16 warps per step.
# ---------------------------------------------------------------------------

const EPS_PHI = 1e-6
function classify(phi)
    isapprox(phi, 0.0; atol=1e-12)                 && return :zero    # Outputs() / bracket-failed
    if     EPS_PHI  <= phi <= pi/2                   return :q1
    elseif -pi/2    <= phi <= -EPS_PHI               return :q2
    elseif pi/2     <= phi <= pi - EPS_PHI           return :q3
    elseif -pi + EPS_PHI <= phi <= -pi/2             return :q4
    end
    return :other
end

# Locate every non-q1 hit in the sim history.
nonq1 = NamedTuple{(:i, :j, :t, :r, :rR, :phi_deg, :quadrant),
                    Tuple{Int, Int, Float64, Float64, Float64, Float64, Symbol}}[]
for i in 1:size(aerostates.phi, 1), j in 1:size(aerostates.phi, 2)
    q = classify(aerostates.phi[i, j])
    if q !== :q1 && q !== :zero
        push!(nonq1, (i=i, j=j, t=tvec[i], r=blade.r[j], rR=blade.r[j]/Rtip,
                       phi_deg=aerostates.phi[i, j] * 180/pi, quadrant=q))
    end
end

println("\n--- Non-q1 solves recovered from aerostates.phi ---")
println("Found $(length(nonq1)) non-q1 hits (should be close to the tally's non-q1 count $(sum(tally[q] for q in (:q2,:q3,:q4,:other))))")

if !isempty(nonq1)
    # Per-section distribution
    println("\nHits by section (aero-node index and radial position):")
    by_section = Dict{Int, Int}()
    for e in nonq1
        by_section[e.j] = get(by_section, e.j, 0) + 1
    end
    for j in sort(collect(keys(by_section)))
        @printf("  j=%2d  r=%6.2f m  r/R=%.3f  hits=%4d\n",
                j, blade.r[j], blade.r[j]/Rtip, by_section[j])
    end

    # Per-timestep clustering
    println("\nHits by timestep (bins of 10 steps = 0.5s):")
    step_bin(i) = fld(i - 1, 10) + 1
    by_bin = Dict{Int, Int}()
    for e in nonq1
        b = step_bin(e.i)
        by_bin[b] = get(by_bin, b, 0) + 1
    end
    for b in sort(collect(keys(by_bin)))
        t0 = (b - 1) * 0.5
        @printf("  t ∈ [%5.2f, %5.2f) s   hits=%4d\n", t0, t0 + 0.5, by_bin[b])
    end

    # Sample rows for spot-checking
    println("\nFirst up-to-20 non-q1 hits:")
    for e in nonq1[1:min(end, 20)]
        @printf("  step=%3d (t=%5.2fs)  j=%2d (r/R=%.3f)  phi=%+6.1f°  quadrant=%s\n",
                e.i, e.t, e.j, e.rR, e.phi_deg, e.quadrant)
    end
end
