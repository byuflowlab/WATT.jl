#=
Step 0 of the GPU-BEMT plan: empirically verify the q1-bracket assumption.

For a GPU-batched BEMT with a fixed-iteration Brent root-find, we want to
search a single bracket per section rather than the four-quadrant fallback
in the CPU solver. The natural candidate is q1 = (eps, pi/2), which is where
the root lives during normal wind-turbine operation.

This script builds a representative NREL 5MW rotor from the OpenFAST inputs
in data/openfast, then sweeps (wind speed, TSR, pitch) across an operating
envelope. For every (section, condition) it asks: which of the four quadrants
(q1..q4) is the FIRST that _CC.firstbracket accepts? It reports:

  - overall hit rate for each quadrant,
  - hit rate broken down by radial station,
  - conditions where q1 fails, so we can see if failures are extreme edge
    cases or scattered across normal operation.

Run with:
    julia --project examples/gpu_bemt_bracket_check.jl

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels
using StaticArrays, StructArrays
using Printf

const _of = OpenFASTTools
const _DS = DynamicStallModels
const _CC = WATT.CCBlade

# --- 1. Build the NREL 5MW rotor (same as examples/aero_only_loads.jl) ------

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
rotor = WATT.Rotor(3, 80.0, true)

# --- 2. Envelope sweep --------------------------------------------------------

rho = 1.225
mu  = 1.464e-5
asound = 343.0

# Operating envelope: covers cut-in through cut-out plus stalled/high-thrust edges.
wind_speeds = collect(3.0:1.0:25.0)          # m/s
tsrs        = collect(3.0:0.5:12.0)          # tip speed ratio (Omega*Rtip / U)
pitches     = collect(-2.0:2.0:25.0) .* (pi/180)  # rad

# Reproduce the CPU solver's quadrant-order selection (from src/bemt.jl).
# Returns a tuple of (label, (phimin, phimax)) pairs in the order the CPU
# solver would try them.
const EPS_PHI = 1e-6
const Q1 = (EPS_PHI, pi/2)
const Q2 = (-pi/2, -EPS_PHI)
const Q3 = (pi/2, pi - EPS_PHI)
const Q4 = (-pi + EPS_PHI, -pi/2)

function quadrant_order(Vx, Vy, theta)
    Vx0 = isapprox(Vx, 0.0, atol=1e-6)
    Vy0 = isapprox(Vy, 0.0, atol=1e-6)
    if Vx0 && Vy0
        return nothing, false
    elseif Vx0
        startfrom90 = false
        if Vy > 0 && theta > 0
            return ((:q1, Q1), (:q2, Q2)), startfrom90
        elseif Vy > 0 && theta < 0
            return ((:q2, Q2), (:q1, Q1)), startfrom90
        elseif Vy < 0 && theta > 0
            return ((:q3, Q3), (:q4, Q4)), startfrom90
        else
            return ((:q4, Q4), (:q3, Q3)), startfrom90
        end
    elseif Vy0
        startfrom90 = true
        if Vx > 0 && abs(theta) < pi/2
            return ((:q1, Q1), (:q3, Q3)), startfrom90
        elseif Vx < 0 && abs(theta) < pi/2
            return ((:q2, Q2), (:q4, Q4)), startfrom90
        elseif Vx > 0 && abs(theta) > pi/2
            return ((:q3, Q3), (:q1, Q1)), startfrom90
        else
            return ((:q4, Q4), (:q2, Q2)), startfrom90
        end
    else
        startfrom90 = false
        if Vx > 0 && Vy > 0
            return ((:q1, Q1), (:q2, Q2), (:q3, Q3), (:q4, Q4)), startfrom90
        elseif Vx < 0 && Vy > 0
            return ((:q2, Q2), (:q1, Q1), (:q4, Q4), (:q3, Q3)), startfrom90
        elseif Vx > 0 && Vy < 0
            return ((:q3, Q3), (:q4, Q4), (:q1, Q1), (:q2, Q2)), startfrom90
        else
            return ((:q4, Q4), (:q3, Q3), (:q2, Q2), (:q1, Q1)), startfrom90
        end
    end
end

function backward_search(phimin, phimax, startfrom90)
    if !startfrom90
        return phimin == -pi/2 || phimax == -pi/2
    else
        return phimax == pi/2
    end
end

# For one (section, Vx, Vy, pitch) tuple, return the label of the first
# quadrant that firstbracket accepts, or :none if all four fail.
function which_quadrant(rotor, blade, section_idx, Vx, Vy, pitch; npts=10)
    twist = blade.twist[section_idx]
    theta = twist + pitch
    airfoil = blade.airfoils[section_idx]
    r = blade.r[section_idx]
    chord = blade.c[section_idx]

    xv = zeros(11)
    xv[1] = r
    xv[2] = chord
    xv[3] = twist
    xv[4] = blade.rhub
    xv[5] = blade.rtip
    xv[6] = Vx
    xv[7] = Vy
    xv[8] = rho
    xv[9] = pitch
    xv[10] = mu
    xv[11] = asound

    pv = (airfoil, rotor.B, rotor.turbine, rotor.re, rotor.mach, rotor.rotation, rotor.tip)

    residual(phi) = _CC.residual_and_outputs(phi, xv, pv)[1]

    order, startfrom90 = quadrant_order(Vx, Vy, theta)
    order === nothing && return :none

    for (label, (phimin, phimax)) in order
        backwards = backward_search(phimin, phimax, startfrom90)
        success, _, _ = _CC.firstbracket(residual, phimin, phimax, npts, backwards)
        success && return label
    end
    return :none
end

# --- 3. Sweep and tally -------------------------------------------------------

# Tally: quadrant -> count, plus per-section quadrant counts, plus a small
# list of (Vx, Vy, pitch, section) failures where q1 didn't win.
tally = Dict{Symbol,Int}(:q1=>0, :q2=>0, :q3=>0, :q4=>0, :none=>0)
per_section = zeros(Int, n, 5)   # columns: q1, q2, q3, q4, none
q_idx = Dict(:q1=>1, :q2=>2, :q3=>3, :q4=>4, :none=>5)
non_q1_examples = Vector{NamedTuple}()

function sweep!(tally, per_section, non_q1_examples, rotor, blade, wind_speeds, tsrs, pitches, rtip, n, q_idx)
    total = 0
    for U in wind_speeds, tsr in tsrs, pitch in pitches
        Omega = tsr * U / rtip
        for j in 1:n
            Vx = U
            Vy = Omega * blade.r[j]
            q = which_quadrant(rotor, blade, j, Vx, Vy, pitch)
            tally[q] += 1
            per_section[j, q_idx[q]] += 1
            total += 1
            if q !== :q1 && length(non_q1_examples) < 30
                push!(non_q1_examples, (U=U, tsr=tsr, pitch_deg=pitch*180/pi,
                                        section=j, r=blade.r[j], rR=blade.rR[j],
                                        Vx=Vx, Vy=Vy, quadrant=q))
            end
        end
    end
    return total
end

total = sweep!(tally, per_section, non_q1_examples, rotor, blade, wind_speeds, tsrs, pitches, rtip, n, q_idx)

# --- 4. Report ----------------------------------------------------------------

println("\n=== GPU-BEMT bracket check: q1 assumption verification ===\n")
@printf("Rotor: NREL 5MW-style, %d aerodynamic sections\n", n)
@printf("Envelope: %d wind speeds x %d TSRs x %d pitches x %d sections = %d cases\n\n",
        length(wind_speeds), length(tsrs), length(pitches), n, total)

println("Overall bracket-hit tally:")
for q in (:q1, :q2, :q3, :q4, :none)
    pct = 100 * tally[q] / total
    @printf("  %-5s  %10d   (%6.3f%%)\n", string(q), tally[q], pct)
end

q1_rate = tally[:q1] / total
@printf("\nq1 hit rate: %.4f%%\n", 100 * q1_rate)

println("\nPer-section q1 hit rate (radius, r/R, q1%, non-q1%):")
for j in 1:n
    q1_pct = 100 * per_section[j, 1] / sum(per_section[j, :])
    nonq1_pct = 100 - q1_pct
    @printf("  j=%2d  r=%6.2f m  r/R=%.3f    q1: %6.2f%%   non-q1: %6.2f%%\n",
            j, blade.r[j], blade.rR[j], q1_pct, nonq1_pct)
end

if !isempty(non_q1_examples)
    println("\nFirst up-to-30 non-q1 hits (for spot-checking whether they're edge cases):")
    for ex in non_q1_examples
        @printf("  U=%4.1f  TSR=%4.1f  pitch=%5.1f deg  j=%2d (r/R=%.2f)  Vx=%6.2f Vy=%7.2f  -> %s\n",
                ex.U, ex.tsr, ex.pitch_deg, ex.section, ex.rR, ex.Vx, ex.Vy, ex.quadrant)
    end
end

# --- 5. Decision rule --------------------------------------------------------

println("\n=== Decision (from plan Step 0) ===")
if q1_rate >= 0.999
    println("q1 hit rate >= 99.9%: proceed with q1-only kernel.")
    println("Non-q1 sections will be flagged success=false; caller can fall back to CPU.")
else
    println("q1 hit rate < 99.9%: v1 kernel should also try q3 as a fallback bracket.")
end
