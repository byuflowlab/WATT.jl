#=
Regenerate the Phase 3 reference simulation `simpleNREL5_100s.jld2`.

Saves only the *data* of the simulation plus the env *scalars* and a path to
the TurbSim file — env itself is not saved, and is rebuilt at test time via
`WATT.environment(turbfile, rho, mu, a, omega, shearexp)`. This keeps the
JLD2 small (~42 MB instead of ~114 MB) and avoids serializing the 3 large
Akima interpolant tables that the TurbulentInflow callable holds.

Saves the result to a NEW filename (`simpleNREL5_100s_v2.jld2`) so the old
file is not touched until results have been verified.

Steps:
  1. Load the old JLD2. Only the data items (blade, assembly, tvec_r, rotor_r,
     aerostates, gxhistory) are used — `env` deserializes as a
     `JLD2.ReconstructedTypes` shell and is discarded after pulling its
     scalars.
  2. Build a fresh env using WATT's new callable-struct constructors. Only
     `env.U` and `env.RS` are queried by the solvers, so the choice of
     `Vinf`, `Udot`, etc. is metadata and does not affect numerical outputs.
  3. Re-run `initialize_sim` + `run_sim!` with the new env.
  4. Compare against the loaded reference aerostates and gxhistory at the
     same RTOL/ATOL used by the test suite. Print a per-field max abs-diff
     report. If any field exceeds tolerance, stop without writing.
  5. If all fields match, save the new JLD2 (scalars + paths, no env object).

Usage:  julia --project=test test/regenerate_reference.jl

After verifying the v2 file, you can move it over the original:
    mv test/simpleNREL5_100s_v2.jld2 test/simpleNREL5_100s.jld2

Adam Cardoza 2026-05-20
=#

using WATT, JLD2, GXBeam, FLOWMath, DelimitedFiles, StaticArrays, LinearAlgebra

const RTOL = 1e-8
const ATOL = 1e-10

localpath = @__DIR__
cd(localpath)

const OLD_FILE = "simpleNREL5_100s.jld2"
const NEW_FILE = "simpleNREL5_100s_v2.jld2"


# ---------------------------------------------------------------------------
# 1. Load the reference data. The deserialized `env` is a closure-based
#    SimpleEnvironment that JLD2 reconstructs as a stub — discard it.
# ---------------------------------------------------------------------------

println("Loading $OLD_FILE …")
data = jldopen(OLD_FILE, "r") do f
    (
        blade      = f["blade"],
        assembly   = f["assembly"],
        tvec_r     = f["tvec_r"],
        rotor_r    = f["rotor_r"],
        aerostates = f["aerostates"],
        gxhistory  = f["gxhistory"],
        env_       = f["env"],     # closure-based; only fluid scalars are usable
    )
end

blade           = data.blade
assembly        = data.assembly
tvec_r          = data.tvec_r
rotor_r         = data.rotor_r
aerostates_ref  = data.aerostates
gxhistory_ref   = data.gxhistory
env_            = data.env_

println("  blade:           ", typeof(blade))
println("  assembly:        ", typeof(assembly))
println("  tvec_r:          length=", length(tvec_r), ", t∈[", tvec_r[1], ", ", tvec_r[end], "]")
println("  rotor_r:         B=", rotor_r.B, ", hubht=", rotor_r.hubht)
println("  aerostates_ref:  size(phi)=", size(aerostates_ref.phi))
println("  gxhistory_ref:   length=", length(gxhistory_ref))


# ---------------------------------------------------------------------------
# 2. Build a fresh env using the new callable structs. Match the exact same
#    Akima fits the original test built (TurbSim file at the standard path).
#    Hard-coded scalars match what the original `simple_NREL5MW.jl` test used.
# ---------------------------------------------------------------------------

# Fluid scalars must match the saved env exactly. Pull them from env_ rather
# than hard-coding — these differed in the original fixture (mu=1.8375e-5,
# a=335.0, shearexp=0.1) from what one would naively expect at sea level.
rho      = env_.rho
mu       = env_.mu
a        = env_.a
shearexp = env_.shearexp
RPM = 11.44
omega = RPM * 2pi / 60
println("\nScalars from loaded env: rho=$rho, mu=$mu, a=$a, shearexp=$shearexp, omega=$omega rad/s")

const TURBFILE = "../data/openfast/TurbSim.dat"
env = environment(TURBFILE, rho, mu, a, omega, shearexp)
println("\nBuilt new env via WATT.environment(turbfile, ...): ", typeof(env))


# ---------------------------------------------------------------------------
# 3. Re-run the simulation against the new env.
# ---------------------------------------------------------------------------

println("\nRunning initialize_sim + run_sim! …")
aerostates_new, gxhistory_new, mesh_new = WATT.initialize_sim(blade, assembly, tvec_r; verbose=false)
WATT.run_sim!(rotor_r, blade, mesh_new, env, tvec_r, aerostates_new, gxhistory_new; verbose=false)
println("  done.")


# ---------------------------------------------------------------------------
# 4. Compare against the reference field-by-field.
# ---------------------------------------------------------------------------

println("\nComparing aerostates fields against reference (RTOL=$RTOL, ATOL=$ATOL):")

aerostate_fields = (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
all_ok = Ref(true)
for f in aerostate_fields
    a_new = getfield(aerostates_new, f)
    a_ref = getfield(aerostates_ref, f)
    diff = abs.(a_new .- a_ref)
    maxabs = maximum(diff)
    maxref = max(maximum(abs.(a_ref)), eps())
    ok = isapprox(a_new, a_ref; rtol=RTOL, atol=ATOL)
    all_ok[] &= ok
    status = ok ? "OK  " : "FAIL"
    println("  $status  $(rpad(string(f), 8))  max|Δ| = $(maxabs)  (max|ref| = $(maxref))")
end

println("\nComparing gxhistory at i ∈ (1, mid, end):")
nt = length(gxhistory_new)
for i in (1, max(1, div(nt, 2)), nt)
    a, b = gxhistory_new[i], gxhistory_ref[i]
    point_u_diff     = maximum(maximum(abs.(pa.u     .- pb.u))     for (pa, pb) in zip(a.points, b.points))
    point_theta_diff = maximum(maximum(abs.(pa.theta .- pb.theta)) for (pa, pb) in zip(a.points, b.points))
    elem_u_diff      = maximum(maximum(abs.(ea.u     .- eb.u))     for (ea, eb) in zip(a.elements, b.elements))
    elem_Fi_diff     = maximum(maximum(abs.(ea.Fi    .- eb.Fi))    for (ea, eb) in zip(a.elements, b.elements))
    println("  i=$i  max|Δ point.u|=$point_u_diff  max|Δ point.θ|=$point_theta_diff  max|Δ elem.u|=$elem_u_diff  max|Δ elem.Fi|=$elem_Fi_diff")
    # Strict per-field tolerance
    for (pa, pb) in zip(a.points, b.points)
        all_ok[] &= isapprox(pa.u,     pb.u;     rtol=RTOL, atol=ATOL)
        all_ok[] &= isapprox(pa.theta, pb.theta; rtol=RTOL, atol=ATOL)
    end
    for (ea, eb) in zip(a.elements, b.elements)
        all_ok[] &= isapprox(ea.u,     eb.u;     rtol=RTOL, atol=ATOL)
        all_ok[] &= isapprox(ea.theta, eb.theta; rtol=RTOL, atol=ATOL)
        all_ok[] &= isapprox(ea.Fi,    eb.Fi;    rtol=RTOL, atol=ATOL)
        all_ok[] &= isapprox(ea.Mi,    eb.Mi;    rtol=RTOL, atol=ATOL)
    end
end


# ---------------------------------------------------------------------------
# 5. Save (if results match).
# ---------------------------------------------------------------------------

if !all_ok[]
    println("\n✗ Results do NOT match the reference within RTOL=$RTOL, ATOL=$ATOL.")
    println("  NOT writing $NEW_FILE. Investigate the deltas above before regenerating.")
    exit(1)
end

println("\n✓ Results match the reference within RTOL=$RTOL, ATOL=$ATOL.")
println("Writing $NEW_FILE …")
#=
Save only the simulation data plus env scalars and the TurbSim file path.
The env is rebuilt at test time via `WATT.environment(turbfile, …)`. This
keeps the JLD2 small (the 3 Akima interpolants the env holds for U/V/W are
~70 MB and the TurbSim source is already in the repo).
Also omit `mesh` — simple_NREL5MW.jl loads but never uses it, and the
SimMesh struct holds the default anonymous pfunc closure.
=#
jldsave(NEW_FILE;
    blade      = blade,
    assembly   = assembly,
    tvec_r     = tvec_r,
    rotor_r    = rotor_r,
    aerostates = aerostates_new,
    gxhistory  = gxhistory_new,
    # Env scalars + TurbSim file path (relative to test/) — test rebuilds env.
    rho        = rho,
    mu         = mu,
    a          = a,
    shearexp   = shearexp,
    RPM        = RPM,
    turbfile   = TURBFILE,
)
println("  done. To adopt: mv test/$NEW_FILE test/$OLD_FILE")
