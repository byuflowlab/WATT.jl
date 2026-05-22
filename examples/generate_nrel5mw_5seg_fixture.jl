#=
One-time fixture generator for the dissertation NREL 5MW aerostructural example.

This script bakes the baseline NREL 5MW PreComp-derived blade geometry plus
the GXBeamCS classical-laminate-theory compliance and mass matrices into
`data/nrel5mw_5seg.jld2`. The companion example
(`examples/aerostructural_dissertation.jl`) then loads that fixture and runs
an aerostructural simulation without needing UnsteadyOpt, GXBeamCS, or the
dissertation `data/` directory at runtime.

Run once, then commit `data/nrel5mw_5seg.jld2` alongside the example.

Requires:
  - UnsteadyOpt and GXBeamCS available in the active environment
  - The Cardoza2026_Efficient_aeroelastic_wind_gradients repo accessible at
    `dissertation_root` below

Adam Cardoza
=#

using UnsteadyOpt, GXBeamCS, CCBlade, JLD2

const uo = UnsteadyOpt

dissertation_root = joinpath(homedir(), "repos", "Cardoza2026_Efficient_aeroelastic_wind_gradients")
yamlpath           = joinpath(dissertation_root, "data", "5MW_PreComp_5seg")
precomppath        = joinpath(dissertation_root, "data", "5MW_PreComp_5seg")
airfoil_interp_path = joinpath(dissertation_root, "data", "airfoils_interpolated")

# --- Blade geometry from PreComp ---
distro, materials, stations = uo.get_precomp_descriptions(precomppath, yamlpath)

rvec     = distro.rvec
cvec     = distro.cvec
twistvec = distro.twistvec
le_loc   = distro.le_loc
nr       = length(rvec)

# Constants matching the dissertation NREL 5MW setup
Rhub    = 1.5
Rtip    = 63.0
precone = 2.5 * pi / 180

# Airfoil-station mapping (kept as data so the example doesn't need UnsteadyOpt).
# Used by the example to pick the *nearest* reference AirfoilInput for DS state
# coefficients; the steady polar at each station is overwritten with the
# interpolated table assembled below.
raf   = [2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25,
         36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333]
afidx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

# --- Per-station interpolated steady polars (matches uo.get_interp_polars usage in optimize.jl) ---
# Each entry is a [alpha cl cd cm] table. Cylinder stations get a 3-row constant
# table (zero cl, constant cd) and are flagged so the example can swap their DS
# model to NoModel() — the BL state update otherwise divides by dCn/dα ≈ 0.
airfoils_ccb  = uo.get_interp_polars(rvec, airfoil_interp_path)
polars        = Vector{Matrix{Float64}}(undef, nr)
cylinder_mask = falses(nr)
for i = 1:nr
    if isa(airfoils_ccb[i], CCBlade.Cylinder)
        cd_const     = airfoils_ccb[i].cd
        polars[i]    = [-pi  0.0  cd_const  0.0;
                         0.0  0.0  cd_const  0.0;
                         pi   0.0  cd_const  0.0]
        cylinder_mask[i] = true
    else
        alpha     = airfoils_ccb[i].alpha
        polars[i] = hcat(alpha, airfoils_ccb[i].cl, airfoils_ccb[i].cd, zeros(length(alpha)))
    end
end

# --- Baseline CLT cross-section analysis ---
num_segments = uo.num_segments5
fvec_base    = ones(5 * nr)                       # baseline (un-scaled) thickness factors
materials_f  = GXBeamCS.Material{Float64}.(materials)
segs_webs    = uo.scale_segments(stations, fvec_base, materials_f, num_segments)

clt_list = [uo.get_clt(
    stations[i].xaf,
    stations[i].yaf,
    cvec[i],
    twistvec[i],
    le_loc[i],
    stations[i].xbreak,
    stations[i].webloc,
    segs_webs[i].segments,
    segs_webs[i].webs,
) for i in 1:nr]

shear_center    = true
compliance_list = [GXBeamCS.compliance_matrix(clt, shear_center)[1] for clt in clt_list]
mass_list       = [GXBeamCS.mass_matrix_clt(clt_list[i]; reference=[cvec[i]*le_loc[i], 0.0, twistvec[i]])[1] for i in 1:nr]

# --- GXBeam discretization (matches optimize.jl) ---
pts    = zeros(nr + 1)
pts[1] = Rhub
for i = 1:nr-1
    pts[i+1] = (rvec[i] + rvec[i+1]) / 2
end
pts[end] = Rtip

points = [[pts[i], 0.0, 0.0] for i in eachindex(pts)]
xp     = [[rvec[i], 0.0, 0.0] for i in 1:nr]
nelem  = length(points) - 1
start  = collect(1:nelem)
stop   = collect(2:nelem+1)

# --- Save ---
@save "nrel5mw_5seg.jld2" rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask

println("Wrote nrel5mw_5seg.jld2")
println("  nr     = $nr")
println("  nelem  = $nelem")
println("  rvec   ∈ [$(rvec[1]), $(rvec[end])]")
