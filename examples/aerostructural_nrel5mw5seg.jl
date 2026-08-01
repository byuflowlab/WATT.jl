#=
Unsteady aerostructural example using the NREL 5MW dissertation blade.

Loads the PreComp-derived geometry and GXBeamCS-derived element compliance
and mass matrices from `data/nrel5mw_5seg.jld2` (produced by
`examples/generate_nrel5mw_5seg_fixture.jl`), builds the dsairfoils from the
OpenFAST `.dat` files in `data/openfast/airfoils/`, drives the rotor with
WATT's existing TurbSim file, and runs a two-way coupled aerostructural
simulation via `initialize_sim` + `run_sim!`.

Plots:
  - Out-of-plane blade-tip deflection vs time
  - Spanwise Fx and Fy at the final time step

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, JLD2
using Plots

const of = OpenFASTTools
const DS = DynamicStallModels

datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")

# --- Load PreComp blade fixture ---
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)

# --- Build dsairfoils from the existing OpenFAST .dat files ---
aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat",  "DU25_A17.dat",  "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]

af_idx    = of.integerfit(raf, afidx, rvec)
afs       = aftypes[af_idx]

dsairfoils = Vector{DS.Airfoil}(undef, nr)
xcp        = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])

    # Overwrite the per-station steady polar with the interpolated table baked
    # into the fixture (matches the dissertation's uo.get_interp_polars usage).
    # Cylinder stations additionally switch to NoModel — their flat polar would
    # send the Beddoes–Leishman state update to Inf in separationpoint.
    if cylinder_mask[i]
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; dsmodel=DS.NoModel(), polar=polars[i])
    else
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
    end
end

# --- Blade, Rotor, Assembly ---
blade = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)

B     = 3
hubHt = 90.0
rotor = WATT.Rotor(B, hubHt, true; tilt=0.0, yaw=0.0)

assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)

# --- Turbulent inflow environment (dissertation values) ---
rho      = 1.225
mu       = 1.837e-5
a        = 335.0
shearexp = 0.2
Vrated   = 11.4
tsr      = 7.55
rotorR   = Rtip * cos(precone)
omega    = Vrated * tsr / rotorR

turbfile = joinpath(ofpath, "TurbSim.dat")
env      = environment(turbfile, rho, mu, a, omega, shearexp)

# --- Simulate ---
tvec = collect(0:0.05:10.0)

aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=true)
WATT.run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; verbose=true)

# --- Plot ---
tip_def = [gx.elements[end].u[3] for gx in gxhistory]

tipplt = plot(tvec, tip_def,
              xlabel="Time (s)", ylabel="Tip deflection (m)",
              title="Out-of-plane tip deflection", label=false, lw=2)

loadplt = plot(blade.r, aerostates.Fx[end, :], label="Fx",
               title="Final-step spanwise loads",
               xlabel="Blade span (m)", ylabel="Force (N/m)", lw=2)
plot!(loadplt, blade.r, aerostates.Fy[end, :], label="Fy", lw=2)

plt = plot(tipplt, loadplt, layout=(2, 1), size=(700, 600))
display(plt)
