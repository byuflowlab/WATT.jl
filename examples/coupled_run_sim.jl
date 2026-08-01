#=
Original coupled aerostructural solution (reference example).

The canonical two-way coupled transient: build the NREL 5MW dissertation blade from the
`data/nrel5mw_5seg.jld2` fixture, drive it with a turbulent inflow, and march the coupled
aero-structural solution from rest with `initialize_sim` + `run_sim!`.

This is the companion/reference to `examples/step_solution_ad.jl`, which restarts from a snapshot
of *this* solution and differentiates a short window with ForwardDiff.

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, JLD2
using Plots

const of = OpenFASTTools
const DS = DynamicStallModels

datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")

# --- Load the PreComp/GXBeamCS blade fixture (geometry + per-element compliance & mass). ---
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)

# --- Dynamic-stall airfoils from the OpenFAST .dat files. ---
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

# --- Blade, rotor, structural assembly. ---
blade    = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)
rotor    = WATT.Rotor(3, 90.0, true; tilt=0.0, yaw=0.0)
assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)

# --- Turbulent inflow environment (dissertation values). ---
rho, mu, a, shearexp = 1.225, 1.837e-5, 335.0, 0.2
Vrated, tsr = 11.4, 7.55
omega    = Vrated * tsr / (Rtip * cos(precone))
turbfile = joinpath(ofpath, "TurbSim.dat")
env      = environment(turbfile, rho, mu, a, omega, shearexp)

# --- Coupled transient from rest. ---
tvec = collect(0.0:0.05:10.0)

aerostates, gxhistory, mesh = initialize_sim(blade, assembly, tvec; verbose=true)
run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; verbose=true)

# --- Report / plot. ---
tip_def = [gx.elements[end].u[3] for gx in gxhistory]
println("Final out-of-plane tip deflection: ", round(tip_def[end]; digits=4), " m")

tipplt = plot(tvec, tip_def; xlabel="Time (s)", ylabel="Tip deflection (m)",
              title="Out-of-plane tip deflection", label=false, lw=2)
loadplt = plot(blade.r, aerostates.Fx[end, :]; label="Fx", lw=2,
               title="Final-step spanwise loads", xlabel="Blade span (m)", ylabel="Force (N/m)")
plot!(loadplt, blade.r, aerostates.Fy[end, :]; label="Fy", lw=2)
display(plot(tipplt, loadplt; layout=(2, 1), size=(700, 600)))
