#=
Aero-only example: NREL 5MW rotor under turbulent inflow.

Builds the `Blade`, `Rotor`, and `SimpleEnvironment` from scratch by reading
the OpenFAST AeroDyn / ElastoDyn / BeamDyn input files and the TurbSim
time-series file in `data/openfast/`. Runs an aero-only simulation
(`initialize_aero` + `simulate!`) and plots blade-spanwise loads at the
first and last time step.

Adam Cardoza
=#

using WATT, FLOWMath, OpenFASTTools, DynamicStallModels
using StaticArrays
using Plots

const _of = OpenFASTTools
const _DS = DynamicStallModels

ofpath = joinpath(@__DIR__, "..", "data", "openfast")

# --- Read OpenFAST inputs ---------------------------------------------------
adblade = _of.read_adblade("sn5_adblade.dat", ofpath)
edfile  = _of.read_edfile("sn5_EDfile.dat", ofpath)

aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat", "DU25_A17.dat", "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [_of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]

af_idx = Int.(adblade["BlAFID"])
afs    = aftypes[af_idx]

chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"]            # degrees from OpenFAST
rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
n    = length(rvec)

airfoils = Vector{_DS.Airfoil}(undef, n)
xcp = Vector{Float64}(undef, n)
for i = 1:n
    airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
end

# --- Build Blade + Rotor ----------------------------------------------------
blade = WATT.Blade(rvec, chordvec, twistvec .* (pi/180), xcp, airfoils; rhub, rtip)

B       = 3
hubht   = 80.0
turbine = true
rotor   = WATT.Rotor(B, hubht, turbine)

# --- Build turbulent inflow environment from TurbSim file -------------------
rho      = 1.225
mu       = 1.464e-5
a        = 343.0
shearexp = 0.1
RPM      = 11.44
omega    = RPM * 2pi / 60

turbfile = joinpath(ofpath, "TurbSim.dat")
env = environment(turbfile, rho, mu, a, omega, shearexp)

# --- Simulate ---------------------------------------------------------------
tvec  = collect(0:0.05:3.0)
pitch = 0.0

aerostates, mesh = WATT.initialize_aero(blade, tvec)
WATT.simulate!(aerostates, mesh, rotor, blade, env, tvec; pitch)

# --- Plot first/last spanwise loads ----------------------------------------
lb, ub = -1000.0, 7000.0

initialplt = plot(blade.r, aerostates.Fx[1, :],   label="Fx", title="Initial loads", ylim=(lb, ub))
plot!(initialplt, blade.r, aerostates.Fy[1, :],   label="Fy")
finalplt   = plot(blade.r, aerostates.Fx[end, :], label="Fx", title="Final loads",   ylim=(lb, ub))
plot!(finalplt,   blade.r, aerostates.Fy[end, :], label="Fy")
plt = plot(initialplt, finalplt, layout=(2, 1), xlabel="Blade span (m)", ylabel="Force (N/m)")
display(plt)
