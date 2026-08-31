# Steady Aeroelastic Analysis

## Fixed-point iteration
A fixed-point iteration on a given inflow condition. It doesn't use the dynamic stall and doesn't produce a time history.

Each solve:

1. Solves BEMT at every aerodynamic station, using the deflections and
   velocities from the previous sweep.
2. Dimensionalizes the coefficients into distributed loads.
3. Integrates those loads onto the structural elements and calls `GXBeam.steady_state_analysis!`.
4. Interpolates the new deflections back onto the aerodynamic nodes.

Repeat until the loads stop moving. Because both sub-problems are solved to tight tolerances internally, the outer iteration converges quickly — see the numbers below.

## Setting up

The steady solver needs the same `blade`, `rotor`, `env`, and `assembly` objects as the transient one, at the same NREL 5MW rated conditions — 10 m/s, TSR 7.55, 25 aerodynamic stations, 24 structural elements.

```julia
using WATT, GXBeam, DynamicStallModels, OpenFASTTools, StaticArrays

DS = DynamicStallModels
of = OpenFASTTools

ofpath = joinpath(pkgdir(WATT), "data", "openfast")
afpath = joinpath(ofpath, "airfoils")

# Airfoils
aftypes = [of.read_airfoilinput(joinpath(afpath, f)) for f in
           ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
            "DU30_A17.dat", "DU25_A17.dat", "DU21_A17.dat", "NACA64_A17.dat")]

adblade = of.read_adblade("sn5_adblade.dat", ofpath)
af_idx  = Int.(adblade["BlAFID"])

made     = [of.make_dsairfoil(aftypes[i]) for i in af_idx]
airfoils = first.(made)
xcp      = last.(made)

# Blade
edfile   = of.read_edfile("sn5_EDfile.dat", ofpath)

rhub     = edfile["HubRad"]
rvec     = adblade["BlSpn"] .+ rhub
rtip     = rvec[end]
chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"] .* (pi/180)

blade = Blade(rvec, chordvec, twistvec, xcp, airfoils; rhub, rtip)

# Rotor
B       = 3
hubht   = 80.0
turbine = true

rotor = Rotor(B, hubht, turbine; tilt = 0.0, yaw = 0.0)

# Structural assembly
bdfile  = of.read_bdfile("sn5_BDfile.dat", ofpath)
bdblade = of.read_bdblade("sn5_BDblade.dat", ofpath)

assembly = of.make_assembly(edfile, bdfile, bdblade)

# Environment
rho      = 1.225
mu       = 1.464e-5
a        = 343.0
shearexp = 0.0

vinf  = 10.0
tsr   = 7.55
omega = vinf*tsr/rtip

env = environment(rho, mu, a, vinf, omega, shearexp)
```

Note that the steady solver takes wind speed and rotor speed from the environment at `t = 0`, so a time-varying environment is legal — it is simply sampled once.

## Running it

```julia
azimuth0 = 0.0    # blade azimuth, held fixed (rad)
pitch    = 0.0    # blade pitch (rad)

aerostates, gxstates, mesh = initialize_static(blade, assembly)

fixedpoint!(aerostates, gxstates, azimuth0, rotor, env, blade, mesh, pitch;
            iterations = 10)
```

[`initialize_static`](@ref) allocates the buffers and builds the coupling mesh;
[`fixedpoint!`](@ref) runs the sweeps and mutates all three in place.

Two things to know about the arguments:

- **`iterations` is a fixed count, not a tolerance.** The `atol` keyword is
  reserved for a future convergence-tracked exit but is not consulted today, so
  the solver always runs exactly `iterations` sweeps. Pick the count from the
  table below.
- **`azimuth0` matters** even though the analysis is steady, because gravity
  enters as `(-g·cos(azimuth0), -g·sin(azimuth0), 0)`. A blade pointing up and a
  blade pointing sideways reach different equilibria.



Measured on NREL 5MW at rated, tracking the relative change in the spanwise
normal load, `norm(Fx - Fx_prev)/norm(Fx)`:

| Iterations | Thrust (kN) | Torque (kN·m) | Tip flap (m) | Relative change in `Fx` |
|---:|---:|---:|---:|---:|
| 1 | 633.08 | −3417.45 | −4.9396 | — |
| 2 | 616.59 | −3218.70 | −4.7336 | 3.9e−2 |
| 3 | 617.30 | −3227.03 | −4.7446 | 2.2e−3 |
| 5 | 617.26 | −3226.63 | −4.7441 | 1.0e−4 |
| 10 | 617.26 | −3226.63 | −4.7441 | 2.6e−7 |
| 20 | 617.26 | −3226.63 | −4.7441 | 8.4e−14 |

Roughly an order of magnitude per sweep. **Ten iterations is a good default** —
it puts the loads well inside any modeling error, and costs little. Twenty
reaches machine precision if you need a reference solution. 

The torque sign is negative here because it follows the structural axis
convention, not the rotor-power convention; take its magnitude for power.

## Reading the results

`aerostates` is a [`StaticAeroStates`](@ref) — nine spanwise vectors, one entry
per aerodynamic station, with no time dimension:

```julia
aerostates.phi      # inflow angle (rad)
aerostates.alpha    # angle of attack (rad)
aerostates.W        # relative velocity magnitude (m/s)
aerostates.Cx       # normal force coefficient
aerostates.Cy       # tangential force coefficient
aerostates.Cm       # moment coefficient
aerostates.Fx       # normal force per unit span (N/m)
aerostates.Fy       # tangential force per unit span (N/m)
aerostates.Mx       # moment per unit span (N·m/m)
```


The structural side is different from the transient path: `gxstates` is a raw
GXBeam state vector, not an `AssemblyState`. To read deflections, build the
assembly state from the system in the mesh:

```julia
using GXBeam

state = GXBeam.AssemblyState(mesh.system, assembly;
                             prescribed_conditions = mesh.prescribed_conditions)

tip_u     = state.points[end].u        # [-0.311, -0.509, -4.744] m
tip_theta = state.points[end].theta    # Wiener-Milenkovic parameters
```

Convert the rotation parameters to Euler angles with `WATT.WMPtoangle` if you
want physical angles rather than Wiener-Milenkovic parameters.
