# API Reference

Every symbol `using WATT` brings into scope, grouped by what it is for. Anything
not listed here is internal — you can still reach it as `WATT.name`, but it may
change without notice.

If you are looking for a place to start rather than a lookup table, read
[Getting Started](gettingstarted.md) first.

```@meta
CurrentModule = WATT
```

## Choosing a solver

WATT offers four ways to run a rotor, in increasing order of what they model:

| Entry point | Aerodynamics | Structures | Use it when |
|---|---|---|---|
| [`initialize_static`](@ref) + [`fixedpoint!`](@ref) | BEMT, steady | GXBeam, steady | You want a trim state or a warm start |
| [`initialize_aero`](@ref) + [`simulate!`](@ref) | BEMT + dynamic stall | rigid | The blade is stiff, or you are isolating aero |
| [`initialize_sim`](@ref) + [`run_sim!`](@ref) | BEMT + dynamic stall | GXBeam, transient | The general case — full two-way coupling |
| [`initialize_sim_surrogate`](@ref) + [`run_sim_surrogate!`](@ref) | BEMT + dynamic stall | learned surrogate | The structural solve is your bottleneck |

For gradients over a short window from a frozen start, see
[Windowed sensitivities](#Windowed-sensitivities) below.

## Types

### Rotor and blade geometry

```@docs
Rotor
Blade
```

### Simulation state

```@docs
AeroStates
StaticAeroStates
```

### Meshes

The mesh types hold the coupling buffers reused at every time step — the
interpolation weights between aerodynamic and structural nodes, the dynamic
stall state, and the GXBeam system. They are built for you by the
`initialize_*` functions; you rarely construct one directly.

```@docs
AbstractSimMesh
SimMesh
AeroMesh
StaticMesh
SurrogateMesh
```

## Solvers

Custom integrators, kept in-package rather than delegating to
DifferentialEquations so that every step stays transparent to ForwardDiff and
ReverseDiff.

```@docs
RK4
BDF1
```

## Environment

Listed per method so each can be linked individually: steady inflow, a whole
TurbSim file, and a TurbSim file from a given time offset.

```@docs
environment(::Number, ::Number, ::Number, ::Number, ::Number, ::Number)
environment(::String, ::Number, ::Number, ::Number, ::Number, ::Number)
environment(::String, ::Number, ::Number, ::Number, ::Number, ::Number, ::Number)
SimpleEnvironment
```

## Steady analysis

```@docs
initialize_static
fixedpoint!
```

## Aerodynamics-only transient analysis

```@docs
initialize_aero
simulate!
```

## Coupled aero-structural transient analysis

```@docs
initialize_sim
run_sim!
run_sim
```

### Windowed sensitivities

These expose a single coupled time step as a callable primitive, so a gradient
can be taken over a short window from a frozen initial state instead of over an
entire march.

```@docs
step_solution!
initialize_from_state
run_from_state!
```

## Structural surrogate

Replaces the GXBeam structural solve with a learned latent model: encode the
initial structural state, march a linear latent step, decode back to
deflections. Implement [`AbstractStructuralSurrogate`](@ref) to plug in your own.

```@docs
AbstractStructuralSurrogate
SurrogatePointState
SurrogateAssemblyState
encode_initial
step_latent
decode
decode!
initialize_sim_surrogate
run_sim_surrogate!
run_sim_surrogate
```

## GPU backend

Device-resident batched types and kernels live on their own page — see
[GPU Backend](gpu.md).

## Post-processing

```@docs
rotorloads
```

## Plotting

`RecipesBase` recipes, so `Plots` is not a dependency of WATT. Load `Plots`
yourself and pass one of these wrappers to `plot`.

```@docs
BladePoints
AssemblyPlot
```

## Internals

Not exported, and not covered by any stability guarantee — reach them as
`WATT.name`. They are documented here because the public docstrings above refer
to them, and because they are where the aero-structural coupling actually
happens. See [For Developers](developers.md) for how they fit together.

```@docs
InterpolationPoint
find_inittype
dualcopy
WMPtoangle
update_mesh!
get_aero_velocities
get_aerostructural_velocities
update_ds_inputs!
extract_ds_loads!
dimensionalize!
build_surrogate_loads
wmp_to_angle_dev
```

## Index

```@index
```
