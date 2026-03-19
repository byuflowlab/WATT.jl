# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

WATT (Wind Aeroelastic Turbine Toolkit) is a Julia package for nonlinear unsteady aeroelastic modeling of wind turbine blades, designed with AD (automatic differentiation) compatibility as a first-class concern. It couples three physics solvers:

- **CCBlade.jl** — Blade Element Momentum Theory (BEMT) aerodynamics
- **DynamicStallModels.jl** — Unsteady aerodynamic (dynamic stall) corrections
- **GXBeam.jl** — Geometrically Exact Beam Theory (GEBT) structural dynamics

AD is supported via ForwardDiff, ReverseDiff, and ImplicitAD throughout the codebase.

## Commands

```julia
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run tests from REPL
using Pkg; Pkg.test("WATT")

# Run a specific test file from the test/ directory
julia --project -e 'include("test/test_bem.jl")'

# Build documentation
julia --project=docs/ docs/make.jl
```

It is normal for warnings to appear during tests — they exercise boundary conditions intentionally.

## Architecture

```
src/
├── WATT.jl               # Package entry, includes all submodules
├── types.jl              # Rotor and Blade data structures
├── environments.jl       # Wind field: SimpleEnvironment with time-varying velocity/RPM
├── bem.jl                # BEMT aerodynamic solver wrapping CCBlade
├── dynamicstallmodels.jl # Wrapper for DynamicStallModels.jl
├── gxbeam.jl             # Structural solver wrapping GXBeam
├── solvers.jl            # Custom ODE solvers: RK4, BDF1
├── mesh.jl               # Aero/structural node interpolation utilities
├── utils.jl              # Rotation transforms, AD-compatible helpers (dualcopy, etc.)
├── aerostructural.jl     # Main transient coupled solver (initialize, simulate, take_step!)
├── aero_only.jl          # Aerodynamics-only transient analysis
├── static.jl             # Steady-state fixed-point aerostructural solver
└── beddoesleishman*.jl   # Legacy dynamic stall implementations (3 variants)
```

### Coupling Strategy

The package uses three simulation modes:

1. **`aerostructural.jl`** — Full two-way coupled transient simulation. Entry points: `initialize()`, `initial_condition()`, `take_step!()`, `simulate()`/`simulate!()`.
2. **`aero_only.jl`** — Aerodynamics-only transient analysis (no structural deformation).
3. **`static.jl`** — Steady-state fixed-point iteration between CCBlade and GXBeam.

### AD Compatibility

Custom ODE solvers (`RK4`, `BDF1` in `solvers.jl`) exist to maintain full AD compatibility. The `dualcopy()` utility in `utils.jl` is used to copy arrays while preserving ForwardDiff/ReverseDiff dual number types. ImplicitAD is used at implicit solve boundaries.

### Key Data Structures (`types.jl`)

- `Rotor` — Rotor-level properties: number of blades, hub height, tilt/yaw angles, correction model flags.
- `Blade` — Per-blade aerodynamic and structural properties: radial stations, twist, chord, airfoil tables, stiffness matrices.
