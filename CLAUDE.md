# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

- Find a project overview in README.md
- This project is written in Julia: Style, efficiency, proficiency see .claude/context/JULIA_GUIDE.md; rules see .claude/rules/julia
- Use julia's testing environment: `julia --project -e 'using Pkg; pkg.test()`
- Use docs.jl to compile documentation: `julia --project=docs/ docs/make.jl`


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
├── aerostructural.jl     # Main transient coupled solver (initialize_sim, run_sim!)
├── aero_only.jl          # Aerodynamics-only transient analysis
└── static.jl             # Steady-state fixed-point aerostructural solver
```

### Coupling Strategy

The package uses three simulation modes:

1. **`aerostructural.jl`** — Full two-way coupled transient simulation. Entry points: `initialize_sim()` and `run_sim!()`. These are the only active Gen 2 API functions; a prior Gen 1 API (`initialize`, `initial_condition!`, `take_step!`, `simulate`, `simulate!`) has been removed.
2. **`aero_only.jl`** — Aerodynamics-only transient analysis (no structural deformation).
3. **`static.jl`** — Steady-state fixed-point iteration between CCBlade and GXBeam.

### AD Compatibility

Custom ODE solvers (`RK4`, `BDF1` in `solvers.jl`) exist to maintain full AD compatibility. The `dualcopy()` utility in `utils.jl` is used to copy arrays while preserving ForwardDiff/ReverseDiff dual number types. ImplicitAD is used at implicit solve boundaries.

### Key Data Structures (`types.jl`)

- `Rotor` — Rotor-level properties: number of blades, hub height, tilt/yaw angles, correction model flags.
- `Blade` — Per-blade aerodynamic and structural properties: radial stations, twist, chord, airfoil tables, stiffness matrices.
