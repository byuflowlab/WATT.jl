# For Developers

How the code is organized, how the models are coupled, and how derivatives get
through it all. Read [Getting Started](gettingstarted.md) first — this page
assumes you know what the user-facing calls look like.

```@meta
CurrentModule = WATT
```

## Code Outline

```
src/
├── WATT.jl                        # Package entry; all exports live here
├── types.jl                       # Rotor, Blade, AeroStates, StaticAeroStates,
│                                  #   SurrogatePointState, SurrogateAssemblyState
├── environments.jl                # Inflow: SimpleEnvironment + callable structs
├── mesh.jl                        # Coordinate transforms, aero↔structural coupling,
│                                  #   the mesh types (AeroMesh/SimMesh/StaticMesh)
├── bemt.jl                        # BEMT solver wrapping CCBlade; ImplicitAD boundary
├── dynamicstallmodels.jl          # Wrapper for DynamicStallModels.jl
├── gxbeam.jl                      # Structural solver wrapping GXBeam; plot recipes
├── solvers.jl                     # AD-transparent ODE solvers: RK4, BDF1
├── utils.jl                       # Rotations, dualcopy, AD-compatible helpers
│
│   # Simulation modes
├── aero_only.jl                   # Aerodynamics only, rigid blade
├── static.jl                      # Steady fixed-point aero-structural
├── aerostructural.jl              # Coupled transient + windowed-AD primitives
├── aerostructural_surrogate.jl    # Coupled transient, learned structural model
│
│   # GPU backends (forward pass only)
├── bemt_gpu.jl                    # Batched BEMT kernel
├── dsmodel_gpu.jl                 # Batched dynamic stall kernel
└── aerostructural_surrogate_gpu.jl # Device-resident coupled march
```

Four entry points, in increasing order of what they model — see the table in the
[API Reference](apireference.md#Choosing-a-solver).

### Where the mesh types fit

The `mesh` argument threaded through every solver is not a geometric mesh. It is
the bundle of coupling state and scratch buffers reused at every time step:
interpolation handles between aerodynamic and structural nodes, the dynamic
stall state and parameters, the GXBeam system, and the deflection and velocity
arrays the aero side reads. There is one type per mode, each carrying only what
that mode needs:

| Type | Mode | Notably carries | Notably omits |
|---|---|---|---|
| [`AeroMesh`](@ref) | aero-only | dynamic stall state | everything GXBeam |
| [`SimMesh`](@ref) | coupled transient | DS state + GXBeam system | — |
| [`StaticMesh`](@ref) | steady | GXBeam system | DS state (`y_ds`, `p_ds`, `xds_idxs`) |
| [`SurrogateMesh`](@ref) | surrogate | DS state + latent surrogate state | GXBeam system |

`StaticMesh` deliberately omits the dynamic stall fields — the steady solver has
no stall model. [`update_mesh!`](@ref) is shared across modes and does not touch
them, which is what makes the omission safe. Keep it that way.

## Reference frames

Four frames, related by rotations built from the rotor and blade geometry:

| Frame | Symbol | Definition |
|---|---|---|
| Global | **G** | Ground-fixed. Where the wind field is defined. |
| Hub-rotating | **HR** | Rotates with the rotor; azimuth, tilt, and yaw removed. |
| Blade-coned | **BC** | HR with precone applied — aligned with the blade reference line. |
| Local element | **L** | The airfoil section frame. Where BEMT and dynamic stall work. |

The transforms are in [mesh.jl](https://github.com/byuflowlab/WATT.jl/blob/master/src/mesh.jl):

```
                 precone
      BC ────────────────────► HR                     WATT.transform_BC_HR
      BC ──── precone → azimuth → tilt → yaw ───► G   WATT.transform_BC_G
      HR ──── sweep → curve → precone ──────────► L   WATT.transform_HR_L
      G  ──── yaw → tilt → azimuth → sweep
                    → curve → precone ───────────► L  WATT.transform_G_L
```

`transform_G_L` is the composite: it is what takes a wind velocity in the ground
frame down to the airfoil section where BEMT can use it. Each step is one of the
`rotate_x`/`rotate_y`/`rotate_z` helpers in `utils.jl`, whose `T` keyword selects
the transpose (inverse) rotation.

Two conventions that trip people up:

- **Rotations are stored as Wiener-Milenkovic parameters**, not Euler angles,
  because that is what GXBeam uses. Convert with `WATT.WMPtoangle`.
- **The moment coefficient is positive about the negative aerodynamic Z axis**,
  so it is negated when mapped onto the structural X axis. See
  [`dimensionalize!`](@ref).

## Current coupling (loose, two-way partitioned)

As in OpenFAST, this is a two-way partitioned loose coupling of explicit
aerodynamic and structural models.

The aerodynamic side is blade element momentum theory plus a Beddoes-Leishman
dynamic stall model. BEMT is well documented in most aerodynamics texts; the
variant here is Ning's, which reduces the problem to a one-dimensional residual
with guaranteed convergence. Environmental flow conditions and static airfoil
data give the inflow angle by solving that residual, as CCBlade does, and the
inflow angle gives the angle of attack.

Because BEMT assumes steady inflow, the angle of attack does not go straight to
a lift/drag lookup. It goes to the dynamic stall model, whose states advance
explicitly from the current inflow conditions and current states — currently
Beddoes-Leishman with Gonzalez's modifications. Lift and drag come from those
dynamic states, and become the loading on the rotor.

The structural side is geometrically exact beam theory via
[GXBeam.jl](https://github.com/byuflowlab/GXBeam.jl), which captures the
geometric nonlinearities exactly while staying cheap for arbitrary section
topology and materials. GXBeam uses constant and linear shape functions, so more
elements are needed for grid independence than a higher-order formulation would
require. Time integration uses GXBeam's undamped Newmark scheme.

### Data flow per time step

```
   env ──► get_aerostructural_velocities ──► (Vx, Vy) in frame L
                     ▲                              │
     deflections,    │                              ▼
     velocities      │                       solve_BEMT  (ImplicitAD boundary)
     interpolated    │                              │
     structural→aero │                              ▼  φ, α, W
                     │                    update_ds_inputs! → DS state march
                     │                              │
                     │                              ▼  cl, cd, cm
                     │                       extract_ds_loads!
                     │                              │
                     │                              ▼
                     │                       dimensionalize!   (Cx,Cy,Cm → Fx,Fy,Mx)
                     │                              │
                     │                              ▼
                     │                       update_forces!    (Gaussian integration
                     │                              │           onto structural elements)
                     │                              ▼
                     │                       GXBeam.step_system!
                     │                              │
                     └──────── update_mesh! ◄───────┘
```

Two-way means each model depends on the other's output. Structure → aero: linear
and angular deflections and linear velocities at every aerodynamic node. Aero →
structure: distributed normal force, tangential force, and axial moment.

The aerodynamic and structural nodes need not be co-located. Deflections are
linearly interpolated from structural to aerodynamic nodes via the
[`InterpolationPoint`](@ref) handles cached in the mesh — the bracket search runs
once at initialization, not every step. Distributed loads go the other way by
Gaussian integration onto the structural elements.

Initialization: the aerodynamic model is brought to steady state first, then the
structural model is initialized from a no-load state.

## AD architecture

Everything is differentiable with ForwardDiff and ReverseDiff, end to end. Three
design decisions make that work.

**1. Custom integrators.** `RK4` and `BDF1` in `solvers.jl` exist so that every
integration step is plain broadcast arithmetic with no opaque internals for duals
to get stuck in. This is the main reason the package does not call
DifferentialEquations.

**2. Type-inferred buffers.** Every `initialize_*` function sizes its buffers
using [`find_inittype`](@ref), which scans representative scalars (chord, twist)
for a `ForwardDiff.Dual` or `ReverseDiff.TrackedReal` and allocates with that
element type. Pass a dual-valued chord and every downstream buffer is wide enough
to carry duals from the start — no retyping mid-march. [`dualcopy`](@ref)
preserves that tracking when arrays are copied.

**3. ImplicitAD at the BEMT root-find.** The inflow-angle solve is wrapped in
`ImplicitAD.implicit` ([bemt.jl:218](https://github.com/byuflowlab/WATT.jl/blob/master/src/bemt.jl)),
so the derivative comes from the implicit function theorem at the converged
solution rather than by propagating duals through every Brent iteration. The
primal root-find runs in `Float64` and only the residual is evaluated with duals.

Profiling confirms this pays off: of the BEMT solve's share of a constraint
Jacobian, the derivative machinery is about 1% while the primal root-find is
about 12%. Note the asymmetry — `BDF1` is *not* wrapped this way and does
propagate duals through its `NLsolve` iterations.

### What is not differentiable

- **No GPU path is differentiable.** Forward pass only.
- **`fixedpoint!` is unverified.** There is no AD test for the steady solver.
  ForwardDiff through `run_sim!` and `run_from_state!` is covered in
  `test/test_ad.jl`; the steady path is not.

### Windowed sensitivities

For optimization, differentiating a full march is often wasteful when the
quantity of interest depends on a short window. [`step_solution!`](@ref) exposes
a single coupled step as a callable primitive, and
[`initialize_from_state`](@ref) / [`run_from_state!`](@ref) march a window from a
frozen starting state — so the gradient covers only the window. See
`examples/step_solution_ad.jl`.

## Future coupling

### Semi-tight coupling

Sub-iterate the aero and structural solves at each time step until they agree,
before advancing time — OpenFAST's tight coupling. The critical detail is that
dynamic stall state must be reset to its start-of-step value at the top of each
sub-iteration, since sub-iteration seeks equilibrium at the current time rather
than advancing stall history. GXBeam's `step_system!` mutates the system in
place, so `system.x` has to be saved and restored too.

The natural AD strategy is not to propagate duals through every sub-iteration but
to converge the loop on stripped `Float64` values and apply the implicit function
theorem at the fixed point — the same `ImplicitAD.implicit` pattern already used
at the BEMT root-find, applied one level up.

### Monolithic coupling

Assemble aerodynamic and structural residuals into a single nonlinear solve per
time step, over a coupled state vector of BEMT inflow angles, dynamic stall
states, and GXBeam displacements and velocities. For NREL 5MW that is roughly 800
unknowns per step, so sparse Jacobians are required.

Neither is implemented. Whether they still should be — given that the structural
surrogate and windowed-AD work reach the same goal by another route — is an open
question tracked in `plan.md`.
