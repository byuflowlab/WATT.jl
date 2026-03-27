# WATT.jl Long-Term Rework Plan

## Deliverables of This Planning Session

These are created once this plan is approved — before any Phase 1 work begins:

1. **`plan.md`** — Save this plan to the project root (`/Users/adamcardoza/.julia/dev/WATT/plan.md`) as the canonical long-term reference
2. **`.claude/context/GXBEAM_STYLE.md`** — Explore the GXBeam.jl source on GitHub; capture struct conventions, dispatch patterns, public/private API separation, and docstring style in a reference file loaded in future sessions
3. **`.claude/commands/watt-subplan.md`** — A skill that, given a phase number, reads `plan.md`, reads the relevant current code, and produces a detailed session-scoped implementation sub-plan with tasks, test steps, and exit criteria for that phase

---

## Context

WATT.jl is a research-grade Julia package for nonlinear unsteady aeroelastic simulation of wind turbine blades. It couples BEMT aerodynamics (CCBlade), dynamic stall (DynamicStallModels), and geometrically exact beam structures (GXBeam) with full AD compatibility (ForwardDiff + ReverseDiff + ImplicitAD). It was built quickly for research and works correctly, but has significant technical debt: dead code, missing tests for main solvers, no exports, naming collisions, and scattered coordinate transforms.

**Goals of this rework:**
1. Clean up technical debt without breaking working functionality
2. Add comprehensive testing (unit, integration, AD, regression)
3. Establish a clean GXBeam-style public API
4. Complete documentation with Literate examples
5. Bring the static fixed-point solver (`fixedpoint!`) fully online
6. Add a fixed-point-in-time tight coupling solver
7. Add a monolithic aero-structural coupled solver

**Guiding principles:**
- Keep original code while building replacements alongside it
- Pattern after [GXBeam.jl](https://github.com/byuflowlab/GXBeam.jl): dispatch over type checks, clean public/private API, structured docstrings
- Follow `.claude/context/JULIA_GUIDE.md` for style; see `.claude/context/GXBEAM_STYLE.md` (created Phase 1) for GXBeam-specific patterns
- AD compatibility must be preserved throughout all refactoring

---

## Phase 1: Foundation — Code Cleanup
**Status:** Not started
**Session estimate:** 1 session

### Goal
Remove all dead code, fix active bugs, clean up debug output, trim unused dependencies. Package should build and pass existing tests. No user-facing API changes.

### Entry criteria
Existing tests are not reliable (written hastily). Instead: the user will provide a known-correct simulation script; run it and save its key outputs (thrust, torque, tip deflection, a few aero state values) as JSON golden files in `test/reference/`. These become the regression baseline for all subsequent phases.

### Exit criteria
All existing tests pass; no `@show`/`println` in hot paths; dead code removed; `Plots` dependency commented out.

### Pre-work: create GXBeam style guide
Before cleanup begins, explore [GXBeam.jl source](https://github.com/byuflowlab/GXBeam.jl) and create `.claude/context/GXBEAM_STYLE.md` capturing:
- How types are structured (parametric structs, field naming, inner constructors)
- How public vs. private API is separated (naming, exports, module layout)
- How multi-dispatch is used instead of `if/isa` branching
- Docstring format and section headings used in that codebase
- How `RecipesBase` / conditional Plots loading is handled (if applicable)

This file is loaded into context for all subsequent sessions.

### Active bugs to fix first (sequential)
| File | Line | Issue |
|------|------|-------|
| `src/bem.jl` | ~194 | `@show typeof(phistar)` is ACTIVE (not commented) — fires on every BEM solve |
| `src/bem.jl` | ~274 | `println("Vx and Vy is zero.")` is ACTIVE |
| `src/aerostructural.jl` | ~371 | References `gxhistory_new` which is never defined; should be `gxhistory` |
| `src/static.jl` | ~30 | `blade.airfoils[1].c` — `Airfoil` has no `.c` field; should be `blade.c[1]` |
| `src/environments.jl` | ~205 | `get_aerostructural_velocities(env, aerov, ...)` calls nonexistent function signature; add `@warn` |

### Dead code to delete (can parallelize)
| File | What to delete |
|------|---------------|
| `src/aerostructural.jl` | All Gen 1 commented blocks (~400 lines): `initialize`, `initial_condition!`, `take_step!`, `simulate`, `simulate!` |
| `src/gxbeam.jl` | `gxbeam_initial_conditions`, `gxbeam_initial_conditions!`, `simulate_gxbeam`, `steady_simulate_gxbeam` (~180 lines) |
| `src/utils.jl` | `prepareextra`, `saveextra`, `readextra`, `plotdshistory`, `sub_brent`, commented Brent block (~100 lines) |
| `src/solvers.jl` | `DBDF1!`, `DBDF2!`, `fixedpointbem`, `secant` (~165 lines). Keep `RK4` and `BDF1`. |
| `src/types.jl` | Commented-out `AeroStates` struct and `get_aerostate` (~75 lines), alternate `Blade` struct |
| `src/bem.jl` | `show_dual_vec` debug function; audit duplicate `solve_BEM!` — keep the no-`phi0` version as canonical, rename the `phi0` version to `solve_BEM_warm!` or delete |

### Dependency cleanup
- **Comment out** `using Plots` and the `plotpoints`/`plotassembly` functions in `src/WATT.jl` temporarily (do not delete — will be refactored in Phase 8 via `RecipesBase` + `Requires.jl`)
- Remove `CurveFit` from `using` and `Project.toml` (never called)
- Keep `NLsolve` (needed by `BDF1`), keep `DelimitedFiles` (used in `environments.jl`)
- Remove `using Revise` from `test/test_bem.jl` (wrong place, breaks CI)

### Tests needed
No new tests — just verify existing test suite still passes after each deletion.

---

## Phase 2: Test Infrastructure
**Status:** Not started
**Session estimate:** 1–2 sessions
**Depends on:** Phase 1

### Goal
Build comprehensive test suite: unit, integration (all three solver modes), and AD compatibility.

### Exit criteria
All new tests pass; ForwardDiff and ReverseDiff verified on at least one end-to-end path; `test_bem.jl` no longer imports `Revise`.

### New test files to create

| File | What it tests |
|------|--------------|
| `test/test_dynamicstall.jl` | `initialize_ds_model`, `update_ds_inputs!`, `update_ds_states!`, `extract_ds_loads!`; AD transparency |
| `test/test_transforms.jl` | Coordinate transforms BC→HR→G→L: identity case, inverse, known-value at azimuth=π/2 |
| `test/test_aero_only.jl` | `initialize_aero` + `simulate!`: shapes, no NaN/Inf, loads match BEM reference at t=0, golden regression values |
| `test/test_aerostructural.jl` | `initialize_sim` + `run_sim!`: structural tip deflection bounded, aero loads have no NaN/Inf, regression golden values |
| `test/test_static.jl` | `initialize_static` + `fixedpoint!`: first iteration matches CCBlade, convergence check (may fail initially — that drives Phase 5) |
| `test/test_ad.jl` | ForwardDiff through `simulate!` (aero-only) and `run_sim!`; ReverseDiff through `simulate!`; compare against finite difference to 1e-4 rtol |

### Notes
- Use NREL 5MW at rated conditions (10 m/s, 12.1 RPM) as the standard test case
- Use short time vectors (5–10 steps) for basic correctness; longer only for regression golden values
- AD tests are expensive — gate with `ENV["WATT_AD_TESTS"] == "true"` so CI can skip them
- Verify `testing/OpenFAST_NREL5MW` data path exists in repo; if not, switch to `data/`

### Existing test fixes
- `test_bem.jl`: remove `using Revise`
- `test_utils.jl`: remove tests for `sub_brent` (being deleted in Phase 1)
- Verify all relative data paths are consistent

---

## Phase 3: API Reorganization
**Status:** Not started
**Session estimate:** 1–2 sessions
**Depends on:** Phases 1 + 2

### Goal
GXBeam-style clean API: proper exports, naming collision fixes, typed structs for `mesh` and `aerostates`, docstrings on all public functions.

### Exit criteria
All symbols callable without `WATT.` prefix after `using WATT`; no ambiguous dispatch; `mesh` is a typed struct.

### Tasks

**3.1 — Resolve naming collisions (do first)**
- `aero_only.jl`: `initialize()` → `initialize_aero()`
- `static.jl`: `initialize()` → `initialize_static()`
- `aerostructural.jl`: keep `initialize_sim`
- Update all call sites in tests and docs

**3.2 — Define `SimMesh` struct** (`src/mesh.jl` or new `src/structs.jl`)
- Replace the `mesh` NamedTuple with a proper parametric struct `SimMesh{...}`
- Create `AeroMesh` variant for aero-only (no GXBeam fields)
- Improves type stability and documents the contract

**3.3 — Define `AeroStates` struct** (`src/types.jl`)
- Replace `aerostates` NamedTuple with an immutable parametric struct `AeroStates{TF}` containing `azimuth`, `phi`, `alpha`, `W`, `Cx`, `Cy`, `Cm`, `Fx`, `Fy`, `Mx`, `xds`
- A regular (immutable) struct is sufficient: the fields hold arrays, and we mutate the array *contents* in-place — not the array references themselves. `mutable struct` would only be needed if we needed to reassign the fields to point at different arrays.
- Add `StaticAeroStates{TF}` for the static solver (1D vectors, no time dimension)

**3.4 — Establish clean exports** (`src/WATT.jl`)
```julia
# Types
export Rotor, Blade, RK4, BDF1
export AeroStates, StaticAeroStates, SimMesh, AeroMesh

# Environment
export environment

# Simulation entry points
export initialize_aero, simulate!
export initialize_sim, run_sim!, run_sim
export initialize_static, fixedpoint!

# Post-processing
export rotorloads, thrusttorque
```
Remove scattered `export` statements from individual files.

**3.5 — Implement `run_sim()` non-mutating wrapper**
```julia
function run_sim(rotor, blade, assembly, env, tvec; kwargs...)
    aerostates, gxhistory, mesh = initialize_sim(blade, assembly, tvec)
    run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; kwargs...)
    return aerostates, gxhistory, mesh
end
```

**3.6 — Docstring audit** (parallel with 3.5)
Priority order: `run_sim!` → `initialize_sim` → `simulate!` → `initialize_aero` → `fixedpoint!` → `initialize_static` → `take_aero_step!` → `Blade`, `Rotor` → `SimpleEnvironment`

Format (GXBeam-style):
```julia
"""
    function_name(arg1, arg2; kwarg=default) -> out1, out2

One-sentence description.

**Arguments**
- `arg1::Type`: description

**Returns**
- `out1::Type`: description

**Notes**
Caveats, AD compatibility, performance notes.
"""
```

---

## Phase 4: Documentation
**Status:** Not started
**Session estimate:** 1 session
**Depends on:** Phase 3

### Goal
All docs pages complete; `docs/make.jl` builds without warnings; all three simulation modes have tutorial coverage.

### Tasks
| File | Action |
|------|--------|
| `docs/src/apireference.md` | Replace stub with `@docs` blocks for all exported symbols, organized by category |
| `docs/src/steady.md` | Write full tutorial for static fixed-point solver (depends on Phase 5 for final numbers) |
| `docs/src/developers.md` | Add: reference frame diagram (BC→HR→G→L), data flow diagram, AD architecture, coupling strategy details |
| `docs/src/gettingstarted.md` | Update: rename `initialize()` → `initialize_aero()`, add section showing `initialize_sim` + `run_sim!`, add brief note on static analysis |
| `docs/make.jl` | Uncomment `steady.md`; add Examples section when Literate examples are added |

---

## Phase 5: Static Solver Integration
**Status:** Not started
**Session estimate:** 1–2 sessions
**Depends on:** Phase 1 bug fixes (can begin in parallel with Phase 3)

### Goal
`fixedpoint!` fully functional, tested, and integrated. Provides warm-start for transient solver.

### Exit criteria
`fixedpoint!` converges on NREL 5MW at rated conditions; `test_static.jl` passes; loads compare within a few percent of time-averaged transient.

### Tasks

**5.1 — Debug and fix `fixedpoint!`** (`src/static.jl`)
Known issues to investigate after Phase 1 fixes:
- Verify `GXBeam.steady_state_analysis!` return signature matches what `fixedpoint!` expects
- Verify `AssemblyState` type from `steady_state_analysis!` is compatible with `update_mesh!`
- The `SimMesh` for static case is missing DS fields (`y_ds`, `p_ds`) — verify `update_mesh!` doesn't require them
- Add convergence tracking: `norm(Fx - Fx_prev) / norm(Fx_prev)` each iteration; return convergence history

**5.2 — Convergence validation** (`test/test_static.jl`)
- First iteration should match CCBlade at undeflected position
- Residual should drop toward tolerance over iterations
- Final loads should match time-averaged transient within a few percent

**5.3 — Integrate as transient initial condition** (`src/aerostructural.jl`)
```julia
function run_sim!(...; static_ic=nothing, ...)
    # if static_ic provided, use (aerostates_s, gxstates_s) to warm-start
end
```

---

## Phase 6: Fixed-Point Iteration in Time (Tight Coupling)
**Status:** Not started
**Session estimate:** 1–2 sessions
**Depends on:** Phase 5

### Goal
Partitioned aero-structural sub-iteration at each time step — iterate until convergence before advancing time. Analogous to OpenFAST tight coupling.

### Exit criteria
`run_sim_tight!` passes convergence tests; AD-compatible; benchmarked.

### Architecture
Implement as separate `run_sim_tight!` function (not a keyword on `run_sim!`) for clean separation during development.

**Sub-iteration structure at each time step:**
```julia
for sub in 1:max_subiter
    take_aero_step!(...)       # aero using current structural state
    update_forces!(...)        # dimensionalize, map to GXBeam elements
    GXBeam.step_system!(...)   # structural response
    update_mesh!(...)          # update aero mesh with new deflections

    if norm(Fx - Fx_prev) / (norm(Fx_prev) + eps()) < sub_atol
        break
    end
    Fx_prev .= Fx
    xds .= xds_old             # reset DS state for next sub-iteration
end
```

**Critical detail:** DS state must be reset to `xds_old` at the start of each sub-iteration (sub-iteration finds equilibrium at current time, not advancing DS history).

**AD strategy for sub-iteration — use `ImplicitAD.implicit`:**
Rather than propagating dual numbers through every sub-iteration (expensive and unnecessary), wrap the converged fixed-point step as an implicit function and apply the implicit function theorem (IFT) for derivatives:
1. Strip duals: converge the sub-iteration loop using `ForwardDiff.value.(x)` (pure Float64, no dual overhead)
2. At the converged state, re-evaluate the single step once with full dual numbers to form the residual Jacobians
3. Apply IFT: `∂x_converged/∂p = -(∂F/∂x)⁻¹ (∂F/∂p)`

`ImplicitAD.implicit(solve_fp, residual_fp, p)` handles this automatically. The sub-iteration `solve_fp` receives stripped Float64 inputs; `residual_fp` is evaluated with duals only at the solution point.

**GXBeam detail:** `step_system!` modifies `system` in-place. Save and restore `system.x` at start of each sub-iteration.

### Tests needed (`test/test_tight_coupling.jl`)
- With rigid blade: 1 sub-iteration should suffice (no structural feedback)
- With flexible blade: 3 sub-iterations should be closer to exact coupling than 1
- Benchmark: 1, 2, 5, 10 sub-iterations vs wall clock
- AD test: ForwardDiff through `run_sim_tight!` with 2 sub-iterations

---

## Phase 7: Monolithic Coupled Solver
**Status:** Not started
**Session estimate:** 2–3 sessions
**Depends on:** Phase 6
**Risk:** High — most technically complex phase

### Goal
Single nonlinear solve assembling both aerodynamic and structural residuals simultaneously per time step. Highest fidelity for gradient-based optimization.

### Architecture
Coupled state vector per time step:
```
x = [phi_1..phi_na,       # BEM inflow angles (na elements)
     xds_1..xds_ns,       # DS model states   (ns = 4*na)
     xgx_1..xgx_ngx,      # GXBeam DOFs       (~18*ne)
     xgxdot_1..xgxdot_ngx # GXBeam velocities (~18*ne)]
```

For NREL 5MW (17 aero nodes, 20 structural elements): ~805 unknowns/step — sparse Jacobians required.

**Residual assembly** (`src/aerostructural.jl` or new `src/coupling.jl`):
```julia
function residual_coupled!(R, x, p, blade, rotor, env, mesh, t, dt, xds_old, xgx_old)
    phi, xds, xgx, xgxdot = unpack_coupled_state(x, ...)
    update_mesh_from_state!(mesh, xgx, ...)
    # BEM residuals via CCBlade.residual_and_outputs(phi[j], ...)
    # DS residuals: xds - RK4_step(xds_old, y_ds_new, p_ds, dt)
    # GXBeam residuals via GXBeam.system_residual(...)
end
```

**Solver:** Use `NonlinearSolve.jl` (better AD support than `NLsolve.jl`, more actively maintained) with sparse Jacobian via `SparseDiffTools.jl`.

**AD:** Wrap with `IAD.implicit_unsteady` treating coupled state as the implicit variable, or hand-code the `rrule`.

### Validation (`test/test_monolithic.jl`)
- For simple test case (no DS dynamics, 2–3 steps): monolithic must match tight partitioned
- Benchmark: monolithic vs. tight vs. loose — wall clock and derivative computation time

### Key risks
- GXBeam residual function may not be directly accessible in the right form — check GXBeam source
- `CCBlade.residual_and_outputs` must be differentiable through; ImplicitAD wrapping already exists at BEM level but not at coupled level
- DS state reset in sub-iteration: same issue as Phase 6 but inside the residual function

---

## Phase 8: Package Polish
**Status:** Not started
**Session estimate:** 1 session
**Depends on:** Phases 3–7

### Tasks
1. **Plots extension via `RecipesBase` + `Requires.jl`** — Restore the commented-out `plotpoints`/`plotassembly` functions using conditional loading:
   ```julia
   # src/WATT.jl
   function __init__()
       @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plots_ext.jl")
   end
   ```
   Move plot functions to `src/plots_ext.jl`; implement as `RecipesBase` recipes so they work with any Plots backend without Plots being a hard dependency.
2. **Precompilation** — Add `PrecompileTools.jl`, run minimal NREL 5MW setup through `@compile_workload` for the three solver modes. Target: load time < 1s after precompile.
3. **Benchmarks** (`benchmarks/bench_aerostructural.jl`) — Time per step breakdown: `take_aero_step!`, `GXBeam.step_system!`, `update_mesh!`; allocation count per step (goal: zero in hot path); loose vs. tight vs. monolithic for 10 revolutions.
4. **Allocation audit** — Known sources: `dualcopy()` in `run_sim!`, `GXBeam.DistributedLoads` construction per element per step, DS state copy for sub-iteration. Fix if possible.
5. **Final dependency audit** — Tighten compat bounds in `Project.toml` for `CCBlade`, `GXBeam`, `DynamicStallModels`, `ImplicitAD`.
6. **Version bump** — Bump to `0.4.0` in `Project.toml` (breaking API changes from Phase 3). Write `CHANGELOG.md`.

---

## Phase Dependency Graph

```
Phase 1 (Cleanup)
    └─> Phase 2 (Tests)
            ├─> Phase 3 (API) ─────────────────> Phase 4 (Docs)
            │
            └─> Phase 5 (Static)  [can start in parallel with Phase 3 after Phase 1]
                    └─> Phase 6 (FP-in-time)
                            └─> Phase 7 (Monolithic)
                                        └─> Phase 8 (Polish)  [partial overlap with Phase 7]
```

---

## Critical Files Reference

| File | Role | Primary Concern |
|------|------|----------------|
| `src/aerostructural.jl` | Main transient API | Active bug (gxhistory_new), home for Phases 6+7 |
| `src/static.jl` | Static fixed-point API | Active bug (blade.airfoils[1].c), Phase 5 target |
| `src/bem.jl` | BEMT solver | Active @show debug, duplicate function, ImplicitAD boundary |
| `src/gxbeam.jl` | Structural wrapper | ~180 lines dead code, Plots dependency link |
| `src/mesh.jl` | Coordinate transforms + coupling | Untested transforms (Phase 2), NamedTuple→struct (Phase 3) |
| `src/dynamicstallmodels.jl` | DS state management | Untested (Phase 2), DS reset for sub-iteration (Phase 6) |
| `src/environments.jl` | Wind field + velocity transforms | Broken function at line 205 |

---

## Testing Robustness Checklist (for any phase)

A well-tested physics package needs tests at these levels:

- [ ] **Unit tests** — each public function in isolation with known inputs/outputs
- [ ] **Physics identity tests** — mathematical invariants: rotation orthogonality, zero-shear inflow, rigid-body limit
- [ ] **Integration tests** — all three solver modes run end-to-end without error on NREL 5MW
- [ ] **Regression tests** — golden values stored for key outputs; fail if results drift
- [ ] **AD tests** — ForwardDiff and ReverseDiff agree with finite differences to 1e-4 rtol on at least one path per solver mode
- [ ] **Convergence tests** — iterative solvers (fixedpoint!, run_sim_tight!, monolithic) converge at expected rate
- [ ] **Boundary/edge case tests** — zero wind speed, zero structural stiffness, single blade node, etc.
- [ ] **Type stability tests** — `@code_warntype` clean or `JET.jl` no errors on hot paths

---

*This plan is the overarching roadmap. Each session should create a focused sub-plan derived from one phase.*
