# WATT.jl Long-Term Rework Plan

> **Status at a glance** — last audited **2026-08-19** against the code, not against commit messages.
>
> | | Phase | State |
> |---|---|---|
> | ✅ | 1 Cleanup · 2 Tests · 3 API · 4 Docs · 11 Windowed AD | Complete |
> | 🟡 | 5 Static (5.3 open) · 8 Polish (2,4,5,6 open) · 9 Surrogate (untested) · 10 GPU (no AD, no tests) | Partial |
> | ⬜ | 6 Tight coupling · 7 Monolithic | Not started — **see the open question at the dependency graph** |
>
> Suite is green: **2156/2156**, zero `@test_broken`. Phases 9–11 were added
> retroactively to describe work that happened off-roadmap from 2026-05-22 onward.

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
**Status:** Completed (2026-05-18)
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
**Status:** Complete (2026-05-19)
**Session estimate:** 1–2 sessions
**Depends on:** Phase 1

### Bugs surfaced (audited 2026-08-19 — see status below)
- ~~`src/aero_only.jl::initialize` calls the outdated 3-positional / kwarg form
  of `initialize_ds_model` and builds a mesh NamedTuple missing `p_ds`.~~
  **RESOLVED in Phase 3** — the `initialize_aero` rename landed together with the
  DS-init signature fix and the missing `p_ds`. The `@test_broken` markers in
  [test/test_aero_only.jl](test/test_aero_only.jl) are now active `@test`.
- AD through `run_sim!` and `fixedpoint!` errors at GXBeam's
  `step_system!` / `steady_state_analysis!`: those write dual values into a
  Float64 `DynamicSystem`/`StaticSystem`.
  **PARTIALLY RESOLVED.** `run_sim!` now differentiates end-to-end —
  [test/test_ad.jl:248](test/test_ad.jl#L248) runs a live
  `ForwardDiff.gradient` over a compliance+mass parameter vector, and
  [test/test_ad.jl:343](test/test_ad.jl#L343) covers the frozen-start window via
  `run_from_state!`. **`fixedpoint!` remains unverified** — `test_ad.jl` has no
  static/steady testset at all, so whether the `steady_state_analysis!` dual
  problem still bites is unknown. Carry this into Phase 5.3.

**The suite currently has zero `@test_broken` and passes 2156/2156 (4m22s).**

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
**Status:** Completed (2026-05-21, commits `39fa6cc` + `3573cc3`; docstring sweep `bea5413`)
**Session estimate:** 1–2 sessions
**Depends on:** Phases 1 + 2

### Audit notes (2026-08-19)
All of 3.1–3.7 verified present in the code. Two deviations worth recording:
- **`thrusttorque` was never written.** The 3.4 export list below names it, but only
  `rotorloads` exists ([src/aero_only.jl:304](src/aero_only.jl#L304)). Either write it
  or drop it from the intended surface — currently it is neither.
- **3.4 has regressed.** The Phase 10 GPU work reintroduced scattered `export`
  statements in six places ([src/bemt_gpu.jl:24](src/bemt_gpu.jl#L24),
  [src/dsmodel_gpu.jl:34](src/dsmodel_gpu.jl#L34), and four in
  [src/aerostructural_surrogate_gpu.jl](src/aerostructural_surrogate_gpu.jl)), taking the
  public surface from ~28 names to 58 — including the implementation-detail constants
  `DS_NSTATES` and `N_BRENT_ITERS_DEFAULT`. Restoring the single centralized export
  block is folded into Phase 4 prerequisites.

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

**3.6 — Refactor `SimpleEnvironment` to drop closure fields**
The current `SimpleEnvironment{TF, F<:Function, G<:Function, ...}` stores eight anonymous closures (`ufun`, `omegafun`, `Vinf`, `RS`, …). Two problems with this:
1. **JLD2 cannot roundtrip closures** — their gensym'd type names don't survive a fresh session. Test fixtures (e.g. [test/simple_NREL5MW.jl](test/simple_NREL5MW.jl)) have to rebuild the env locally rather than load it.
2. **Type-parameter explosion** — eight `<:Function` parameters means every distinct env spawns fresh method specializations, hurting TTFX and invalidating downstream methods.

Refactor each closure into a named callable struct, then parameterize on the concrete struct type:
```julia
struct Constant{T}; val::T; end
(c::Constant)(t) = c.val

struct TurbulentInflow{A}
    Ufit::A; Vfit::A; Wfit::A   # A is e.g. FLOWMath.Akima
end
(t::TurbulentInflow)(s) = SVector(t.Ufit(s), t.Vfit(s), t.Wfit(s))

struct SimpleEnvironment{TF, F, G, ...}
    rho::TF; mu::TF; a::TF; shearexp::TF
    ufun::F          # callable struct, not anonymous closure
    omegafun::G
    # ...
end
```
Performance is unchanged (still a direct dispatch through a concrete field type — same cost as the closure version). Bonus: JLD2 can serialize the whole env now, which lets us save/load reference simulations end-to-end instead of patching env on load.

Migration: the three `environment(...)` constructors in [src/environments.jl](src/environments.jl) get updated to wrap their interpolants/constants in the new callable structs instead of returning closures. No call-site changes needed — `env.ufun(t)` still works.

**3.7 — Docstring audit** (parallel with 3.5)
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
**Status:** Completed (2026-08-19)
**Session estimate:** 1 session
**Depends on:** Phase 3

### What the audit found (2026-08-19)
The docs were not merely incomplete — **the build was red and the tutorial did not run.**
- `apireference.md` was a 9-line stub whose only `@docs` entry was `simulate`, deleted
  in Phase 1. This failed the build.
- `steady.md` was **0 bytes**.
- `developers.md` had two section headers with empty bodies.
- `gettingstarted.md` had six defects making its code non-runnable, including
  `StructArray{DS.Airfoil}` (obsoleted by `5a5f0b4`) and `Rotors.SimpleEnvironment`
  (pre-rename package name).
- `docs/Project.toml` pinned `Documenter = "0.27"`.

### What was done
- **Docstrings:** wrote the 4 missing ones (`BDF1` + its call method, `rotorloads`,
  `DS_NSTATES`, `N_BRENT_ITERS_DEFAULT`) plus 3 internals referenced by public
  docstrings (`dimensionalize!`, `find_inittype`, `wmp_to_angle_dev`). Converted
  `solvers.jl` from `###` headers to the bold-header house style. **All 55 exported
  symbols are now documented**, enforced by `checkdocs = :exports`.
- **Exports:** removed the six scattered `export` statements the GPU work added and
  consolidated them into `src/WATT.jl`, un-exporting `DS_NSTATES` and
  `N_BRENT_ITERS_DEFAULT`. Restores the Phase 3.4 invariant; surface is 55 names.
- **Bug fixed en route:** `rotorloads(rhub, rtip, rvec, loads...)` accumulated with
  `+=` into `undef` memory. Now `zeros`.
- **Pages:** rewrote `apireference.md` (grouped `@docs` + an Internals section),
  `steady.md`, `gettingstarted.md`, `developers.md`, `index.md`; added `gpu.md`
  (split out so no page exceeds Documenter's size threshold). Added three
  walkthroughs — `aeroonly.md`, `gradients.md`, `assembly.md` (building a
  `GXBeam.Assembly` without OpenFAST inputs).
- **Figures:** 7 plots generated from real runs into `docs/src/assets/` and
  embedded in the tutorials. Regenerate with the script recorded in the session
  scratchpad if the physics changes.
- **Also exported** `BladePoints` and `AssemblyPlot` — documented plot recipes
  that had no way to be reached without a `WATT.` prefix. Surface is now 57.
- **Documenter** bumped to 1.x via `Pkg.compat`.

### Deployment (fixed 2026-08-19)
The docs had **never deployed** — all 8 historical Documentation workflow runs
failed at *Install dependencies* in ~20s, and no `gh-pages` branch existed.
Cause: `DynamicStallModels` (a hard dep) is not in the General registry and
`*Manifest.toml` is gitignored, so a clean runner could not resolve it.
- `.github/workflows/documentation.yaml` now `Pkg.add`s it by URL before
  `instantiate`; actions bumped to `checkout@v4` / `setup-julia@v2` (the old ones
  hit the Node 20 removal in Sept 2026) and Julia to 1.11, plus `cache@v2`.
- `deploydocs` gained `push_preview = true`, which `cleanup.yaml` already assumed.
- `gh-pages` created and Pages enabled: **http://flow.byu.edu/WATT.jl/** (the
  domain comes from the org, matching GXBeam.jl — no CNAME file needed).
- **Still required:** the workflow only triggers on `master`, so this must be
  merged there before anything deploys.

### Verification
- `julia --project=docs/ docs/make.jl` — **clean, zero warnings.**
- Every code block in `steady.md` and `gettingstarted.md` was **executed** against the
  `test/fixtures/nrel5mw.jl` setup, not just read. The convergence table and tip
  deflection figures in `steady.md` are measured output.
- `Pkg.test()` — **2156/2156 pass** after the export and `rotorloads` changes.

### Deferred
Literate.jl examples (mentioned in the original goals) were not added. The runnable
tutorials plus the scripts in `examples/` cover the same ground for now.

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
**Status:** Partial — 5.1 and 5.2 complete; **5.3 not started**
**Session estimate:** 1–2 sessions (≈0.5 remaining)
**Depends on:** Phase 1 bug fixes (can begin in parallel with Phase 3)

### Audit notes (2026-08-19)
- **5.1 done.** `fixedpoint!` converges to ~1e-14 by iteration 10 on NREL 5MW rated.
- **5.2 done.** [test/test_static.jl](test/test_static.jl) covers the shape/eltype
  contract, the iteration-1 match against direct CCBlade at the undeflected position
  (rtol 1e-12), and iteration-10 convergence.
- **5.3 not started.** `static_ic` appears nowhere in `src/` — the converged static
  solution still cannot warm-start the transient solver.
- **Carried in from Phase 2:** AD through `fixedpoint!` has no test coverage. Verify it
  as part of 5.3, since a static warm-start that breaks the gradient chain is worse than
  no warm-start.

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
**Status:** Partial — items 1 and 3 landed early, out of order
**Session estimate:** 1 session (≈0.6 remaining)
**Depends on:** Phases 3–7

### Audit notes (2026-08-19)
- **Item 1 (plot recipes) — done, and done differently than planned.** Rather than the
  `Requires.jl` + `__init__` route below, the recipes went in as plain `RecipesBase`
  recipes in [src/gxbeam.jl:324](src/gxbeam.jl#L324) and `:375`, landing during Phase 1.
  `RecipesBase` is a hard dep but a featherweight one, and `Plots` is not a dep at all —
  this is the better outcome; the `Requires.jl` text below is superseded.
- **Item 3 (benchmarks) — done.** Benchmark ladder in `8e10189`; the BEMT profiling
  result and the negative warm-start finding are written up in
  [BEMT_PERF_HANDOFF.md](BEMT_PERF_HANDOFF.md) and
  [WARMSTART_BRENT_FINDINGS.md](WARMSTART_BRENT_FINDINGS.md).
- **Still open:** item 2 (`PrecompileTools`), item 4 (allocation audit), item 5 (compat
  bounds), item 6 (version bump — still `0.3.0`, no `CHANGELOG.md`).

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

## Phases 9–11: Work that happened off-roadmap

*Added 2026-08-19. Phases 1–8 were written 2026-05-18. From 2026-05-22 onward the work
went somewhere the roadmap never anticipated — 13 of the 18 commits since, and five new
source files. These phases are recorded retroactively so the roadmap describes the
package that actually exists.*

---

## Phase 9: Structural Surrogate (Koopman)
**Status:** Working, under-tested
**Landed:** 2026-05-25 → 2026-07-30 (`8bec1cc`, `384b902`, `c62f1d7`, `ddac5c4`)

Replaces the GXBeam structural solve with a learned Koopman-style latent model:
encode the initial structural state, march a linear latent step, decode back to
deflections. Cuts the per-step structural cost to a matrix multiply, which is what makes
the long marches in the downstream optimization harness affordable.

- **Code:** [src/aerostructural_surrogate.jl](src/aerostructural_surrogate.jl);
  `SurrogateMesh`, `SurrogatePointState`, `SurrogateAssemblyState`
  ([src/types.jl:300-318](src/types.jl#L300-L318)); the
  `AbstractStructuralSurrogate` interface — `encode_initial`, `step_latent`,
  `decode`, `decode!`.
- **Entry points:** `initialize_sim_surrogate`, `run_sim_surrogate!`, `run_sim_surrogate`.
- **Gap:** there is **no `test/test_surrogate.jl`**. The surrogate path is exercised only
  through `examples/`, so nothing in CI catches a regression in it. Highest-value
  outstanding work in this phase.
- **Also:** p-conditioned variant demoed in
  [examples/aerostructural_nrel5mw5seg_surrogate_pcond.jl](examples/aerostructural_nrel5mw5seg_surrogate_pcond.jl).

---

## Phase 10: GPU Backends
**Status:** Forward pass validated on hardware; no AD
**Landed:** 2026-07-10 → 2026-07-17 (`cc0b495` → `100949c`)

KernelAbstractions ports of the aero stack, for batched evaluation across many
simulations at once.

- **Code:** [src/bemt_gpu.jl](src/bemt_gpu.jl) (conditional q1+q3 kernel, fixed-iteration
  Brent), [src/dsmodel_gpu.jl](src/dsmodel_gpu.jl) (ADG port),
  [src/aerostructural_surrogate_gpu.jl](src/aerostructural_surrogate_gpu.jl) (coupled
  forward pass with the Phase 9 surrogate as the structural model).
- **Validated:** GPU-BEMT on an H200 cluster node; DS model against the CPU backend;
  three-way coupled comparison in
  [examples/gpu_aerostructural_benchmark.jl](examples/gpu_aerostructural_benchmark.jl).
- **Backend-generic** via a generic `similar_type`, so no per-script glue is needed.
- **Gaps:** forward pass only — no AD through any GPU path. Exports need narrowing
  (see the Phase 3 audit note). No GPU tests in the suite; validation lives in
  `examples/` and SLURM output files.

---

## Phase 11: Windowed AD / Single-Step Primitive
**Status:** Complete, with AD coverage
**Landed:** 2026-07-29 (`75984b6`)

Exposes one coupled time step as a callable primitive, so a gradient can be taken over a
short window from a frozen start rather than over an entire march — the enabling piece
for cheap windowed sensitivities in optimization.

- **Code:** `step_solution!`, `initialize_from_state`
  ([src/aerostructural.jl:508](src/aerostructural.jl#L508)), `run_from_state!`.
- **Tests:** [test/test_step_solution.jl](test/test_step_solution.jl) gates
  equivalence against the full march; [test/test_ad.jl:343](test/test_ad.jl#L343)
  covers `ForwardDiff` through the frozen-start window.
- **Note:** this quietly did much of what Phase 6 was meant to enable. See the open
  question below.

---

## Phase Dependency Graph

```
Phase 1 (Cleanup) ✓
    └─> Phase 2 (Tests) ✓
            ├─> Phase 3 (API) ✓ ───────────────> Phase 4 (Docs) ✓
            │
            └─> Phase 5 (Static)  [5.1 ✓ 5.2 ✓ | 5.3 open]
                    └─> Phase 6 (FP-in-time)     ← not started
                            └─> Phase 7 (Monolithic)  ← not started
                                        └─> Phase 8 (Polish)  [1 ✓ 3 ✓ | 2,4,5,6 open]

Off-roadmap track (2026-05-22 →):
    Phase 9 (Koopman surrogate) ──> Phase 10 (GPU backends)
    Phase 11 (Windowed AD)  [independent]
```

### ⚠ Open question — is the Phases 6→7 route still wanted?

Phases 6 (tight coupling) and 7 (monolithic) were designed as the route to cheap,
accurate gradients for optimization. Since then the off-roadmap track went after the
same goal by a different road: Phase 9 makes the structural solve cheap enough that
sub-iteration may not be the bottleneck it was assumed to be, Phase 11 delivers windowed
sensitivities without a monolithic residual, and Phase 10 attacks throughput instead of
per-step cost.

Phase 7 was flagged **high risk / 2–3 sessions** when written. Before spending that,
decide explicitly:

1. **Still needed** — the surrogate is an approximation; monolithic coupling is the
   high-fidelity reference the surrogate must be validated against. Keep 6→7.
2. **Superseded** — the surrogate + windowed-AD route is the real path forward. Demote
   6→7 to "someday", and redirect the effort into Phase 9's missing test coverage and
   AD through the GPU path.
3. **Partially** — do Phase 6 (tight coupling, the cheaper half) as the validation
   reference, and drop Phase 7.

*Not resolved as of 2026-08-19. This blocks nothing currently in flight, but it should
be answered before the next multi-session commitment.*

---

## Critical Files Reference

*Refreshed 2026-08-19 — the Phase 1 bugs listed here originally are all fixed, and
`bem.jl` was renamed `bemt.jl` in `bea5413`.*

| File | Role | Current concern |
|------|------|----------------|
| `src/aerostructural.jl` | Main transient API + windowed-AD primitive | Home for Phases 6+7 if pursued; hosts Phase 11 |
| `src/static.jl` | Static fixed-point API | Phase 5.3 target; AD through `fixedpoint!` unverified |
| `src/bemt.jl` | BEMT solver | ImplicitAD boundary; primal Brent is the remaining hot spot (see `WARMSTART_BRENT_FINDINGS.md`) |
| `src/gxbeam.jl` | Structural wrapper | Also holds the `RecipesBase` plot recipes (Phase 8.1) |
| `src/mesh.jl` | Coordinate transforms + coupling | Typed meshes done; BC→HR→G→L chain needs documenting (Phase 4) |
| `src/dynamicstallmodels.jl` | DS state management | DS reset for sub-iteration, if Phase 6 goes ahead |
| `src/environments.jl` | Wind field + callable-struct env | Clean since Phase 3.6 |
| `src/aerostructural_surrogate.jl` | Koopman structural surrogate | **No test coverage** (Phase 9) |
| `src/bemt_gpu.jl`, `src/dsmodel_gpu.jl`, `src/aerostructural_surrogate_gpu.jl` | GPU backends | Scattered exports; no AD; no tests (Phase 10) |

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
