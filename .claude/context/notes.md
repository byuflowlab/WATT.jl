# WATT Development Notes

---

## 2026-03-19 — Gen 1 API Removal and Codebase Cleanup

### Context
The codebase had accumulated a significant amount of dead code from an earlier "Gen 1" simulation API that had been superseded by the current "Gen 2" API (`initialize_sim` / `run_sim!`). A full codebase survey was performed and results were compiled into `CLEANUP.md`, where decisions were made on each item before executing deletions.

### What was done

**Identified and deleted Gen 1 API (aerostructural.jl)**
The original simulation API — `initialize()`, `initial_condition!()`, `take_step!()`, `simulate()`, `simulate!()` — was removed. These were replaced by `initialize_sim()` and `run_sim!()`, which use `GXBeam.initialize_system!()` + `GXBeam.step_system!()` (Newmark-Beta) instead of the Gen 1 `GXBeam.take_step()` (backward Euler on the DAE). The Gen 1 functions had already been commented out pending removal.

**Deleted legacy GXBeam functions (gxbeam.jl)**
- `gxbeam_initial_conditions()` and `gxbeam_initial_conditions!()` — initial condition helpers that only the Gen 1 API called
- `simulate_gxbeam()` — standalone structural-only transient simulation, subsumed by `run_sim!()`
- `steady_simulate_gxbeam()` — standalone structural-only steady simulation, subsumed by `static.jl`

**Deleted legacy utility functions (utils.jl)**
- `prepareextra()`, `saveextra()`, `readextra()` — file I/O helpers from an early development phase, never called
- `plotdshistory()` — debug visualization helper, never called
- `sub_brent()` — custom Brent root-finder translated from OpenFAST, superseded by FLOWMath's Brent solver, never called

**Deleted legacy dynamic stall files**
- `src/beddoesleishman.jl`
- `src/beddoesleishman_aerodyn.jl`
- `src/beddoesleishman_gonzalez.jl`

These were local implementations of three Beddoes-Leishman dynamic stall variants from before the package adopted `DynamicStallModels.jl`. They were never `include()`d in `WATT.jl` and thus already inactive.

**Updated CLAUDE.md**
- Removed `beddoesleishman*.jl` from the file tree
- Updated `aerostructural.jl` description to reflect the active Gen 2 API

### Deferred items (see CLEANUP.md)
Several items were flagged but not acted on:
- `solve_BEM!()` has two implementations with slightly different signatures (with/without `phi0`); needs investigation before consolidating
- `thrusttorque()` / `ccthrusttorque()` in `bem.jl` — dead but deferred
- `derivative_me()` / `mat_derivative()` in `utils.jl` — pre-AD finite-difference utilities, deferred
- `get_aerostructural_velocities()` first overload in `environments.jl` — broken/dead, deferred
- `initialize()` naming collision across `static.jl` and `aero_only.jl` — deferred rename
- `initialize_ds_model()` return count — flagged for verification (Gen 2 code unpacks 4 values)
- `run_sim()` empty stub in `aerostructural.jl` — marked `fix`, not yet implemented

---

## Session — 2026-05-18

- Phase 1 sub-plan focus: green up the existing test suite and stand up a full-simulation regression harness. Working plan at `~/.claude/plans/we-are-in-phase-rustling-eich.md`.
- All five legacy test files now run green: `test_utils.jl`, `test_types.jl`, `test_mesh.jl`, `test_bem.jl`, `test_environments.jl`, plus `test_gxbeam.jl`. Fixed systemic issues: `using OpenFASTsr` → `using OpenFASTTools` everywhere, added missing `using StructArrays`, removed broken `sub_brent` testset, removed `using Revise` from `test_bem.jl`.
- Built `test/simple_NREL5MW.jl` (renamed from `full_test.jl` by user): full-simulation golden-value regression suite at `rtol=1e-8, atol=1e-10` covering `initialize_sim` shape/eltype, `run_sim!` aerostates fields, `gxhistory` spot-checks, no-NaN/Inf, integrated thrust/torque, and idempotency. Wired into `runtests.jl`. Reference JLD2 stores `env` as a `ReconstructedStatic` (SimpleEnvironment holds closures) so the test rebuilds env via `WATT.SimpleEnvironment(...)` from a TurbSim file path.
- Added Phase 3 task 3.6 to project `plan.md`: refactor `SimpleEnvironment` to hold callable interpolant structs (e.g. `Akima`, `Constant`, `TurbulentInflow`) instead of anonymous closures. Same dispatch cost, fixes JLD2 roundtrip, reduces type-parameter explosion.
- `test_bem.jl`: bugs uncovered and fixed in the test itself, not in WATT — endpoint nodes (`rvec[1]==rhub`, `rvec[end]==rtip`) coincide with CCBlade's hub/tip short-circuit so restricted comparison to interior nodes; `twistvec` passed to `CCBlade.Section` in degrees instead of radians made CCBlade converge to a non-physical root. Now asserts numerical field-by-field agreement to `rtol=1e-8`. Also switched to canonical no-`phi0` `solve_BEM!` (the `phi0` overload is dead code — used by nothing internal; flagged for Phase 1 cleanup).
- `test_environments.jl`: stale debug-print scaffolding (`# println`, `# @show`, `# Vxo, Vyo`) and stale failure markers (`#Todo: Failing`, `#, rtol=0.02`) all removed; renumbered cases with explanatory headers (Case 1–12). Two precone blocks were silently passing the *wrong* test — Blade was built without `precone` kwarg while CCBlade reference used it. Added new `@testset "Aerostructural Velocities"` driven by direct `GXBeam.time_domain_analysis`, replacing the long-commented block. Includes root-identity check that `get_aerostructural_velocities` reduces to `get_aero_velocities` at the prescribed (fixed) root.
- `test_gxbeam.jl`: expanded from one trivial test to 197 assertions covering the WATT↔GXBeam interface used by `run_sim!`: `retrieve_eulerangles`, `WMPtoangle`, `get_bladelength_vector`, `update_forces!` (including a linearity check at `rtol=1e-10`).
- `test_mesh.jl`: resurrected commented-out "Interpolation Functions" testset. Fixed undefined `precone`, unqualified `GXBeam.DistributedLoads`/`PrescribedConditions`, strict `==` on SVectors, and `SVector ≈ scalar` comparisons. Dropped the `steady=true` block (deprecated kwarg + fragile assertion). Added `@test_logs (:warn, r"not found")` assertions around the out-of-bounds `find_point_indices` calls — codifies current warn-and-clamp behavior so a future flip to error-on-OOB surfaces as a test failure.
- API observations worth a future session: `WATT.run_sim!` uses `GXBeam.DynamicSystem` + `initialize_system!` + `step_system!`, **not** `GXBeam.time_domain_analysis`. Several tests still use `time_domain_analysis` for setup (it's convenient and produces the same `AssemblyState` type). User OK'd leaving as-is for now; could revisit for path-consistency later.
- `SimpleEnvironment` closure-based design is the underlying reason JLD2 reference data can't roundtrip env directly. The fix (callable structs) is now scheduled in Phase 3.
- New plan files: `~/.claude/plans/we-are-in-phase-rustling-eich.md` (this Phase 1 sub-plan; iteratively revised across the session).
