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

---

## Session — 2026-05-18 (Phase 1 finish)

- Closed out the remaining Phase 1 cleanup. Plan at `~/.claude/plans/alright-let-s-make-a-calm-nest.md`. Full test suite (incl. the 100 s `simple_NREL5MW` regression) passes — no physics drift.
- Active bugs deleted in `src/bem.jl`: `@show typeof(phistar)` after the IAD.implicit call, `println("Vx and Vy is zero.")` short-circuit, and the entire `show_dual_vec` debug helper.
- `src/static.jl:30`: `find_inittype(blade.airfoils[1].c, …)` → `find_inittype(blade.c[1], …)`. `Airfoil` never had a `.c` field; this path was latent-dead until now.
- `src/environments.jl`: deleted the unused 9-arg `get_aerostructural_velocities(env, aerov, t, r, azimuth, precone, tilt, yaw, hubht)` overload. Verified zero callers; the rotor/blade-first overload at the bottom of the file is the only live one.
- `src/solvers.jl`: removed `DBDF1!`, `DBDF2!`, `fixedpointbem`, `secant` (~165 lines). `RK4` and `BDF1` retained — `BDF1` is why `NLsolve` stays a dep.
- `src/types.jl`: removed commented-out alternate `Blade` struct, the `AeroStates` struct, and `get_aerostate`. A commented `Mesh` struct further down was left in place (out of scope for the approved plan; flag for next pass).
- **Plots → RecipesBase**: dropped `Plots` from `Project.toml`, added `RecipesBase`. `using Plots` in `src/WATT.jl` replaced with `using RecipesBase`. `plotpoints` and `plotassembly` rewritten as `@recipe` blocks on new wrapper types `BladePoints` and `AssemblyPlot` (in `src/gxbeam.jl`), both exported from `WATT`. Wrapper-types pattern was chosen explicitly to avoid type piracy on `GXBeam.Assembly`. User call pattern: `using Plots; plot(AssemblyPlot(assembly))`.
- Recipe smoke test was deferred — neither recipe is exercised by the regression suite. Add to Phase 2 backlog if Plots compatibility matters.
- `solve_BEM!` phi0-variants (bem.jl ~lines 47–206) intentionally kept per user decision; rename to `solve_BEM_warm!` deferred.
- The user ran `Pkg.rm("Plots"); Pkg.add("RecipesBase")` themselves (foreground task in their REPL); I background-ran it once but stopped it on their request — Pkg ops should run in their REPL, not via subprocess.

---

## Session — 2026-05-19 (Phase 2: Test infrastructure)

- Phase 2 sub-plan at `~/.claude/plans/2-expressive-bunny.md`. Full test suite green (plain `Pkg.test()` and `WATT_AD_TESTS=true`). 2 `@test_broken` markers persist for the aero-only `simulate!` path — recorded in `plan.md` Phase 2 as deferred.
- New test files: `test/fixtures/nrel5mw.jl` (shared NREL 5MW blade/rotor/env/assembly fixture for new tests; legacy tests keep their inline setup), `test/test_dynamicstall.jl` (240 assertions on the DS wrapper layer), `test/test_aero_only.jl` (mostly broken/skipped, documents the broken state), `test/test_aerostructural.jl` (short-tvec `initialize_sim`/`run_sim!` coverage + idempotency), `test/test_static.jl` (49 assertions on `initialize` and `fixedpoint!`), `test/test_ad.jl` (BEM-level FwdDiff/RevDiff + run_sim! gradient over compliance/mass — gated by `ENV["WATT_AD_TESTS"]`).
- `src/aero_only.jl:71` narrow fix: `find_inittype(blade.airfoils[1].c, …)` → `find_inittype(blade.c[1], …)`. Same `.c` bug as `static.jl` in Phase 1.
- Two pre-existing baseline issues surfaced and resolved by user: (1) `Statistics` stdlib missing from `test/Project.toml` — added via `Pkg`. (2) `WATT/Manifest.toml` (project root) pinned `DynamicStallModels` to an old git rev (`adam` branch, v0.1.0) that predated `Discrete`; this shadowed the dev path inside `Pkg.test()`'s sandbox even though direct `--project=test` runs worked. User dev-linked DS at the WATT root.
- **AD pattern landed (this is the load-bearing finding)**: end-to-end AD through `run_sim!` works only via the GXBeam parameter-vector approach mirrored from Cardoza2026's `fatigue_analysis!`. Pack per-element compliance + mass into `p`, define `pfunc(p,t)` that rebuilds a dual-typed `Assembly` + `distributed_loads` on every step, and a `prepp!(p, Fx, Fy, Mx)` callback that writes the dual-typed aero outputs into `p`'s load slot. Pass the non-dual skeleton `assembly` to `initialize_sim` for preallocation only; the dual assembly comes from `pfunc`. Also promote `blade.c` / `blade.twist` to the dual eltype so WATT's `find_inittype` allocates dual `gxhistory`/`aerostates` storage. Naively calling `ForwardDiff.derivative(pitch -> run_sim!(...).Fx, 0.0)` errors at `step_system!` because `system::DynamicSystem{Float64}` rejects duals — this is misleading; it's a usage error, not a WATT/GXBeam limitation. Plan.md Phase 2 deferred-bugs section should be updated next session to reflect this.
- AD verification approach: BEM unit-level uses `FiniteDifferences.central_fdm(5,1)` against a single scalar; coupled `run_sim!` uses a **directional derivative** check (random unit vector over compliance+mass region, scaled per-entry by `|p0|` because values span ~10 orders of magnitude). Per-entry spot checks were flaky due to symmetry of compliance/mass and dramatic scale variation; the directional approach integrates those out into one FD call. Tolerance kept wide (`rtol=1e-2`) — 3-step coupled FD is noisy.
- Test env deps added via `Pkg`: `Statistics`, `ForwardDiff`, `ReverseDiff`, `Random`. `FiniteDifferences` was tried and then removed by the user from WATT main deps (it lives only in test env).
- Cleanup: removed accidental `test/test/` directory I created by passing `--project=test` from inside `/test/` — interpreted as relative `./test/`.
- `fixedpoint!` is **not** broken (plan.md said it was). Converges to L2 relchange ~2e-7 by iter 10 on NREL 5MW rated; first iteration matches direct CCBlade to `rtol=1e-10`. Phase 1 likely fixed it incidentally.
- Final test counts: `Pkg.test()` = 1233 pass + 2 broken; `WATT_AD_TESTS=true Pkg.test()` = 1236 pass + 4 broken (those 4 reduce after the user's local cleanup of the test_ad broken markers — Phase 2 plan file `2-expressive-bunny.md` exit checklist all met).
