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
