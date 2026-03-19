# WATT Cleanup Review

Go through each item and mark the **Decision** column with one of:
- `delete` — remove the code
- `keep` — intentionally retained, no action needed
- `fix` — needs to be corrected rather than deleted
- `move` — belongs somewhere else
- `defer` — revisit later

---

## aerostructural.jl

### 1. Gen 1 commented-out block — exports
**Lines:** 6–7
**What:** `#= export initialize, initial_condition, take_step!, simulate, simulate! =#`
**Why:** Replaced by `initialize_sim` / `run_sim!`. Already commented out pending deletion.
**Decision:** delete

---

### 2. Gen 1 commented-out block — `initialize()`
**Lines:** 55–162
**What:** Original initialize function, wrapped in `#= ... =#`
**Why:** Superseded by `initialize_sim()`. Different DS model call signature (no `p_ds`), different GXBeam state management. Already commented out pending deletion.
**Decision:** delete

---

### 3. Gen 1 commented-out block — `initial_condition!()`
**Lines:** 180–306
**What:** Original initial condition solver, wrapped in `#= ... =#`
**Why:** Uses `GXBeam.gxbeam_initial_conditions!()` (backward Euler) and `extract_ds_loads!` without `p_ds`. Replaced by inline logic in `run_sim!()`. Already commented out pending deletion.
**Decision:** delete

---

### 4. Gen 1 commented-out block — `take_step!()`
**Lines:** 308–430
**What:** Single-step function using `GXBeam.take_step()` (backward Euler on DAE), wrapped in `#= ... =#`
**Why:** Replaced by inline `GXBeam.step_system!()` logic (Newmark-Beta) in `run_sim!()`. Already commented out pending deletion.
**Decision:** delete

---

### 5. Gen 1 commented-out block — `simulate()` and `simulate!()`
**Lines:** 432–508
**What:** Both wrapped in `#= ... =#`. `simulate()` calls non-existent non-bang functions — broken. `simulate!()` calls the Gen 1 `take_step!` / `initial_condition!`.
**Why:** Both superseded by `run_sim!()`. Already commented out pending deletion.
**Decision:** delete

---

### 6. Empty stub — `run_sim()`
**Lines:** ~828–830
**What:** `function run_sim() end` — empty body, no implementation
**Why:** Appears to be a placeholder for a non-mutating wrapper around `run_sim!()` that was never written.
**Decision:** fix

---

### 7. `checkforwarnings()`
**Lines:** 14–36
**What:** Validates that rotor/blade inputs are in radians and within reasonable ranges
**Why:** Was used in Gen 1 `initialize()` (now dead). Commented out even there. Never called anywhere currently active.
**Decision:** keep

---

### 8. `find_inittype()`
**Lines:** 40–52
**What:** Determines the numeric type (Float64 vs ForwardDiff.Dual vs ReverseDiff.TrackedReal) from input variables
**Why:** Called by both the dead `initialize()` and the live `initialize_sim()`. Still needed.
**Decision:** keep

---

## bem.jl

### 9. Duplicate `solve_BEM!()` — first version (with `phi0` parameter)
**Lines:** 53–206
**What:** `solve_BEM!(rotor, blade, env, phi0, idx, Vx, Vy, pitch, xv; ...)`
**Why:** Has an extra `phi0` initial guess parameter not present in the second version. The second version (item 10) is called by the active code (`run_sim!` → `take_aero_step!`). Unclear if this version is called anywhere.
**Decision:** defer

---

### 10. Duplicate `solve_BEM!()` — second version (without `phi0`)
**Lines:** 235–370
**What:** `solve_BEM!(rotor, blade, env, idx, Vx, Vy, pitch, xv; ...)`
**Why:** Appears to be the actively used version. Companion to item 9.
**Decision:** keep

---

### 11. Non-bang wrapper `solve_BEM()`
**Lines:** 47–51
**What:** `solve_BEM(...)` — just calls `solve_BEM!(...)` and returns result
**Why:** Trivial passthrough. Exists for API symmetry but adds no value.
**Decision:** defer

---

### 12. `thrusttorque()` and `ccthrusttorque()`
**Lines:** 388–417
**What:** Two functions that integrate blade loads into rotor thrust/torque. Both TODO-commented as possibly duplicating CCBlade functionality.
**Why:** Neither is called anywhere in the package.
**Decision:** defer

---

### 13. `show_dual_vec()`
**Lines:** ~208–212
**What:** Debug helper that prints values of a ForwardDiff Dual array
**Why:** Only used for debugging; not called in active code paths.
**Decision:** keep

---

## environments.jl

### 14. `get_aerostructural_velocities()` — first overload (broken)
**Lines:** 205–218
**What:** `get_aerostructural_velocities(env, aerov, t, r, azimuth, precone, tilt, yaw, hubht)`
**Why:** Calls `get_aero_velocities()` with a signature that doesn't exist. Would error at runtime. The function itself has a TODO asking "What is this doing? Is it needed?" Never called anywhere.
**Decision:** defer

---

## gxbeam.jl

### 15. `plotpoints()` and `plotassembly()`
**Lines:** 240–323
**What:** Visualization helpers for GXBeam beam assemblies
**Why:** Never called anywhere in the package. Likely developed for manual debugging.
**Decision:** keep

---

### 16. `gxbeam_initial_conditions()` — non-bang version
**Lines:** 324–353
**What:** Non-mutating version that returns initial GXBeam state
**Why:** Never called anywhere. The bang version (`gxbeam_initial_conditions!`) was used only in the now-dead Gen 1 `initial_condition!()`.
**Decision:** delete

---

### 17. `gxbeam_initial_conditions!()` — bang version
**Lines:** 354–440
**What:** Mutating version of the above; sets up GXBeam initial conditions using the Gen 1 backward-Euler API
**Why:** Only referenced in the now-commented-out Gen 1 `initial_condition!()`. The active `run_sim!()` uses `GXBeam.initialize_system!()` directly instead.
**Decision:** delete

---

### 18. `simulate_gxbeam()`
**Lines:** 498–598
**What:** Standalone GXBeam-only transient simulation (no aero coupling)
**Why:** Never called in any active code path. Functionality subsumed by `run_sim!()`.
**Decision:** delete

---

### 19. `steady_simulate_gxbeam()`
**Lines:** 599–644
**What:** Standalone GXBeam steady-state simulation
**Why:** Never called in active code. Functionality subsumed by `static.jl`.
**Decision:** delete

---

### 20. `get_blade_weight()`
**Lines:** 645–657
**What:** Computes total blade weight from GXBeam assembly
**Why:** Never called anywhere in the package.
**Decision:** keep

---

## utils.jl

### 21. `derivative_me()`, `mat_derivative()`
**Lines:** 64–87
**What:** Manual finite-difference derivative utilities
**Why:** Never called. Likely predates AD integration.
**Decision:** defer

---

### 22. `prepareextra()`, `saveextra()`, `readextra()`
**Lines:** 89–134
**What:** File I/O helpers for saving/loading extra simulation state
**Why:** Never called. Likely from an early development phase before the current data storage approach.
**Decision:** delete

---

### 23. `plotdshistory()`
**Lines:** 136–149
**What:** Plots dynamic stall state history
**Why:** Never called. Likely a debug visualization from early development.
**Decision:** delete

---

### 24. `sub_brent()`
**Lines:** 279–343
**What:** Custom Brent's method root-finder
**Why:** Commented out in bem.jl (line 352). Superseded by FLOWMath's Brent solver. Never called.
**Decision:** delete

---

## aero_only.jl

### 25. `rotorloads()` — both overloads
**Lines:** 342–383
**What:** Two overloads that integrate sectional aerodynamic loads into rotor thrust/torque
**Why:** Neither overload is called anywhere internally. May be user-facing API, but not tested.
**Decision:** keep

---

## Legacy files — beddoesleishman*.jl

### 26. `beddoesleishman.jl`
**What:** Local implementation of the original Beddoes-Leishman dynamic stall model
**Why:** Not `include()`d in WATT.jl. Package now uses the `DynamicStallModels.jl` external package instead.
**Decision:** delete

---

### 27. `beddoesleishman_aerodyn.jl`
**What:** AeroDyn variant of Beddoes-Leishman dynamic stall
**Why:** Not `include()`d in WATT.jl. Same situation as item 26.
**Decision:** delete

---

### 28. `beddoesleishman_gonzalez.jl`
**What:** Gonzalez variant of Beddoes-Leishman dynamic stall
**Why:** Not `include()`d in WATT.jl. Same situation as item 26.
**Decision:** delete

---

## Naming / namespace issues (not deletion, but worth noting)

### 29. `initialize()` defined in three places
**Files:** `static.jl` (line 14), `aero_only.jl` (line 49), `aerostructural.jl` (dead Gen 1 code)
**What:** All three files define an `initialize()` function for their respective simulation modes
**Why:** Julia's multiple dispatch handles this if signatures differ, but the names are ambiguous to users. Consider `initialize_static()`, `initialize_aero()`, etc.
**Decision:** defer

---

## Possible signature mismatch (needs verification)

### 30. `initialize_ds_model()` return count
**File:** `dynamicstallmodels.jl` lines 3–24 vs `aerostructural.jl` line 516
**What:** `initialize_sim()` unpacks 4 return values `(xds, xds_idxs, y_ds, p_ds)` but the Gen 1 code unpacked only 3. The actual definition in `dynamicstallmodels.jl` should be verified to confirm it returns 4 values.
**Why:** If the function only returns 3, `run_sim!()` would error on initialization.
**Decision:** defer
