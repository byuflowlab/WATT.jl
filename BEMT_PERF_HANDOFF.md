# WATT BEMT / airfoil performance handoff

**Author of request:** Adam Cardoza · **Branch:** `rewrite` · **Date:** 2026-07-31

You are working in the WATT dev package (`~/.julia/dev/WATT`). This doc asks for two
small, low-risk performance fixes to the surrogate-driven BEMT march, plus a third
follow-up experiment. All three came out of profiling the constraint-Jacobian of a
downstream optimization harness (`Cardoza2026_FastBladeOpt/old_optimization_improvement/
surrogate_objcon.jl`) — a colored ForwardDiff Jacobian over a 2001-step
`run_sim_surrogate!` march.

## Where the time goes (profiler result, for context)

Of the whole constraint Jacobian:

- **~78%** is the surrogate march (`run_sim_surrogate!`).
- Inside the march, per timestep: **`take_aero_step!` ≈ 37%** (full CCBlade BEMT +
  dynamic stall), Koopman `decode` ≈ 23%, `step_latent` ≈ 17%.
- Inside the aero step, **`solve_BEMT` ≈ 21% of the entire Jacobian**, split roughly:
  - `IAD.implicit` = 13.9% — of which **~12.1% is the primal Brent root-find** and only
    **~1% is the derivative machinery** (`jvp` + `drdy_forward` + linear solve).
  - `blade.airfoils[idx]` StructArray getindex = **3.6%** (+ related `createinstance`
    frames, ~8% cumulatively).
  - `firstbracket` (bracket search) = **~2.6%**.

**Key finding already established:** AD is *not* being propagated through Brent.
`IAD.implicit` (`ImplicitAD._implicit`, `~/.julia/dev/ImplicitAD/src/nonlinear.jl:82`)
strips duals with `xv = fd_value(x)`, runs `solve` (Brent) on Float64 primals, and
recovers the derivative via the implicit function theorem. That part is correct and
should NOT be changed. The primal Brent cost is genuine work. The two fixes below target
the *avoidable* overhead around it; the Brent cost itself is the follow-up experiment.

---

## Fix 1 — Store airfoils as a plain `Vector` in `Blade` (kill StructArray reconstruction)

### Problem
`Blade.airfoils` is currently populated (by callers) with a
`StructArray{DS.Airfoil}`. `DS.Airfoil` is an **immutable** struct whose fields are all
static (polar table, cl/cd/cm spline fits, curve slopes, separation function — see
`~/.julia/dev/DynamicStallModels/src/airfoils.jl`). But StructArray stores those fields
as **separate column arrays**, so every `airfoils[idx]` **reconstructs a fresh
`DS.Airfoil`** (gather each field + call the constructor). `solve_BEMT` /
`extract_ds_loads!` / the DS functor all index `airfoils[i]` in per-section,
per-timestep loops → ~17 × 2001 ≈ 34k reconstructions per march (× coloring passes under
the Jacobian). That is the 3.6% + `createinstance` overhead.

### Why it's safe
- `DS.Airfoil` is immutable and time-invariant. The dynamic-stall *state* lives in the
  separate `states/y_ds/p_ds` vectors, NOT in the airfoil struct
  (`dynamicstallmodels.jl:134` calls `airfoils(states_old, states_new, ...)`, passing
  state as arguments).
- Every airfoil use-site already works on `AbstractVector{<:DS.Airfoil}`:
  - `src/bemt.jl:77`, `src/bemt.jl:249`: `airfoil = blade.airfoils[idx]`
  - `src/dynamicstallmodels.jl` loops: `for j in eachindex(airfoils)` then `airfoils[j]`
  - the DS functor `(airfoils::AbstractVector{<:Airfoil})(...)`
    (`~/.julia/dev/DynamicStallModels/src/solve.jl:117,133`) — itself just loops and
    indexes `airfoils[i]`.
  None of them require StructArray column storage. A plain `Vector{DS.Airfoil{...}}`
  satisfies all of them, and `airfoils[i]` becomes an O(1) pointer load.
- The airfoils are design-independent Float64 (built from fixed polar data in the
  harness), so they are constant across the entire optimization — building the vector
  once is always valid.

### What to change (`src/types.jl`)
1. The `Blade` field is currently abstractly typed (a mild type instability too):
   ```julia
   struct Blade{TF, TF2}
       ...
       airfoils::AbstractVector{<:DS.Airfoil}   # abstract field
   end
   ```
   Parametrize on the airfoil-vector type so the field is concrete:
   ```julia
   struct Blade{TF, TF2, TAF<:AbstractVector{<:DS.Airfoil}}
       ...
       airfoils::TAF
   end
   ```
   (Update any `Blade{...}` type annotations elsewhere accordingly — check with
   `grep -rn "Blade{" src/`. Most code uses plain `::Blade`, which still works.)

2. In the convenience constructor `Blade(span, chord, twist, xcp, airfoils; ...)`
   (`src/types.jl:119`), **materialize once** so callers passing a StructArray still get
   a plain contiguous vector stored:
   ```julia
   airfoils = collect(airfoils)   # StructArray -> Vector{DS.Airfoil{...}} (concrete eltype), 1 reconstruction each
   ```
   Put this near the top of the constructor (after the length check). `collect` on a
   StructArray yields a `Vector` with the concrete element type. This alone hoists the
   34k reconstructions down to 17 (once per `Blade` construction).

### Optional (bigger win, harness-side, do only if you also touch the harness)
The `Blade` is rebuilt every constraint eval, so even with `collect` you pay 17
reconstructions per eval. If you additionally change the *harness* to build a
`Vector{DS.Airfoil}` once (instead of `StructArray{DS.Airfoil}(undef, nr)`) and reuse it,
the `collect` becomes a near-noop and reconstruction happens exactly once for the whole
optimization. This is out of scope for the WATT change but worth noting to the harness
owner. **Do not do the harness change from the WATT session** unless asked.

### Verify
- `grep -rn "StructArray" src/` — confirm nothing in WATT *depends* on the airfoils being
  a StructArray (it shouldn't; the harness constructs it).
- Run the existing WATT test suite / any `run_sim_surrogate!` smoke test. Loads must be
  **bitwise identical** before/after (immutable, same data) — assert exact equality, not
  just `isapprox`.

---

## Fix 2 — Bracket search on primals, not duals (`src/bemt.jl`)

### Problem
There are two `solve_BEMT` methods:
- `src/bemt.jl:75` (`...phi0, idx...; newbounds`) — `firstbracket` at **line 183**
- `src/bemt.jl:246` (`...idx...; epsilon`) — `firstbracket` at **line 355**

Both do:
```julia
success, phiL, phiU = CCBlade.firstbracket(phi -> residual(phi, xv, pv), phimin, phimax, npts, backwardsearch)
```
`xv` (the length-11 scratch, `mesh.xcc`) is **dual-typed under AD** — its eltype is
inferred from `blade.c` via `find_inittype` (`src/aero_only.jl:75,89`). `firstbracket`
does up to `npts` (=10) residual evaluations purely to find a **sign change** (a
bracket), but because `xv` is dual it computes full ForwardDiff partials on every probe
and **throws them away**. That's the ~2.6%.

The subsequent `IAD.implicit(solve, residual, xv, pv)` (lines 211 / 365) legitimately
needs the dual `xv` (for the IFT `jvp`), so we can't just strip `xv` globally — we strip
it **only for the bracketing closure**.

### What to change
In **both** methods, replace the `firstbracket` call so the residual it evaluates uses a
primal (value-stripped) copy of `xv`. `ForwardDiff.value` is identity on `Float64`, so
this is correct in both the primal and AD paths:

```julia
# once, before the quadrant loop (after update_BEMT_variables! fills xv):
xv_val = ForwardDiff.value.(xv)          # Float64 copy for sign-only bracketing
...
# inside the quadrant loop, bracket on primals:
success, phiL, phiU = CCBlade.firstbracket(phi -> residual(phi, xv_val, pv),
                                           phimin, phimax, npts, backwardsearch)
```
Leave `IAD.implicit(solve, residual, xv, pv)` and the final
`residual_and_outputs(phistar, xv, pv)` on the **dual** `xv` — those need the partials.

Notes:
- `ForwardDiff` is already a dep of WATT (used throughout); confirm it's imported in the
  module namespace used by `bemt.jl` (or qualify as `WATT.ForwardDiff` / whatever is in
  scope). `grep -n "ForwardDiff" src/bemt.jl src/WATT.jl`.
- `xv_val = ForwardDiff.value.(xv)` allocates one length-11 vector per `solve_BEMT` call.
  That's cheap relative to the ~10 dual residual evals it saves, but if you want zero
  alloc you can carry a second Float64 scratch of length 11 on the mesh (`mesh.xcc_val`)
  and `map!(ForwardDiff.value, xv_val, xv)`. Preallocation is optional — start with the
  simple version and measure.
- `phiL`, `phiU` come out as `Float64` (they're the primal bracket bounds), which is what
  Brent wants. `phistar` still gets its derivative from `IAD.implicit`. No change to
  output sensitivities.

### Verify
- Loads / gradients must match pre-change to tight tolerance (bracket bounds are found
  from primal residual signs, which equal the dual residual's value part — so the bracket
  is identical; `phistar` and all outputs unchanged).
- Cross-check one BEMT gradient via finite difference on a single section.

---

## Benchmark protocol (for Fixes 1 & 2)

The downstream harness already has a timing/profiling harness
(`Cardoza2026_FastBladeOpt/old_optimization_improvement/surrogate_objcon.jl`,
`RUN_FULL` / `RUN_JAC` / `RUN_PROFILE` flags). But from **within WATT** you can measure in
isolation with a single-march benchmark:

- Build a `Blade` + surrogate, run `run_sim_surrogate!` once to warm up, then
  `@benchmark` a single march (`BenchmarkTools`, `samples=10 evals=1`).
- For the AD path, wrap the march in a `ForwardDiff.jacobian` (or reuse the harness's
  colored Jacobian) so the dual-width cost is represented — the fixes only matter under
  duals.
- Report median before/after for: (a) primal march, (b) one dual march. Expected combined
  gain from Fixes 1+2 is **~6% of the constraint Jacobian** — real but modest. Confirm the
  numbers didn't change (exact for Fix 1, tight-tol for Fix 2).

---

## Follow-up (SEPARATE, after 1 & 2 land) — Warm-start Brent

The dominant remaining cost is the **primal Brent root-find running every section every
timestep** (~12% of the Jacobian, the biggest single genuine cost). Adam has tried
warm-starting before and doesn't recall a big gain, so approach this **systematically and
measure**, don't assume:

- **Idea:** the converged inflow angle `φ*` changes slowly between timesteps, so the
  previous step's `φ*` (per section) is a good initial bracket/guess. Options:
  1. Seed `firstbracket` with a narrow bracket around last step's `φ*` (fall back to the
     full quadrant search on failure). Fewer bracketing probes.
  2. Seed Brent itself with last step's `φ*` (Brent needs a bracket, not a point, so this
     mostly helps if paired with a tight bracket).
- **State plumbing:** you'd need to persist per-section `φ*` across timesteps
  (`mesh` already carries per-section history — `mesh.cchistory` is referenced in commented
  code at `src/aero_only.jl:21-23`; check whether a `φ` history field exists or add one).
- **Measure:** count average Brent iterations per solve before/after (instrument
  `FLOWMath.brent` or wrap it), not just wall time. If iteration count barely drops, the
  bracket/quadrant search is the cost, not Brent's convergence — which changes the fix.
- **Watch correctness:** warm-starting must not change the converged `φ*` (same root, same
  tolerance). The IFT derivative is unaffected (it depends only on the converged point),
  so gradients should be unchanged — assert this.

Report iteration-count and wall-time deltas so we can decide whether warm-starting is
worth keeping. If the gain is again marginal, document that and stop — the real lever may
be reducing BEMT call frequency (e.g. sub-stepping the aero less often than the structure)
or surrogating the aero loads, which are larger design changes for later.

---

## Scope guardrails
- Do **not** change the ImplicitAD path — the implicit-function-theorem handling is
  correct and is ~1% of cost.
- Do **not** commit/push unless explicitly asked.
- Keep Fix 1 and Fix 2 as separate commits so each can be benchmarked and reverted
  independently.
- Preserve exact numerics: Fix 1 is bitwise-identical; Fix 2 is identical up to the
  bracketing being done on the value part of the dual (which is exact).
