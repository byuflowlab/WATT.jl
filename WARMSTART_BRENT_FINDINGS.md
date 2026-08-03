# Warm-start Brent in BEMT — findings (spike, do-not-pursue)

**Author:** Adam Cardoza (with Claude) · **Repo:** `~/.julia/dev/WATT` · **Date:** 2026-08-01
**Status:** investigated, measured, **reverted** — warm-start is not worth pursuing here.

## TL;DR

The dominant *remaining* cost in the surrogate-driven BEMT march (after the Vector-airfoils
and primal-bracketing fixes landed) is the **primal Brent root-find run every section every
timestep** — ~12% of the constraint Jacobian in the earlier profile. The natural idea is to
warm-start the root-find from the previous timestep's converged inflow angle φ*, which moves
slowly step-to-step. We built it behind an opt-in flag, instrumented the root-find, and
measured. **Conclusion: the gain is within measurement noise (~1–3% of the aero march), and a
strictly-correct version is defeated by poles in CCBlade's residual. Not worth it.** This
matches Adam's prior recollection that warm-starting "didn't help much." All spike code was
reverted; only the instrumentation approach is documented below for reuse.

## What we measured

Instrumented `solve_BEMT` (method B, the one every march uses) with an opt-in counter (same
zero-overhead `Ref{Union{Nothing,Dict}}` model as the existing `BEMT_QUAD_TRACE`), capturing
per solve: Brent iterations and function-calls (from `FLOWMath.brent`'s returned
`info.iter`/`info.fcalls`, normally discarded via `phistar, _ = ...`) and `CCBlade.firstbracket`
residual probes (by wrapping the residual closure with a probe counter).

**Workload:** aero-only transient march on the golden NREL-5MW fixture
(`test/simpleNREL5_100s.jld2`), 25 sections × 400 steps = **10,000 solves**. Warm-start seed
was the previous step's converged φ* for that section (`aerostates.phi[i-1, j]`, which is
already stored — it just was never read back as a seed).

### Baseline root-find cost (`:off`, per solve)

| Metric | Value |
|---|---|
| Brent iterations | **7.25** |
| Brent function-calls | **8.25** |
| `firstbracket` probes | **3.11** (worst case is `npts=10`) |
| short-circuit hits | **2 / 10,000** |

**The bracket search is already cheap** (~3 probes, not 10) — the residual's sign change is
found in the first couple of probes almost every time. **Brent convergence dominates** the
residual-evaluation budget (~8 fcalls vs ~3 probes). This is the crux: warm-start can only
cheaply skip the *bracket search*, and the bracket search isn't where the time goes.

## Two warm-start modes tried

1. **`:bracket` (strict, intended to be exact):** seed Brent with a narrow bracket
   `[φ*−δ, φ*+δ]` (δ = 5°) when it straddles the root; Brent still converges to full tolerance
   (`atol=2e-12`), so — *in principle* — the converged φ* and all IFT-derived gradients are
   unchanged. Fall back to the full quadrant search on a miss.
   - Hit rate was high (93–98% of solves), and it did drop probes 3.11→2.1 and Brent iters
     7.25→7.06. But it was **NOT numerically identical** — the march's AD derivative came out
     0.5–11% wrong, with per-element load errors up to ~60×.
   - **Root cause (important for anyone reusing CCBlade):** the BEMT residual has **interior
     poles**, not just the obvious ones at φ=0 and ±π/2. Where the axial induction factor
     a→1 (the momentum / empirical-region boundary), CCBlade's residual flips sign
     discontinuously. A narrow window around φ* can straddle such a pole and bracket a
     **spurious** sign change; Brent then converges to a non-physical "root." Clamping the
     window to φ*'s quadrant does **not** fix it — the pole is *inside* the quadrant. A
     correctness-safe narrow-bracket warm-start would need to detect the physical region
     (a > a_crit), which is more machinery than the payoff justifies.

2. **`:shortcircuit` (aggressive):** if `|residual(φ*)| < 1e-6`, return the loads at φ* without
   bracketing or Brent (this is the logic the dead `phi0` method of `solve_BEMT` already had).
   - Fired **2 times in 10,000**. φ* moves enough between steps that its previous value is
     essentially never a 1e-6 root, so this is a no-op — and it actually costs one extra
     residual probe per solve, making it marginally *slower*. Perturbs loads ~1e-6 when it
     does fire.

## Wall time (median ms, aero-only march, tracing off)

| Mode | Primal | Dual (ForwardDiff) |
|---|---|---|
| `:off` | 196–203 | 231–236 |
| `:bracket` | 193–195 | 223–227 |
| `:shortcircuit` | 204–209 | 231–234 |

`:bracket` looks ~2–4% faster, but `:off` alone varied 196→203 ms run-to-run — the difference
is inside the noise, and the `:bracket` numbers are *tainted* anyway because it converges to
some wrong roots (possibly faster precisely because they're spurious). `:shortcircuit` is not
faster. Net: **no reliable, correct speedup.**

## Why the ceiling is low (independent of the bugs)

Even a *perfect* warm-start can only touch the bracket-search portion of `solve_BEMT`. The
earlier profile put the bracket search at ~2.6% of the whole Jacobian, and this spike confirms
it's already near-minimal (~3 probes/solve). Brent's ~7 iterations are set by its
tolerance and the residual's local shape, not by the initial bracket width — a tighter bracket
saved <0.2 iterations. So the achievable win is a fraction of ~2.6% of the Jacobian: real work
for a change that is, at best, in the noise.

## Recommendation

**Do not pursue warm-starting the BEMT root-find.** The real remaining levers are the larger
design changes the original performance handoff flagged, which reduce *how often* the BEMT is
solved rather than shaving each solve:

- **Sub-step the aero less often than the structure** — hold BEMT loads across several
  structural steps (the loads change slowly relative to the structural integrator's steps).
- **Surrogate the aero loads** — replace the per-section BEMT + DS solve with a learned model
  (the direction the Koopman/surrogate work is already exploring).

Both are proper design changes, not micro-optimizations, and each removes far more than the
~2.6% ceiling warm-start is bounded by.

## Reusable takeaways for the other project

- `FLOWMath.brent` returns `(xstar, info)` with `info.iter` / `info.fcalls` — cheap to
  instrument iteration counts; the call sites currently throw `info` away.
- `CCBlade.firstbracket` does up to `npts` probes and returns early on the first sign change;
  it does **not** report probe count — wrap the closure to count.
- **CCBlade's BEMT residual has interior sign-flips at the a→1 boundary.** Any scheme that
  hands Brent a bracket not produced by the quadrant `firstbracket` scan (warm-start, custom
  seeding, sub-domain solves) must guard against bracketing across that discontinuity, or it
  will converge to spurious roots. This is the single most important gotcha here.
- If you *do* want a correct warm-start elsewhere, the safe form is: keep the full quadrant
  `firstbracket` to establish the *true* bracket, and only use φ* to seed an interior guess to
  a solver that accepts one (Brent does not) — but given the measured ~3-probe bracket cost,
  the upside is minimal.
