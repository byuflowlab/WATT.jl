Different from the static sweep. Under real unsteady coupled loading:

```
q1:      3650   (99.401%)
q3:        22   ( 0.599%)
q2/q4:      0
```
22 transient solves land in q3 — the (π/2, π) quadrant. Physically that's inflow angles between 90° and 180°, i.e. flow going the "wrong way" through a section (heavy stall, structural velocity dominating rotational, or a fluke of the turbulent gust interacting with a heavily deformed blade). The static sweep never hit these because it can't reproduce structural feedback or turbulent gusts.

Consequence for the current GPU kernel: those 22 (t, section) events would be flagged bracket-failed and get zero loads — different from what the CPU sim just did (it took the q3 fallback and returned real loads there).

Three options:

Accept the miss. 0.6% of transient solves get zeroed. Whether this matters depends on which sections/times get hit — a couple of tip-section hits during a gust are probably harmless; a bunch of hits at a load-carrying section during peak loading are not.
Extend the GPU kernel to also try q3 as a second sub-bracket scan + Brent. Doubles the residual evaluations per section (still cheap), but recovers those 22 solves.
First instrument to see where the q3 hits happen. I can extend the trace to record (t_step, section_idx, quadrant) tuples, then we look at which sections and time steps fall through. If they cluster at outer aero nodes during high-turbulence gusts, they may not matter for the surrogate use case.
I'd suggest (3) before deciding — 22 events is small enough to inspect one by one. Want me to extend the trace that way?

Also — should I lengthen tvec to something more like 30–60s to make sure we've sampled a wider slice of the turbulent inflow? The current 5s window is only ~7 rotor revolutions and might undercount q3.

Very clean picture. Every non-q1 hit is at the root cylinder region.

Summary
Where: All 53 q3 hits are at aero nodes j=1..6, all with r/R < 0.037. Zero hits anywhere else on the blade.

Phi values: The reported inflow angles are 90.1° – 93.8° — literally 0.1° to 4° past the q1/q3 boundary at π/2. These aren't wild solutions in the middle of q3; they're right at the fence.

When: Three clusters in the 30s window:

t ∈ [0, 1)s: 22 hits — the sim's startup transient before the structural response settles
t ≈ [15, 16)s: 19 hits — a turbulent gust in the TurbSim record
t ≈ [21, 22)s: 12 hits — a second gust
Physically why: at r=1.8m, Vy = Ω·r ≈ 0.27 m/s while Vx ≈ 11 m/s → geometric phi ≈ atan(11/0.27) ≈ 88.5°, sitting right at the q1 ceiling. A tiny perturbation from structural velocity or a turbulent v-component tips it past 90°. Cylinder sections also have cl ≈ 0 in the polar, so the residual has essentially no lift contribution to pull the root into q1.

Warp clustering on GPU: Sections j=1..6 are 6 consecutive thread-lanes, which fit inside a single warp (of 32) for any n_sims layout. So on a GPU launch, only 1 warp per timestep would ever need q3 — the other ~5–15 warps (for n_sims=20–50) never touch it. If we did option 2 conditional q3, the divergence penalty is tiny.

What I'd recommend now
Option 1 (accept the miss) is fine for this rotor.

Reason: those 6 sections are root cylinders. Their aerodynamic contribution to total thrust/torque is essentially nothing:

No lift (cylinders are cl ≈ 0 across all α)
Only cd × 0.5ρW²·c drag, at a chord that's small relative to mid-span
Moment arm r/R < 0.04, so their contribution to torque is negligible
Blade-root bending moment gets some contribution but it's dominated by mid-span loads
Setting them to zero loads (q1 fallback) is physically defensible — arguably more so than trusting the BEMT residual at cylinder sections in the first place, which is why WATT's coupled solver already tolerates CCBlade.Outputs() returns for un-solvable sections.

If you want belt-and-suspenders, option 2 with conditional q3 is cheap here (1 warp per step), and would be worth adding only if you later see a design where cylinder loads matter (e.g. blade-root fatigue calcs that care about that segment).