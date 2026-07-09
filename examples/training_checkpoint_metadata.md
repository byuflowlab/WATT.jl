# Add normalization & convention metadata to saved checkpoints

## Goal

The trained Neural Koopman checkpoints (`JLD2.jldsave(...)` calls at the end of
the curriculum loop *and* at the final model save) currently omit the
preprocessing constants the inference side needs. That forces inference code
(e.g. `WATT/examples/aerostructural_nrel5mw5seg_surrogate.jl`) to maintain a
parallel copy of `NORM_SCALES`, `FNORM`, `N_REGRID_NODES`, and
`ALL_STATE_TYPES` — which is exactly the source of every "scaling mismatch"
debugging session we've had.

Augment **both** `jldsave` calls in the training files (the per-phase
`checkpoint_*` save inside the curriculum loop, and the final `model_*` save
after the loop) to include every preprocessing constant the inference side
needs to fully reconstruct an input or invert a prediction.

## Fields to add to every `jldsave` call

Add the following keyword args to each `JLD2.jldsave(...)` invocation
(both the per-phase checkpoint and the final model save). Put them next to
the existing `s_elem`, `x_mean`, `x_std`, `state_indices` entries:

```julia
# --- state preprocessing ---
norm_scales_vec      = collect(scale_vec),                        # length 18 (Float64), per-component scale, ALREADY in the order matching state_indices
norm_scales_dict     = NORM_SCALES,                               # full per-type dict (Dict{String, Vector{Float64}}); redundant with norm_scales_vec but useful for human inspection
all_state_types      = collect(ALL_STATE_TYPES),                  # ("u","udot","theta","thetadot","V","Vdot","Omega","Omegadot","F","M")
n_regrid_nodes       = N_REGRID_NODES,                            # 6 — radial regrid grid size

# --- force preprocessing ---
fnorm                = fnorm,                                     # 1e4 — force scale divisor applied before regrid
```

### Field reference

| field | type | meaning |
|---|---|---|
| `norm_scales_vec` | `Vector{Float64}` length 18 | Per-component scale, in the same order as `state_indices`. The inference side multiplies inputs by this before the encoder and divides outputs by this after the decoder. |
| `norm_scales_dict` | `Dict{String, Vector{Float64}}` | Human-readable version of `NORM_SCALES` — maps each state-type name to its xyz scale triple. Inference can ignore this; it's stored so a person inspecting the JLD2 can see the convention. |
| `all_state_types` | `Vector{String}` length 10 | The canonical state-type names corresponding to GXBeam's 30-component point state, in groups of 3. Lets inference decode `state_indices` into human-readable labels. |
| `n_regrid_nodes` | `Int` | The number of radial nodes the state/force are regridded to before entering the network (currently 6). |
| `fnorm` | `Float64` | Force-side normalization divisor applied before regridding (currently 1e4). |

## Inference-side contract (informational — no code change needed here)

Inference code (e.g. `aerostructural_nrel5mw5seg_surrogate.jl`) will then load
those fields verbatim from the checkpoint and stop hard-coding them. The
parallel constants currently defined in the inference file —
`NORM_SCALES_VEC`, `FNORM`, `N_REGRID`, and `N_SEL` — will be derived from
the checkpoint:

```julia
ckpt = JLD2.load(model_path)
norm_scales_vec = ckpt["norm_scales_vec"]
fnorm           = ckpt["fnorm"]
n_regrid        = ckpt["n_regrid_nodes"]
n_sel           = length(ckpt["state_indices"])
```

## Where to make the edits

In each of the training files we have on disk
(`examples/x_jordan_lux_f64_seed.jl` and any sibling training scripts):

1. **Per-phase checkpoint save inside the curriculum loop.** Look for the
   `JLD2.jldsave("checkpoint_*.jld2"; ...)` call inside the
   `for n_steps in TIME_STEPS_SCHEDULE` loop.
2. **Final model save after the loop.** Look for the
   `JLD2.jldsave("model_*.jld2"; ...)` call near the bottom of the file
   (right after the prediction plot block).

Add the same set of fields to both. No other logic in the training script
needs to change — the constants being saved are already defined globally
(`scale_vec`, `NORM_SCALES`, `ALL_STATE_TYPES`, `N_REGRID_NODES`, `fnorm`).

## Validation

After the change, running the training file end-to-end should produce a JLD2
that, when loaded, reports the new fields:

```julia
julia> JLD2.load(model_path)["norm_scales_vec"]
18-element Vector{Float64}: …
julia> JLD2.load(model_path)["fnorm"]
10000.0
```

The inference example will be updated separately to consume these and
delete the hard-coded duplicates.
