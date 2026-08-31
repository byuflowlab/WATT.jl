# GPU Backend

Device-resident, batched implementations of the aerodynamic stack, built on
[KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl). They
exist to evaluate **many simulations at once** — a population of designs, a
sweep of inflow conditions — rather than to make a single simulation faster.

```@meta
CurrentModule = WATT
```

## Status and limitations

Read this before building on the GPU path:

- **Forward pass only.** Automatic differentiation does not propagate through
  any GPU kernel. If you need gradients, use the CPU path — see
  [Windowed sensitivities](apireference.md#Windowed-sensitivities).
- **The structural model is the surrogate**, not GXBeam. The coupled GPU march
  pairs GPU aerodynamics with the learned structural surrogate; there is no
  device-side geometrically exact beam solve.
- **Validation lives in `examples/`, not in the test suite.** The BEMT kernel was
  checked against the CPU path on an H200 node, and the dynamic stall port
  against the CPU backend, but nothing in CI guards these.
- **Backend-generic.** The kernels dispatch through a generic `similar_type`, so
  no per-backend glue is required — the same code runs on a CPU backend, which
  is how the ports were validated.

## Batched types

Each mirrors a CPU type, but stores its fields as device arrays with a batch
dimension so one kernel launch covers every section, blade, and simulation.

```@docs
BladeGPU
RotorGPU
DSAirfoilGPU
AeroGeometryGPU
InterpGPU
GPUBEMTOutputs
GPUPointStates
DSHistory
```

## Kernels

### Blade element momentum theory

The inflow-angle root-find cannot branch per thread on a convergence test
without wrecking warp coherence, so it runs a fixed iteration count for every
section instead — `WATT.N_BRENT_ITERS_DEFAULT`.

```@docs
solve_BEMT_gpu!
```

### Dynamic stall

```@docs
march_ds_gpu!
ds_init_step_gpu!
ds_step_gpu!
```

### Coupling

```@docs
aero_velocities_gpu!
ds_inputs_gpu!
ds_loads_dimensionalize_gpu!
build_surrogate_loads_gpu!
update_feedback_gpu!
```

## Coupled march

The GPU analogue of [`initialize_sim_surrogate`](@ref) and
[`run_sim_surrogate!`](@ref).

```@docs
GPUSurrogateMesh
GPUSurrogateHistory
initialize_sim_surrogate_gpu
run_sim_surrogate_gpu!
```

## Examples

Runnable scripts, including the cluster benchmarks these kernels were validated
with, are in the `examples/` directory:

| Script | What it does |
|---|---|
| `gpu_bemt_benchmark.jl` | BEMT kernel throughput |
| `gpu_bemt_cpu_backend_check.jl` | BEMT kernel vs. the CPU path |
| `gpu_dsmodel_validate.jl` | Dynamic stall port vs. the CPU reference |
| `gpu_dsmodel_benchmark.jl` | Dynamic stall throughput |
| `gpu_aerostructural_surrogate.jl` | Full coupled GPU march |
| `gpu_aerostructural_benchmark.jl` | Three-way CPU / surrogate / GPU comparison |
