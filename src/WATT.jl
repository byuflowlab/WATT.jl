module WATT

using FLOWMath, LinearAlgebra, StaticArrays, NLsolve, UnPack
using CCBlade, GXBeam, DynamicStallModels
using ImplicitAD, ForwardDiff, ReverseDiff
using DelimitedFiles #Needed for reading Turbsim files. 
using RecipesBase

# Types
export Rotor, Blade
export AeroStates, StaticAeroStates
export AbstractSimMesh, SimMesh, AeroMesh, StaticMesh, SurrogateMesh
export SurrogatePointState, SurrogateAssemblyState, AbstractStructuralSurrogate
export encode_initial, step_latent, decode, decode!
export RK4, BDF1

# Environment
export environment, SimpleEnvironment

# Simulation entry points
export initialize_aero, simulate!
export initialize_sim, run_sim!, run_sim
export step_solution!, initialize_from_state, run_from_state!
export initialize_sim_surrogate, run_sim_surrogate!, run_sim_surrogate
export initialize_static, fixedpoint!

# GPU backend
# Device-resident batched types and kernels. These mirror the CPU path above but
# are declared here, not in their own files, so the public surface stays in one
# place. Implementation details (DS_NSTATES, N_BRENT_ITERS_DEFAULT, the enums,
# the internal kernel helpers) stay unexported — reach them as `WATT.name`.
export BladeGPU, RotorGPU, GPUBEMTOutputs, solve_BEMT_gpu!
export DSAirfoilGPU, DSHistory, march_ds_gpu!, ds_init_step_gpu!, ds_step_gpu!
export AeroGeometryGPU, InterpGPU, GPUPointStates
export aero_velocities_gpu!, ds_inputs_gpu!, ds_loads_dimensionalize_gpu!
export build_surrogate_loads_gpu!, update_feedback_gpu!
export GPUSurrogateMesh, GPUSurrogateHistory
export initialize_sim_surrogate_gpu, run_sim_surrogate_gpu!

# Post-processing
export rotorloads

# Plotting (RecipesBase recipes — no hard Plots dependency)
export BladePoints, AssemblyPlot

DS = DynamicStallModels
IAD = ImplicitAD

### Structs, solvers, and base models. 
include("./solvers.jl")
include("./utils.jl")
include("./types.jl")
include("./environments.jl")
include("./mesh.jl")
include("./bemt.jl")
include("./bemt_gpu.jl")
include("./dynamicstallmodels.jl")
include("./dsmodel_gpu.jl")
include("./gxbeam.jl")


### Couplings
include("./aero_only.jl")
include("./static.jl")
include("./aerostructural.jl")
include("./aerostructural_surrogate.jl")
include("./aerostructural_surrogate_gpu.jl")

end # module
