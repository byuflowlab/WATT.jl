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
export encode_initial, step_latent, decode
export RK4, BDF1

# Environment
export environment, SimpleEnvironment

# Simulation entry points
export initialize_aero, simulate!
export initialize_sim, run_sim!, run_sim
export initialize_sim_surrogate, run_sim_surrogate!, run_sim_surrogate
export initialize_static, fixedpoint!

# Post-processing
export rotorloads

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

end # module
