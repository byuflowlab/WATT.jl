module WATT

using FLOWMath, LinearAlgebra, StaticArrays, NLsolve, UnPack
using CCBlade, GXBeam, DynamicStallModels
using ImplicitAD, ForwardDiff, ReverseDiff
using DelimitedFiles #Todo: Why is this here?
using RecipesBase
#Todo: Add the solution parameters as arguments to the user functions.

# Types
export Rotor, Blade
export AeroStates, StaticAeroStates
export AbstractSimMesh, SimMesh, AeroMesh, StaticMesh
export RK4, BDF1

# Environment
export environment, SimpleEnvironment

# Simulation entry points
export initialize_aero, simulate!
export initialize_sim, run_sim!, run_sim
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
include("./bem.jl")
include("./dynamicstallmodels.jl")
include("./gxbeam.jl")


### Couplings
include("./aero_only.jl")
include("./aerostructural.jl")

# include("../dev/indicialgxbeam.jl")

include("./static.jl")  


end # module
