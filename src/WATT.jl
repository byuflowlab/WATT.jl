module WATT

using FLOWMath, LinearAlgebra, StaticArrays, NLsolve, UnPack
using CCBlade, GXBeam, DynamicStallModels
using ImplicitAD, ForwardDiff, ReverseDiff
using DelimitedFiles #Todo: Why is this here?
using RecipesBase
#Todo: Add the solution parameters as arguments to the user functions.

export BladePoints, AssemblyPlot

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
