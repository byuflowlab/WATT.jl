module WATT

using FLOWMath, LinearAlgebra, StaticArrays, CurveFit, NLsolve, UnPack #Todo: Am I really using CurveFit? 
using CCBlade, GXBeam, DynamicStallModels
using ImplicitAD, ForwardDiff, ReverseDiff
using DelimitedFiles #Todo: Why is this here? 
using Plots 
#Todo: Why do I have plots as a dependency? -> So I can plot mid simulation. I wonder if there is a way to pass plots in or something. 
#Todo: Why do I have NLsolve as a dependency? -> I think for the BDF solvers... but I don't know if I need that. 
#Todo: Add the solution parameters as arguments to the user functions. 

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
