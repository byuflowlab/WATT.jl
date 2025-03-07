module Rotors

using FLOWMath, LinearAlgebra, StaticArrays, CurveFit, NLsolve, UnPack
using CCBlade, GXBeam, DynamicStallModels
using ImplicitAD, ForwardDiff, ReverseDiff
using Infiltrator
using DelimitedFiles
using Plots
#Todo: Why do I have plots as a dependency? -> So I can plot mid simulation. I wonder if there is a way to pass plots in or something. 
#Todo: Why do I have NLsolve as a dependency? -> I think for the BDF solvers... but I don't know if I need that. 
#Todo: Add the solution parameters as arguments to the user functions. 

DS = DynamicStallModels
IAD = ImplicitAD
# using DifferentialEquations
# DE = DifferentialEquations #Todo. I'm not sure that I need this as a dependency

### Structs, solvers, and base models. 
include("./solvers.jl")
include("./utils.jl")
include("./types.jl")
include("./environments.jl")
include("./mesh.jl")
include("./bem.jl")
include("./dynamicstallmodels.jl")
include("./gxbeam.jl")

### Dynamic Stall Model
# include("./riso.jl") 
# include("./beddoesleishman_aerodyn.jl")
# include("./beddoesleishman_gonzalez.jl")
# include("./beddoesleishman.jl")



### Structures
# include("../dev/gxbeam.jl") #Todo: Need to get the relavent functions out. Then remove the file. -> Or just simplify the file. 
#Todo: I need to unify the structural and aerodynamic reference frames. 
# include("./gxbeam.jl")


### Couplings
include("./aero_only.jl")
include("./aerostructural.jl")

# include("../dev/indicialgxbeam.jl")

include("./static.jl")  


end # module
