module Rotors

using DifferentialEquations, FLOWMath, CCBlade, GXBeam, LinearAlgebra, StaticArrays, Plots, CurveFit, NLsolve, OpenFASTsr, dynamicstallmodels

#Todo: Do I need OpenFASTsr? 

DS = dynamicstallmodels
DE = DifferentialEquations

### Structs, solvers, and whatnot. 
include("./types.jl")
include("./solvers.jl")
include("./utils.jl")


### Base model
include("./blades.jl")
include("./environments.jl")


### Dynamic stall models. 
include("./riso.jl") #Todo: Need to move relavent functions to this file, and main body functions to dynamicstallmodels.jl


### Structures
include("../dev/gxbeam.jl") #Todo: Need to get the relavent functions out. Then remove the file. 


### Couplings
include("./loosely.jl")
include("./coupled.jl")

include("./static.jl") #Todo: Need to update to not use BEM struct and Riso struct and the like. 


end # module
