module Rotors

using DifferentialEquations, FLOWMath, CCBlade, GXBeam, LinearAlgebra, StaticArrays, Plots, CurveFit, NLsolve, DynamicStallModels

DS = DynamicStallModels
DE = DifferentialEquations

### Structs, solvers, and whatnot. 
include("./types.jl")
include("./solvers.jl")
include("./utils.jl")


### Base model
include("./blades.jl")
include("./environments.jl")


### Dynamic stall models. 
include("./riso.jl") 
include("./beddoesleishman_aerodyn.jl")
include("./beddoesleishman_gonzalez.jl")
include("./beddoesleishman.jl")



### Structures
# include("../dev/gxbeam.jl") #Todo: Need to get the relavent functions out. Then remove the file. -> Or just simplify the file. 
#Todo: I need to unify the structural and aerodynamic reference frames. 
include("./gxbeam.jl")


### Couplings
include("./aero_only.jl")
include("./aerostructural.jl")

# include("../dev/indicialgxbeam.jl")

include("./static.jl") #Todo: Need to update to not use BEM struct and Riso struct and the like. 


end # module
