# Rotors.jl Documentation

## Summary

Rotors.jl is an all Julia multi-physics wind turbine model patterned after OpenFAST. Similar to OpenFAST it uses a BEM (implemented here by the package CCBlade.jl) with a dynamic stall model (implemented here by the package DynamicStallModels.jl ) to model aerodynamics. Structural dynamics are modeled by the package GXBeam.jl, an implementation of geometrically exact beam theory. 

Each of these packages accomplish similar feats to OpenFAST's rotor analysis package with these key exceptions:
- CCBlade.jl allows for analysis of propellers, and helicopter blades in addition to wind turbine blades. 
- More dynamic stall models are available through dynamicstallmodels.jl.
- GXBeam uses constant and linear shape elements instead of Legendre-spectral finite elements (sic).
- All of the code, and packages used are in julia, allowing for easy implementation of automatic derivatives (AD). 

Since all of the code can use AD, this package will be improved to provide derivatives of design variables with respect to outputs, throughout a time domain simulation. 

Some function to easily parse OpenFAST files can be found through the package OpenFASTsr.jl, which allows for easy creation of inputs for this package. 


## Installation
```julia
pkg> add https://github.com/byuflowlab/Rotors.jl.git
```

## Documentation
- To start, check out the quick start tutorial. 
- Then take a look at the examples.
- If more information is required, checkout the API in the reference. 



## Developers
I suggest you first work through the quick start tutorial and examples, then checkout the Developers section for a outline of the code and basic theory description. 