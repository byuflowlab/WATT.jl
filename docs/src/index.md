# WATT.jl Documentation - Wind Aeroelastic Turbine Toolkit

A toolkit for nonlinear unsteady aeroelastic modeling of wind turbine blades, specifically designed for derivative computation. Our model couples blade element momentum theory ([CCBlade](https://github.com/byuflowlab/CCBlade.jl.git)), a dynamic stall model ([DynamicStallModels](https://github.com/byuflowlab/DynamicStallModels.jl.git)), and geometrically exact beam theory ([GXBeam](https://github.com/byuflowlab/GXBeam.jl.git)). 

## Installation
```julia
pkg> add https://github.com/byuflowlab/WATT.jl.git
```

## Documentation
- To start, check out the quick start tutorial. 
- Then take a look at the examples.
- If more information is required, checkout the API in the reference. 


## Capabilities
- Nonlinear unsteady aerostructural analysis
- AD compatibility
- Nonlinear steady aerostructural analysis


## Overview
The model is based on [OpenFAST](https://github.com/OpenFAST/openfast) from the National Renewable Energy Laboratory. This implementation is not intended to replace OpenFAST, but rather serves as a research platform for rapidly prototyping and evaluating differentiation techniques.

Key differences:
- WATT.jl is compatible with mature algorithmic differentiation packages including ForwardDiff, ReverseDiff, and ImplicitAD. 
- CCBlade uses a slightly different implementation of Brent's method which makes small differences which accumulate overtime. The implementation that CCBlade uses converges to tighter tolerances. 
- GXBeam uses constant property linear elements with extended Milenković parameters. This formulation produces an exceptionally robust solver that avoids excessive quadrature, enabling tight convergence of structural states at each time step.



## Developers
I suggest you first work through the quick start tutorial and examples, then checkout the Developers section for a outline of the code and basic theory description. 