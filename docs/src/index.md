# WATT.jl Documentation - Wind Aeroelastic Turbine Toolkit

A toolkit for nonlinear unsteady aeroelastic modeling of wind turbine blades, specifically designed for derivative computation. Our model couples blade element momentum theory ([CCBlade](https://github.com/byuflowlab/CCBlade.jl.git)), a dynamic stall model ([DynamicStallModels](https://github.com/byuflowlab/DynamicStallModels.jl.git)), and geometrically exact beam theory ([GXBeam](https://github.com/byuflowlab/GXBeam.jl.git)). 

## Installation
```julia
pkg> add https://github.com/byuflowlab/WATT.jl.git
```

## Documentation
- To start, see [Getting Started](gettingstarted.md).
- [For Developers](developers.md) covers the code layout, reference frames, coupling, and AD architecture.
- [API Reference](apireference.md) documents every exported symbol; the [GPU Backend](gpu.md) has its own page.


## Capabilities
- Nonlinear unsteady aerostructural analysis
- Nonlinear steady aerostructural analysis
- Aerodynamics-only transient analysis
- AD compatibility — ForwardDiff and ReverseDiff, end to end
- Windowed sensitivities from a frozen starting state, for surrogate training
- Optional learned structural surrogate in place of the beam solve
- Batched GPU aerodynamics for evaluating many simulations at once (forward pass only)


## Overview
The model is based on [OpenFAST](https://github.com/OpenFAST/openfast) from the National Laboratory of the Rockies (previously the National Renewable Energy Laboratory). This implementation is not intended to replace OpenFAST, but rather serves as a research platform for rapidly prototyping and evaluating differentiation techniques.

Key differences:
- WATT.jl is compatible with mature algorithmic differentiation packages including ForwardDiff, ReverseDiff, and ImplicitAD. 
- CCBlade uses a slightly different implementation of Brent's method which makes small differences which accumulate overtime. The implementation that CCBlade uses converges to tighter tolerances. 
- GXBeam uses constant property linear elements with extended Milenković parameters. This formulation produces an exceptionally robust solver that avoids excessive quadrature, enabling tight convergence of structural states at each time step.
