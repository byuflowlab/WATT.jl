# WATT.jl (Wind Aeroelastic Turbine Toolkit)
A toolkit for nonlinear unsteady aeroelastic modeling of wind turbine blades, specifically designed for derivative computation. Our model couples blade element momentum theory ([CCBlade](https://github.com/byuflowlab/CCBlade.jl.git)), a dynamic stall model ([DynamicStallModels](https://github.com/byuflowlab/DynamicStallModels.jl.git)), and geometrically exact beam theory ([GXBeam](https://github.com/byuflowlab/GXBeam.jl.git)). 

### Capabilities
- Nonlinear unsteady aerostructural analysis
- AD compatibility
- Nonlinear steady aerostructural analysis


### Overview
The model is based on [OpenFAST](https://github.com/OpenFAST/openfast) from the National Renewable Energy Laboratory. This implementation is not intended to replace OpenFAST, but rather serves as a research platform for rapidly prototyping and evaluating differentiation techniques.

Key differences:
- WATT.jl is compatible with mature algorithmic differentiation packages including ForwardDiff, ReverseDiff, and ImplicitAD. 
- CCBlade uses a slightly different implementation of Brent's method which makes small differences which accumulate overtime. The implementation that CCBlade uses converges to tighter tolerances. 
- GXBeam uses constant property linear elements with extended Milenković parameters. This formulation produces an exceptionally robust solver that avoids excessive quadrature, enabling tight convergence of structural states at each time step.





<!-- ### Aeroelastic wind turbine primer
The core of the model is a two-way coupled, partitioned model that converges across time. It uses blade element moment theory (BEMT) to predict an angle of attack and an inflow velocity. The model uses the angle of attack and inflow angle to predict time-varying loads via the Beddoes-Leishman dynamic stall model. Geometrically exact beam theory (GEBT) uses those loads to predict the dynamic response of the wind turbine blade. The linear and angular deflections and velocities are used to augment the inflow conditions seen by the BEMT for the next time step. -->


