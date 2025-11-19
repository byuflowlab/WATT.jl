# For Developers
Here we describe how the code is organized (for easier debugging and additions to the code) and how the models are coupled together. 

## Code Outline



## Current Coupling (Loosely Coupled)
As is done in OpenFAST, this package is a two-way partitioned, loose coupling of explicit aerodynamic and structural models. The aerodynamic models are a blade element momentum theory (BEM) and the Beddoes-Leishman dynamic stall model. Blade element momentum theory is a common method that is well documented and can be found in most common aerodynamics textbooks. The method assumes steady inflow and relies on the integral forms of the linear and angular momentum equations. Ning presents a variant of the BEM that uses a one-dimensional residual equation which allows for guaranteed convergence. 

Environmental flow conditions and static airfoil data is used to calculate the inflow angle of the rotor (by solving the residual formed by the blade element and momentum equations as done by CCBlade.). That inflow angle is used to calculate the angle of attack of the airfoil. 

Since the BEM assumes steady inflow, we couple the BEM with various dynamic stall models, currently the most developed is the Beddoes-Leishman dynamic stall model with Gonzalez’s modifications. In this implementation of the Beddoes-Leishman model, the new states of the dynamic stall model are explicitly calculated based on the inflow conditions (angle of attack and inflow velocity) and current states.

Instead of looking up the lift and drag to calculate the loadings, that angle of attack is passed to the dynamic stall model. The states of the dynamic stall model are updated (based on the angle of attack and the environmental conditions). Based on the dynamic states, the lift and drag are used to calculate the loading on the rotor. Those loadings are fed to the structural model to simulate the structural response for a given time step. This process is repeated throughout time. 

For the structural model, we use geometrically exact beam theory (GEBT). This model is especially useful for computationally efficient modeling of beams with arbitrary section topology and materials, while exactly capturing all the geometric non-linearities. We use GXBeam.jl, a Julia implementation of GEBT, to calculate the dynamic response of the wind turbine blade. GXBeam also has finite element analysis capabilities to model the cross-sectional structural parameters (mass and stiffness matrices). It is important to note that GXBeam uses constant and linear shape functions, which increases the number of elements required to achieve grid-independent results. We use the integrator provided with GXBeam, which is an undamped Newmark scheme, to integrate the states across time.

The phrase two-way coupling means that both the aerodynamic and structural components are dependent on outputs from the other model. The inputs to the aerodynamic model from the structural model are the linear and angular deflections and linear velocities at every aerodynamic node. The inputs to the structural model from the aerodynamic model are the distributed loads, including the normal forces, tangential forces, and the axial moment. Since the aerodynamic and structural nodes are not required to be co-located, the deflections are linearly interpolated from the structural to aerodynamic nodes. The distributed loads are integrated onto the structural elements via Gaussian integration. The aerodynamic model is initialized to steady state before initializing the structural model. The structural model is initiated from a no-load state. 
<!-- Intialization of the dynamic stall states can be chosen to follow the suggested states by the individual model, or by using the steady state solution of the undeflected turbine. The structural states can be initialized at steady state of either no loading, gravity only, or with steady state aerodynamic loads. -->




## Future Coupling (Semi-tight Coupling)