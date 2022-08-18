# For Developers


## Code Outline

Here I'll describe the how the code is put together for easier debugging and additions to the code.



Here I'll describe how the models are coupled together. 


## Current Coupling (Loosely Coupled)
As is done in OpenFAST, this package provides a loose coupling between the aerodynamic models, and the structural model. Environmental flow conditions and static airfoil data is used to calculate the inflow angle of the rotor (by solving the residual formed by the blade element and momentum equations as done by CCBlade.). That inflow angle is used to calculate the angle of attack of the airfoil. 

Instead of looking up the lift and drag to calculate the loadings, that angle of attack is passed to the dynamic stall model. The states of the dynamic stall model are updated (based on the angle of attack and the environmental conditions). Based on the dynamic states, the lift and drag are used to calculate the loading on the rotor. Those loadings are fed to GXBeam to simulate the structural response for a given time step. This process is repeated throughout time. 

Intialization of the dynamic stall states can be chosen to follow the suggested states by the individual model, or by using the steady state solution of the undeflected turbine. The structural states can be initialized at steady state of either no loading, gravity only, or with steady state aerodynamic loads. 

The aerodynamic states are the states of the given dynamic stall model, which are given in dynamicstallmodels.jl. The structural states are given by GXBeam.jl. 


## Future Coupling (Semi-tight Coupling)