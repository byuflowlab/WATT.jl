# Rotors.jl
A coupling between CCBlade, dynamic stall models, and GXBeam. 


### Note for Developers
The general idea here is that we're going to couple together a couple of models, then DifferentialEquations is going to solve the states of the models. 

So for each model you're going to need a function to calculate the residuals such that R(dx, x, y, p, t, model), where dx is the state rates, x is the states, y is the time varying inputs, p is the constant inputs, t is the time, and model is... well, the model. Those residuals could be state rates or constraint equations. 

Then there can be different couplings. For each coupling you'll need: 
- A function that calculates the time varying inputs. Y(dx, x, p, t) This where the coupling between models happens. 
- A function that creates the DAE function to be solved by DifferentialEquations.jl. 
- A function that returns the differential variables. 
- A function to parse the solution output of DifferentialEquations.jl.
