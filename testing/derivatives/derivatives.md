# Derivatives

The goal is to pass derivatives through Rotors via AD. Here is the game plan.

- Set up the file to run a single step at a time
- Set up the wind turbine to have 20 stations
- Decide an objective 
  - Maybe root bending moment... or damage equivalent moment? 
- Determine design variables 
    - Set twist as my design variable
- Take derivatives of a single step by finite difference
- Take derivatives of the same step by AD. 
- Compare

I took the analyze_simpleturbine file and cleaned it up. Then I converted it to run an objective function. I decreased the discretization both spatially and temporally for the beam. Note that the temporal discretization can't get much larger that dt=0.05. The aerostructural coupling for OpenFAST has a cow.... I don't know where Rotor's starts to break down. (Especially since I don't have much in the ways of checking if the coupling is doing well). 

For my objective I decided the average bending moment across time. It isn't perfect, but hey it was something quick. I'm sure that whatever I set my objective as has a large effect on what the derivative value is. It might make more sense to convert to damage equivalent moment (DEM), just because that seems like a more applicable time varying quantity. I did set my DVs as twist.... in degrees apparently... whatever. I don't think it'd make that big of a difference for now. 

I have the finite difference gradient. I don't know how accurate it is. I remember Dr. Ning's book mentioning a mathematical way of checking the derivative... but I need to look into that. 
-> The test I was thinking of was called the dot-product test. And I think that is only really useful for checking adjoint formulations, not the derivatives itself. So I think I'm kind of stuck not using that. But the book does suggest using complex-step and AD, if possible. 