using Revise, DifferentialEquations, StaticArrays, FLOWMath, GXBeam, Plots, CurveFit, BenchmarkTools, LinearAlgebra, DelimitedFiles

#=
Test my function that creates a beam based on radial nodes (for the turbine description).
=#

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")



### Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
B = 3.0
pitch = 0.0
precone = 0.0 
yaw = 0.0*pi/180
tilt = 0.0 #5.0*pi/180
azimuth = 0.0
hubht = 90.0
shearexp = 0.0

vinf = 10.0
tsr = 7.55
rotorR = rtip*cos(precone)
omega = vinf*tsr/rotorR
pitch = 0.0
azimuth = 0.0*pi/180
rho = 1.225
mu=18.13e-6
a=343.0



thickvec = chordvec.*0.08 #Assume 8 percent thickness




### Create models
p, xp, xe = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

gxmodel = gxbeam(xp, xe)

assembly = create_gxbeam_assembly(gxmodel, p)

plt = plotassembly(assembly)
# display(plt)

nothing
