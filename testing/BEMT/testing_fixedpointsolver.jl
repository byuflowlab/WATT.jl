using Rotors, Plots, StaticArrays, LinearAlgebra, CCBlade, FLOWMath, CurveFit




#=
Adam Cardoza 10/07/22

Test the fixed point iteration solver I wrote to work with CCBlade. 
=#





function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end

path = dirname(@__FILE__)
cd(path)



#### Define variables. 
## Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [2.8667, 5.6, 8.333300000000001, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333, 62.9999]
chordvec = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.01, 2.764, 2.518, 2.313, 2.086, 1.419, 1.419]
twistvec =[ 0.23226841685540536, 0.23226841685540536, 0.23226841685540536, 0.23226841685540536, 0.20036379812894903, 0.17736035858766377, 0.15727161889720903, 0.136048415192958, 0.11421434625050891, 0.093567101199416, 0.07309438907352252, 0.0545415391248228, 0.0404741853537485, 0.02663372438543347, 0.015062191444711064, 0.006457718232379019, 0.0018500490071139892, 0.0018500490071139892]
B = 3.0
hubht = 90.0


## Airfoil Constants
A = [0.3, 0.7] #Aerodyn defaults
b = [0.13, 0.53]
Tp = 1.7
Tf = 3.0


## Turbine Control variables
pitch = 0.0
precone = 0.0 #2.5*pi/180 #Todo: !!!! I need to work in a way to include precone
yaw = 0.0*pi/180
tilt = 0.0 #5.0*pi/180
azimuth = 0.0

## Environmental variables
vinf = 10.0
tsr = 7.55
rotorR = rtip*cos(precone)
omega = vinf*tsr/rotorR
frequency = 1.0
amplitude = 0.0
rho = 1.225
mu = 18.13e-6
a = 343.0
shearexp = 0.0


ccpath = "../data/polars"

### Prep the ASD rotor and operating conditions 
aftypes = Array{AlphaAF}(undef, 8) 
aftypes[1] = AlphaAF(ccpath*"/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF(ccpath*"/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF(ccpath*"/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF(ccpath*"/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF(ccpath*"/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF(ccpath*"/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF(ccpath*"/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF(ccpath*"/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8, 8]


airfoils = aftypes[af_idx] 


############### Run CCBlade by itself  ###########
rotor = CCBlade.Rotor(rhub, rtip, B; precone=precone, turbine=true, tip=nothing)

sections = CCBlade.Section.(rvec, chordvec, twistvec, airfoils)

### Create Operating Point
operatingpoints = CCBlade.windturbine_op.(env.Vinf(0.0), env.RS(0.0), pitch, rvec, 0.0, 0.0, 0.0, 0.0, hubht, shearexp, env.rho)

outs = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

outsme = Rotors.fixedpointbem.(Ref(rotor), sections, operatingpoints, outs.phi) 

#=

=#





phiplt = plot(xaxis="Radius (m)", yaxis="Inflow Angle (rads)")
plot!(rvec, outs.phi, lab="CCBlade", markershape=:x)
plot!(rvec, outsme.phi, lab="Fixed-Point")
display(phiplt)


nothing