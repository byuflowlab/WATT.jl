using DifferentialEquations, FLOWMath, CCBlade, GXBeam, LinearAlgebra, Plots, StaticArrays, CurveFit, NLsolve

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")
include("../src/riso.jl")
include("../src/bem-riso.jl")
include("../src/gxbeam.jl")
include("../src/coupled.jl")
include("../src/aeroonly.jl")
include("../src/solvers.jl")
include("../src/extra.jl")
include("../src/static.jl")




#### Define variables. 
## Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0


rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

thickvec = chordvec.*0.08
B = 3
hubht = 90.0


## Airfoil Constants
# A = [0.3, 0.7] #Aerodyn defaults
A = [0.29, 0.33]
b = [0.13, 0.53]
Tp = 1.7
Tf = 3.0


## Turbine Control variables
pitch = 0.0
precone = 0.0 #2.5*pi/180 #TODO: !!!! I need to work in a way to include precone
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

#### Define solution
tspan = (0.0, 4.0)
dt = 0.001

### Prep the ASD rotor and operating conditions 
aftypes = Array{AlphaAF}(undef, 8)
aftypes[1] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = [3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]


# create airfoil array
airfoils = aftypes[af_idx]

n = length(rvec)
afs = Array{Airfoil}(undef, n)

for i = 0:n-1

    localpolar = hcat(airfoils[i+1].alpha, airfoils[i+1].cl, airfoils[i+1].cd)

    afs[i+1] = complexairfoil(localpolar; A=A)


end

rR = rvec./rtip
blade = Blade(rhub, rtip, rR, afs)

dsmodel = Riso()
bemmodel = createbem(;shearexp=shearexp)

function rampingU(rho, mu, a, vinf, Omega, ramptime)
    function U(t)
        if t<ramptime
            return t*vinf/ramptime + 0.1
        else
            return vinf
        end
    end

    function dU(t)
        if t<ramptime
            return vinf/ramptime
        else
            return 0.0
        end
    end

    ufun(t) = SVector(U(t), 0.0, 0.0)
    omegafun(t) = SVector(0.0, 0.0, 0.0)
    udotfun(t) = SVector(dU(t), 0.0, 0.0)
    omegadotfun(t) = SVector(0.0, 0.0, 0.0)
    Vinf(t) = U(t)
    RS(t) = Omega
    Vinfdot(t) = dU(t)
    RSdot(t) = 0.0
    return Environment(rho, mu, a, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end

env_steady = environment(rho, mu, a, vinf, omega)
env = rampingU(rho, mu, a, vinf, omega, 0.5)





t0 = tspan[1]
cchistory, dshistory, tvec = simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env_steady, tspan, dt, B; solver=DiffEQ(), verbose=false, coupling=ThreeWay())

#=
8/1/22 Adam Cardoza
Current it's converging to a solution, but the inboard nodes are super negative and the outboard nodes aren't nearly large enough in magnitude. 

Additionally, the initial solve is doing some interesting things... with the exact same problem. 

-> Well... If I run the steady environment instead of the ramping environment, then the loading at least starts normal (Then it appears to blow up at some point). 

-> Well.... I was doing Vy super wonky. So I fixed that. But it still appears to blow up.

8/2/22 
- Well bummer. I was checking to see if the Riso model is fine, and it appears to solve just fine. ARgh. Which is a pain. Well.... Maybe I revisit the tightly coupled aero-only solution and see if I can't get that going for Teagan. 

=#




nt = length(tvec)



dsstates = dshistory[1]

dsstate = dsstates[end]

alfa = collect(0:.1:14)
alpha = alfa.*(pi/180)
nalf = length(alfa)


static = afeval.(Ref(aftypes[8]), alpha, Ref(0.0), Ref(0.0))
dynamic = afeval.(Ref(dsstate), alpha, Ref(0.0), Ref(0.0))

clstatic = [static[i][1] for i in 1:nalf]
cldyn = [dynamic[i][1] for i in 1:nalf]

# dsclplt = plot(xaxis="Angle of Attack (deg)", yaxis="Coefficient of Lift")
# plot!(alfa, clstatic, lab="Static")
# plot!(alfa, cldyn, lab="Dynamic")
# display(dsclplt)

#Note: The lift doesn't change with alpha when A1 + A2 = 1. So I started using the coefficients from Hansen's paper. Now they change with alpha, somewhat, but it appears to be a linear change... which I guess makes sense. Judging by the equation for alpha_e. --> I really need to settle how I'm coupling the BEM and dynamic stall model. I'll need to sit down with Dr. Ning. 


# dsplt = plotdshistory(dshistory, tvec, 14)
# display(dsplt)

# for i = 1:n #Todo: I'm afraid I don't understand enough on the dynamic stall model why this is failing. 
#     dsplt = plotdshistory(dshistory, tvec, i; titletext="Station $i")
# display(dsplt)
# end



# tidx = 2
# ti = tvec[tidx]
# clplt = plot(xaxis="Blade radius (m)", yaxis="Distributed Load (N/m)", title="t = $ti")
# plot!(rvec, cchistory[tidx].Np, lab="Normal")
# plot!(rvec, cchistory[tidx].Tp, lab="Tangent")
# display(clplt)

#=
- Checking time index 2, it appears to blow up instantly. Which means either the time step is too large... or there is something broken. -> Well, changing the time step didn't do jack squat. Which means that there is something wrong between the first and second solve. Which makes sense, because that's when I first use the dynamic stall states to call CCBlade, and I don't believe that I've tested getting the lift and drag. 

=#



nothing