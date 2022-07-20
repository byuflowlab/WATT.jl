using DifferentialEquations, FLOWMath, CCBlade, GXBeam, LinearAlgebra, Plots, StaticArrays, CurveFit

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")
include("../src/riso.jl")
include("../src/bem-riso.jl")
include("../src/gxbeam.jl")
include("../src/coupled.jl")
include("../src/solvers.jl")
include("../src/extra.jl")
include("../src/static.jl")




#### Define variables. 
## Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
# rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
# chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
# twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
# thickvec = chordvec.*0.08
# B = 3.0
# hubht = 90.0

# rvec = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
# chordvec = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
# twistvec = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

thickvec = chordvec.*0.08
B = 1.0
hubht = 90.0


## Airfoil Constants
A = [0.3, 0.7] #Aerodyn defaults
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
tspan = (0.0, 20.0)
dt = 0.05

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
# af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

# create airfoil array
airfoils = aftypes[af_idx]

n = length(rvec)
afs = Array{Airfoil}(undef, n)
p_a = zeros(7*n)

for i = 0:n-1

    localpolar = hcat(airfoils[i+1].alpha, airfoils[i+1].cl, airfoils[i+1].cd)

    afs[i+1] = complexairfoil(localpolar)

    p_ccblade = [rvec[i+1], chordvec[i+1], twistvec[i+1], pitch, rhub, rtip, hubht]

    p_a[1+(7*i):7+(7*i)] = p_ccblade

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


p_s, points, elements = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)
gxmodel = gxbeam(points, elements)

p = vcat(p_a, p_s)

## Create system indices  
start = 1:gxmodel.ne
stop = 2:gxmodel.np

## Create assembly
assembly = create_gxbeam_assembly(gxmodel, p_s, start, stop) 

outs, state, system, assembly, prescribed_conditions, converged, iters, resids, distributed_load = fixedpoint(bemmodel, gxmodel, env_steady, blade, p; maxiterations = 100, verbose = false, tolerance= 1e-12, g=0.0)


t0 = tspan[1]
cchistory, gxhistory, tvec = simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, assembly, tspan, dt; verbose=false, coupling=TwoWay(), b=5000)





# outs, state, system, assembly, prescribed_conditions, converged, iters, resids, distributed_load = fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations = 100, verbose = true, tolerance= 1e-12) #Todo The static solution is mucho different from the time domain simulation... even though I thought I was passing in the same assembly... but obviously not. That is .... troublesome. It makes sense that the fixedpoint solution is near the initial time domain solution. They are both steady state solutions. -> Aight... See below (not the one directly below, but the one below that). 

#Todo. Another thing that is interesting is that this FP simulation is very different from the one in the other file.... -> Aha, the fixed point solution uses the environment struct at t=0... which the velocity is super low... which makes sense. Okay.. 

# Todo. Now the time domain solution is pretty far off from the steady solution. Well... right from the outset, the solution is off. -> Alright. So I was using a follower load in the fixed point iteration, which changes the solution... right. Then the difference at that point is the fact that I'm plotting the solution after several iterations, and the scale I'm looking at shows the difference between just a one-way coupling and a fixed point iteration. 
# -> Oh dang, also... I was showing the end element's deflection versus the end point's deflection... so yeah.. that's the real reason why it's off (I had forgot that a single iteration of the fixed point iteration gets it right to the fourth decimal point (at least for this case).).

nt = length(tvec)

tipdef = [gxhistory[i].points[end].u[2] for i in 1:nt]

tipplt = plot(tvec, tipdef, leg=:topright, xaxis="time (s)", yaxis="deflection (m)", lab="dynamic")
hline!([state.points[end].u[2]], lab="static")
# hline!([0.1339370376174049]) #On top of the static solution. (First iteration from static solution.)
display(tipplt) 

#Alright, the tip is not displacing any further because the beam is starting from steady state. Which is what we were looking for. 

# x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

# deflection = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]

# bodydefplt = plot(x, deflection, xaxis="beam length (m)", yaxis="Deflection (m)", lab="simulate")
# display(bodydefplt)

# anim = @animate for i in 1:nt
#     state = gxhistory[i]
#     x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

#     deflection = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]
#     plot(x, deflection, leg=false, xaxis="Beam length (m)", yaxis="Deflection (m)", ylims= (-0.05, 0.65), xlims=(0, rtip))
# end
# gif(anim, "bodydeflections_threeway_071822.gif", fps = 100)

# anim = @animate for i in 1:nt
#     ccout = cchistory[i]
#     state = gxhistory[i]
#     x = [assembly.elements[ipoint].x[1] + state.elements[ipoint].u[1] for ipoint = 1:length(assembly.elements)]
    
#     plot(x, ccout.Np, leg=:topleft, xaxis="Beam length (m)", yaxis="Distributed Loading (N/m)", ylims= (-500, 8000), xlims=(0, rtip), lab="Normal")
#     plot!(x, ccout.Tp, lab="Tangent")
# end
# gif(anim, "bodyforces_threeway_071822.gif", fps = 100)



# nothing