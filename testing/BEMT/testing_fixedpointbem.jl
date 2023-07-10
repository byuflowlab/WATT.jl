using Plots, StaticArrays, LinearAlgebra, CCBlade, FLOWMath, CurveFit

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/fixedpointbem.jl")


#=
Adam Cardoza 8/4/22

Test the fixed point iteration BEM. 
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

# create airfoil array
airfoils = aftypes[af_idx]

n = length(rvec)
afs = Array{Airfoil}(undef, n)
p_a = zeros(7*n)

for i = 0:n-1
    localpolar = hcat(airfoils[i+1].alpha, airfoils[i+1].cl, airfoils[i+1].cd)

    afs[i+1] = complexairfoil(localpolar)
end

blade = Blade(rhub, rtip, rvec, afs)

env = environment(rho, mu, a, vinf, omega)


############### Considering the end of the blade, looking at input phi versus output phi and considering convergence.

# index = n
# phix = -0.3:0.01:0.3 #-pi/4:0.01:pi/4
# phiy = zero(phix)
# Fvec = zero(phix)

# state = nothing
# Vx = env.Vinf(0)
# Vy = env.RS(0)*rvec[index]

# for i = 1:length(phix)
#     phiy[i], Fvec[i] = update_inflowangle(phix[i], Vx, Vy, state, rvec[index], chordvec[index], twistvec[index], blade.airfoils[index], pitch, rhub, rtip, B, true, Prandtltiphub(), NoStall())
# end

   
# phiplt = plot(xaxis="Input ϕ (rads)", yaxis="Output ϕ (rads)",leg=:bottomright)
# plot!(phix, phix, lab="ϕ=ϕ")
# plot!(phix, phiy, lab="function")
# display(phiplt)

#=
BAHAHAHAHA, well... that didn't turn out so well. For the tip, we got that there are 7 unique solutions that the solver could potentially converge to... Right? Well... It can't converge to all of these. If it's coming from the left and it is below the line y=x, then it will move away from the solution. Similarly, if it is coming from the right and it is above the line x=y it should move away from the solution. Right? Because let's say we're at 1 (which the function value is above the line), then we'd get a higher phi value (3 to be exact), which would move us up to phi = 3. Handily phi = 3 sends us to phi = 0, which should converge to a solution... but it converges to the nearby solution, not to the correct solution. Plus, it can't converge to some of these that look like discontinuities. 

Oh look, you can see the discontinuities at pi/2. Duh... 

Yeah, zooming it to where the solution looks like it'll converge (near zero), I notice a couple of things. One, where the correct solution occurs there is a discontinuity from negative to positive. Two, the algorithm I've written will converge not to where the function crosses the line x=y, but to the point where the function crosses from positive to negative, right? 

No, that can't be right. It just has to be a point where it crosses x=y (but smoothly crosses). For instance at 0.2, on the right the solution will be forced left since it is below the line. On the left, the solution will be forced right because it is above the line. That's the convergence criteria. Above the line on the left of a crossing point, and below the line on the right. 

=#






phi, Fvec, converged = fixedpointbem(rvec, chordvec, twistvec, blade, pitch, rhub, rtip, B, env; turbine=true, tipcorrection=Prandtltiphub(), stallmodel=NoStall(), tolerance=1e-12, maxiterations=100, verbose=true)


#=
Currently it is just converging to a different phi then CCBlade.

-> There is a difference in the correction factor. Which is interesting, since it's a copy and paste. -> But that might be from the fact that it's a different phi.

I wonder if my tolerance criteria is bad. But.... this is what it is converging to. So if I used something like trying to converge the residual... I'm not confident that it would converge. 


    Ideas
    - I could compare induction velocities. 


I had copied and pasted the residual code from Github for CCBlade, then I realized that Dr. Ning had changed things and so some of the corrections were different than what my current CCBlade version was using. 
=#











############### Run CCBlade by itself  ###########
rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=true)

sections = CCBlade.Section.(rvec, chordvec, twistvec, airfoils)

### Create Operating Point
operatingpoints = CCBlade.windturbine_op.(env.Vinf(0.0), env.RS(0.0), pitch, rvec, 0.0, 0.0, 0.0, 0.0, hubht, shearexp, env.rho)

outs = CCBlade.solve.(Ref(rotor), sections, operatingpoints)






########### Compare a phi input to my function versus the CCBlade function. 


##### Correct phi = 0.03166366199021948
index = n

phi_in = 0.05

state = nothing
Vx = env.Vinf(0)
Vy = env.RS(0)*rvec[index]

rc, outc = CCBlade.residual(phi_in, rotor, sections[index], operatingpoints[index])


phim, Fm = update_inflowangle(phi_in, Vx, Vy, state, rvec[index], chordvec[index], twistvec[index], blade.airfoils[index], pitch, rhub, rtip, B, true, Prandtltiphub(), NoStall())

#### If phi = 0.05, index = n (so greater than the solution. )
# @show outc.phi, phim #Of course these are going to be different, one is the one we passed in and the other is the one we just calculated. 

# @show outc.F, Fm #The correction factor looks correct. 

# @show outc.alpha #These are the same. (I mean, they're negatives of each other, but CCBlade returns the negative of alpha that I'm showing (inside the function))

# @show outc.cn, outc.ct #Todo. Interesting these are different. 

#=
(cn, ct) = (0.1174041786922641, -0.000267587028855076) ## Me
(outc.cn, outc.ct) = (0.7558213391796139, 0.03251722528752242) ## CCBlade


By a lot. Mine is looking up the negative angle of attack and CCBlade is looking up the positive angle of attack (I checked)...

alpha = -0.04814995099288601
outc.alpha = 0.04814995099288601

julia> afs[end].cl(0.04814995099288601) #CCBlade's lookup (Is close to the Cn and ct. )
0.7565019432245693

julia> afs[end].cl(-0.04814995099288601) # My lookup (far from cn and ct)
0.11785963507499256

But I'm confused. CCBlade says it reports the negative alpha and cn. 

From CCBlade ...
    if rotor.turbine
        return R, Outputs(-Np, -Tp, -a, -ap, -u, -v, phi, -alpha, W, -cl, cd, -cn, -ct, F, G)
    else
        return R, Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)
    end

from my code: 

julia> rotor.turbine
true

There isn't any rotation, mach, or reynolds number corrections applied. 

But I did just find something important... that .... could have very well been screwing me up whenever I tried my own BEM residual.. .

From CCBlade...
    # airfoil cl/cd
    if rotor.turbine
        cl, cd = afeval(af, -alpha, Re, Mach)
        cl *= -1
    else
        cl, cd = afeval(af, alpha, Re, Mach)
    end

Which is just great.. So we need a negative alpha lookup and negative cl. Which would give me the correct alpha when do the airfoil lookup. 

Look at that. That matches CCBlade. Now to go see if the solutions will match. 

Wow... so much better. So so much better. Now I'm sure the phi solution space looks much different. But I should finish checking the intermediary values. 
Well... I don't calculate N, T, u, or v. So I guess a and ap? 

sigma and k are correct... 

=#

# @show outc.a, outc.ap

#= 
Empirical region
(a, ap) = (-0.9406090298691404, -0.5964024597948097)
(outc.a, outc.ap) = (0.6590249571151328, 0.39246748085918975)

Okay... so those bad boys are quite different as well. Maybe F is different? Nope. Those are the same. ... Plus I already checked those... So it is in the calculation of a and ap.

a and ap are returned negative in ccblade, so I expect the values to be negative in my code. But the values are still off. 

Ahh... Dr. Ning updated his code and I copied and pasted the updated code from GitHub, but my version was still from two years ago. But that still only gets a to match, ap is still off. 

Aha, I had the equation for ct wrong. Now let's check what happens. 
=#


# println(outs.F.-Fvec)


phiplt = plot(xaxis="Radius (m)", yaxis="Inflow Angle (rads)")
plot!(rvec, outs.phi, lab="CCBlade", markershape=:x)
plot!(rvec, phi, lab="Fixed-Point")
display(phiplt)

#= 
The first 3 stations are cylinders. 

=# 

nothing