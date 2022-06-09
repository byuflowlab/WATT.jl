using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations, CCBlade, FLOWMath, CurveFit, NLsolve

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")
include("../src/riso.jl")
include("../src/bem-riso.jl")
include("../src/gxbeam.jl")
include("../src/static.jl")
include("../src/bem-gxbeam.jl")


function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end



#### Define variables. 
## Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
thickvec = chordvec.*0.08
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

#### Define solution
tspan = (0.0, 0.5)
dt = 0.01

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
p_a = zeros(7*n)

for i = 0:n-1

    localpolar = hcat(airfoils[i+1].alpha, airfoils[i+1].cl, airfoils[i+1].cd)

    afs[i+1] = complexairfoil(localpolar)

    p_ccblade = [rvec[i+1], chordvec[i+1], twistvec[i+1], pitch, rhub, rtip, hubht]

    p_a[1+(7*i):7+(7*i)] = p_ccblade

end

blade = Blade(afs)

bemmodel = bem(;shearexp=shearexp)
dvb = differentialvars(bemmodel, n)

env = environment(rho, mu, a, vinf, omega, 0.0, 0.0)
gxmodel = gxbeam(n)
dvs = differentialvars(gxmodel)

_, p_s = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

p = vcat(p_a, p_s)

diffvars = vcat(dvb, dvs)

fun = create_bemgxbeamfun(rvec, chordvec, twistvec, rhub, rtip, bemmodel, blade, env)
fun2 = create_explicitbemgxfun(rhub, rtip, bemmodel, blade, env)


#### Initialize states 
println("Initialize states")
println("")
x0, dx0 = initialize_bemgx_states(bemmodel, gxmodel, env, blade, p)

x02 = x0[n+1:end]
dx02 = dx0[n+1:end]

fakeouts = zero(x0)
fo2 = zero(x02)

# fun(fakeouts, dx0, x0, p, 0.0) #Note: I expect this to be equal to zero... since the initialize_bemgx_states function uses the fixed point iteration static solve. So we should be at steady state... and shouldn't need any further convergence... unless my read states function is garbage (for GXBeam). But.... I don't think it is... since I used that function to validate my GXBeam wrapper. 

# println("Call explicit function")
# println("")
fun2(fo2, dx02, x02, p, 0.0)

# @show fakeouts
@show fo2


println("Create the DAE")
# probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars) #This runs... forever. And looking at phi, it looks like it is struggling to find a consistent set of states.

probdae = DifferentialEquations.DAEProblem(fun2, dx02, x02, tspan, p, differential_vars=dvs) #phi is also getting shoved past 1... which means that the twist that I'm adding is... a lot? I'm not sure. Something is off with the twist. It must be. -> Well it is erroring with the tip correction... which is not handy. Is the solution just happen to through the tip correction off? What if I don't do a tip correction. 


## Solve
println("Solve the DAE")
# sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=dt, initializealg=NoInit())





## Step through the simulation
# integrator = init(probdae, DABDF2(); force_dtmin=true, dtmin=dt) #Well.... actually, we might be stalling out at this step. I'm not sure why though. The function goes pretty fast. -> Fixed the velocities and now this works.... but solve doesn't. It goes to a negative number. 

# integrator = init(probdae, DABDF2(); force_dtmin=true, dtmin=dt, initializealg=NoInit()) #This didn't work either.
# integrator = init(probdae, DABDF2(); force_dtmin=true, dtmin=dt, initializealg=ShampineCollocationInit()) #Changes both x0 and dx0 to find a consistent set of states. 

# step!(integrator) #This stalls out even (waited for a half hour). So it can't even take a single time step. 

# iterations_me = 1
# while integrator.t<tspan[2]
#     try
#         step!(integrator)
#     catch
#         println("Failed after $iter iterations")
#     end
#     global iterations_me += 1
# end

# lb_riso = fill(0.01, 4*n) #zeros(4*n)
# for i = 1:n
#     idx = 4*(i-1)
#     lb_riso[idx+1:idx+3] .= -Inf #I don't think states 1-3 of Riso matter if they go negative
# end
# lb_bem = zeros(n)
# lb_struc = fill(-Inf, 12+(18*n))

# lowbounds = vcat(lb_riso, lb_bem, lb_struc)
# upbounds = fill(Inf, length(x0))

# result = static_solve(fun, x0, p, 0.0, lowbounds, upbounds; iterations=10000) 
#Todo: It runs but doesn't converge. I think that since I'm forcing dx to be equal to zero, It can't solve because some of the states would have to be changing... although... moving doesn't mean accelerating. But moving does mean position is changing... so if there is any sort of change of position... then we'd have a problem. Although, it might be moving but not deflecting. 

nothing