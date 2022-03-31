using Revise, DifferentialEquations, StaticArrays, FLOWMath, CCBlade, Plots, CurveFit

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")


### Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
B = 3.0
pitch = 0.0
precone = 2.5*pi/180 #Todo: !!!! I need to work in a way to include precone
yaw = 0.0*pi/180
tilt = 0.0 #5.0*pi/180
azimuth = 0.0
hubht = 90.0
shearexp = 0.0

tspan = (0.0, 10.0)
vinf = 10.0
tsr = 7.55
rotorR = rtip*cos(precone)
omega = vinf*tsr/rotorR
pitch = 0.0
azimuth = 0.0*pi/180
rho = 1.225
mu=18.13e-6
a=343.0

A = [0.3, 0.7] #Aerodyn defaults
b = [0.13, 0.53]
Tp = 1.7
Tf = 3.0

center = 0.0
frequency = 0.0
omegadot = 0.0

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
    afs[i+1] = simpleairfoil(localpolar)

    p_ccblade = [rvec[i+1], chordvec[i+1], twistvec[i+1], pitch, rhub, rtip, hubht]

    p_a[1+(7*i):7+(7*i)] = p_ccblade
end

blade = Blade(afs)

model = bem(;shearexp=shearexp)

env = environment(rho, mu, a, vinf, omega, 0.0, 0.0)

bfun = create_bemfun(model, blade, env)
bdiffvars = differentialvars(model, n)

x0 = twistvec
dx0 = zeros(n)
tspan = (0.0, 2.0)

probdae = DifferentialEquations.DAEProblem(bfun, dx0, x0, tspan, p_a, differential_vars=bdiffvars)

sol = DifferentialEquations.solve(probdae)

# residsplt = plot(sol) #States are steady, as expected. 
# display(residsplt)

phi, N, T, Thrust, Torque = parsesolution(model, blade, env, p_a, sol)

#### Compare to CCBlade
B = 1
rotor = Rotor(rhub, rtip, B, precone=precone, turbine=true)

sections = Section.(rvec, chordvec, twistvec, airfoils)

op = windturbine_op.(vinf, omega, pitch, rvec, precone, yaw, tilt, 0.0, hubht, shearexp, rho)

out = CCBlade.solve.(Ref(rotor), sections, op)


pathname = "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/bladeopt/coupling/coupling/mycoupling/figures/bem/"


Nplt = plot(xaxis="Blade radius", yaxis="Normal Force (N)", legend=:bottomright)
plot!(rvec, N, lab="DAE solve")
plot!(rvec, out.Np, lab="CCBlade")
display(Nplt) #They match
# savefig(pathname*"normalforce_33022.png")

Tplt = plot(xaxis="Blade radius", yaxis="Tangential Force (N)", legend=:bottomright)
plot!(rvec, T, lab="DAE solve")
plot!(rvec, out.Tp, lab="CCBlade")
display(Tplt) #They match
# savefig(pathname*"tangentforce_33022.png")

nothing
