using Plots, GXBeam, StaticArrays, LinearAlgebra, CCBlade, FLOWMath, CurveFit, BenchmarkTools

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")
# include("../src/riso.jl")
# include("../src/bem-riso.jl")
include("../src/gxbeam.jl")
# include("../src/aerostructural.jl")
include("../src/static.jl")



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
B = 1.0
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


env = environment(rho, mu, a, vinf, omega, 0.0, 0.0)
gxmodel = gxbeam(n)


_, p_s = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

p = vcat(p_a, p_s)




outs, state, system, assembly, prescribed_conditions, converged, iters, resids = fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations = 4, verbose = true, tolerance= 1e-12)

# @btime fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations = 4, verbose = false, tolerance= 1e-12) ### 12.281 ms (21900 allocations: 9.92 MiB)
# state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

## Obtain the deflections
def_x = [state.elements[ielem].u[1] for ielem = 1:n]
def_y = [state.elements[ielem].u[2] for ielem = 1:n]
def_z = [state.elements[ielem].u[3] for ielem = 1:n]

x = [assembly.elements[ielem].x[1] + state.elements[ielem].u[1] for ielem = 1:n]





### Run CCBlade by itself
rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bemmodel.turbine)

sections = CCBlade.Section.(rvec, chordvec, twistvec, airfoils)

### Create Operating Point
operatingpoints = CCBlade.windturbine_op.(env.U(0.0), env.Omega(0.0), pitch, rvec, 0.0, 0.0, 0.0, 0.0, hubht, bemmodel.shearexp, env.rho)

outs_ccblade = CCBlade.solve.(Ref(rotor), sections, operatingpoints)








defplt = plot(xaxis="Radial Location (m)", yaxis="Element Deflection (m)", legend=:topleft)
plot!(x, def_x, lab="X deflection")
plot!(x, def_y, lab="Y deflection")
plot!(x, def_z, lab="Z deflection")
display(defplt)

loadsplt = plot(xaxis = "Radial Location (m)", yaxis = "Element Load (N)", legend=:topleft)
plot!(x, outs.Np, lab="Flapwise")
plot!(x, outs.Tp, lab="Lead-Lag")
display(loadsplt)

nothing