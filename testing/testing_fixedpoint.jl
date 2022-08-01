using Plots, GXBeam, StaticArrays, LinearAlgebra, CCBlade, FLOWMath, CurveFit, OpenFASTsr

include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")
# include("../src/riso.jl")
# include("../src/bem-riso.jl")
include("../src/gxbeam.jl")
# include("../src/aerostructural.jl")
include("../src/coupled.jl")
include("../src/static.jl")

#=
Adam Cardoza 8/1/22

Test the fixed point iteration solution (Steady state coupling between CCBlade and GXbeam that iterates back and forth until the solution converges to a given level.). Use the stiffness and mass matrices from OpenFAST's NREL 5MW wind turbine. 

Dr. Ning isn't sure that this will return the same solution as the steady state of a coupling between CCBlade, a dynamic stall model, and GXBeam. 

Note that this solution might be off from the OpenFAST solution... if the polars aren't the same. 

Note: You need to run the OpenFAST simulation in the openfast directory before running this file (you need the out files). In order to run the OpenFAST simulation, in the terminal within the directory run ``openfast NREL5MW_input.fst``
=#


of = OpenFASTsr


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
rvec = [1.5, 2.8667, 5.6, 8.333300000000001, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333, 62.9999]
chordvec = [3.542, 3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.01, 2.764, 2.518, 2.313, 2.086, 1.419, 1.419]
twistvec =[0.23226841685540536, 0.23226841685540536, 0.23226841685540536, 0.23226841685540536, 0.23226841685540536, 0.20036379812894903, 0.17736035858766377, 0.15727161889720903, 0.136048415192958, 0.11421434625050891, 0.093567101199416, 0.07309438907352252, 0.0545415391248228, 0.0404741853537485, 0.02663372438543347, 0.015062191444711064, 0.006457718232379019, 0.0018500490071139892, 0.0018500490071139892]
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
af_idx = [1, 1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8, 8]

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




ofpath = "./OpenFAST_NREL5MW" 

# adfile = of.read_adfile("NREL5MW_ADfile.dat", ofpath)
# adblade = of.read_adblade("NREL5MW_adblade.dat", ofpath)
edfile = of.read_edfile("NREL5MW_edfile.dat", ofpath)
# bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
bdblade = of.read_bdblade("NREL5MW_bdblade.dat", ofpath)

assembly = of.make_assembly(edfile, bdblade)
ne = length(assembly.elements)

outs, state, system, converged, iters, resids = fixedpoint(rvec, chordvec, twistvec, pitch, rhub, rtip, assembly, env, blade, B; maxiterations=100, tolerance=1e-3, verbose=true)






######## Parse the solution ##############
## Obtain the deflections
def_x = [state.elements[ielem].u[1] for ielem = 1:ne]
def_y = [state.elements[ielem].u[2] for ielem = 1:ne]
def_z = [state.elements[ielem].u[3] for ielem = 1:ne]

x = [assembly.elements[ielem].x[1] + state.elements[ielem].u[1] for ielem = 1:ne]


x_assembly = [assembly.elements[ielem].x[1] for ielem = 1:ne]








############### Run CCBlade by itself  ###########
rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=true)

sections = CCBlade.Section.(rvec, chordvec, twistvec, airfoils)

### Create Operating Point
operatingpoints = CCBlade.windturbine_op.(env.Vinf(0.0), env.RS(0.0), pitch, rvec, 0.0, 0.0, 0.0, 0.0, hubht, shearexp, env.rho)

outs_ccblade = CCBlade.solve.(Ref(rotor), sections, operatingpoints)



########## Read in the OpenFAST simulation data ##########



############# Visualize the loads and deflections. ############
defplt = plot(xaxis="Radial Location (m)", yaxis="Element Deflection (m)", legend=:topleft)
plot!(x, def_x, lab="X deflection")
plot!(x, def_y, lab="Y deflection")
plot!(x, def_z, lab="Z deflection")
display(defplt)

loadsplt = plot(xaxis = "Radial Location (m)", yaxis = "Element Load (N)", legend=:topleft)
plot!(rvec, outs.Np, lab="Flapwise")
plot!(rvec, outs.Tp, lab="Lead-Lag")
plot!(rvec, outs_ccblade.Np, lab="Flapwise - CCBlade") 
plot!(rvec, outs_ccblade.Tp, lab="Lead-Lag - CCBlade")
display(loadsplt)

#=
It looks like most of the convergence occurs in the GXBeam solution. The max difference between the initial CCBlade solution and the fixed point iteration is like .16 N... which is tiny. 
=#

nothing