using Revise, DifferentialEquations, StaticArrays, FLOWMath, CCBlade, Plots, CurveFit, NLsolve

#=
Test that the ODE, DAE, time-marching ODE, and time-marching linearized Riso models are working. (All the same model, just implemented different ways). 

Adam Cardoza 8/2/22

TODO: I might want to delete this sucker. I think this is probably more of a DSM test. 
=#


include("../src/blades.jl")
include("../src/environments.jl")
include("../src/riso.jl")
include("../src/coupled.jl")
include("../src/solvers.jl")

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

### Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
B = 3.0
pitch = 0.0
precone = 2.5*pi/180 
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

frequency = 1.0
amplitude = 0.0

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
p_a = zeros(n)

for i = 0:n-1

    localpolar = hcat(airfoils[i+1].alpha, airfoils[i+1].cl, airfoils[i+1].cd)

    afs[i+1] = complexairfoil(localpolar)

    p_a[i+1] = chordvec[i+1]

end

blade = Blade(rhub, rtip, rvec/rtip, afs)

model = Riso()

env = environment(rho, mu, a, vinf, omega)

rfun = create_risofun(twistvec, blade, env, frequency, amplitude)
rdiffvars = differentialvars(model, n)

x0 = zeros(4*n)
dx0 = zeros(4*n)
tspan = (0.0, 20.0)

probdae = DifferentialEquations.DAEProblem(rfun, dx0, x0, tspan, p_a, differential_vars=rdiffvars)

sol = DifferentialEquations.solve(probdae)

t, Cl, Cd = parsesolution(model, blade, env, p_a, sol, twistvec, frequency, amplitude)

rcl = [blade.airfoils[1].cl(twistvec[1])]
tcl = [blade.airfoils[end].cl(twistvec[end])]


rode = create_risoODE(twistvec, blade, env, frequency, amplitude)
probode = DifferentialEquations.ODEProblem(rode, x0, tspan, p_a)
solode = DifferentialEquations.solve(probode)

tode, Clode, Cdode = parsesolution(model, blade, env, p_a, solode, twistvec, frequency, amplitude)




rode3 = createrisoode(blade)
p_aero3 = vcat(chordvec, twistvec, twistvec.*2, env.Vinf(0).*ones(n), zeros(n), pitch)
solver = BDF1() #It appears that BDF1 is the shortest solve, and only has a slight delay on coming to steady state. 

nt = length(tode)
x3 = zeros(nt, 4*n)
x3[1,:] = x0  

for i = 2:nt
    local dt = tode[i]-tode[i-1]
    x3[i,:] = solver(rode3, x3[i-1,:], p_aero3, tode[i-1], dt)
end

#=
It looks like it solves most of the states correctly, but some are blowing up. Like the third to last state looks like it exploded. .... So is it the method... or the solver? 

Well, after changing the solver to DiffEQ, all the entries in a simple print are a lot closer. Like, essentially the same... which tells me that it isn't the model, but is my solver... bummer. Because the DifferentialEquations solver is so much slower. 

-> So I have a couple of options. I could either use a better algorithm. Like, who knows how hard it is to write my own TSIT algorithm, or I could see if I could write a more efficient interface with DifferentialEquations. 

#TODO: Write a better interface with DifferentialEquations to leverage their superior algorithms
 
=#

cl3 = zeros(nt, n)
cd3 = zeros(nt, n)
for i = 1:nt
    for j = 1:n
        y = [env.Vinf(0), 0.0, 0.0, 0.0, twistvec[j], 0.0]
        idx = 4*(j-1)
        cl3[i,j], cd3[i,j] = riso_coefs(x3[i,idx+1:idx+4], y, chordvec[j], blade.airfoils[j])
    end
end








pathname = "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/bladeopt/coupling/coupling/mycoupling/figures/riso/"

clplt = plot(xaxis="time (s)", yaxis="Coefficient of Lift", legend=:bottomright, title="Root Analysis") 
plot!(t, Cl[:,1], lab="DAE solution")
plot!(tode, Clode[:,1], lab="ODE solution")
plot!(tode, cl3[:,1], lab="Time-Marching")
hline!(rcl, lab="Steady value")
display(clplt)
# savefig(pathname*"rootcl_steady_33022.png")

cltipplt = plot(xaxis="time (s)", yaxis="Coefficient of Lift", legend=:bottomright, title="Tip Analysis") 
plot!(t, Cl[:,end], lab="DAE solution")
plot!(tode, Clode[:,end], lab="ODE solution")
plot!(tode, cl3[:,end], lab="Time-Marching")
hline!(tcl, lab="Steady value")
display(cltipplt)
# savefig(pathname*"tipcl_steady_33022.png")





# pathname = "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/bladeopt/coupling/coupling/mycoupling/figures/riso/"
# anim = @animate for i = 1:length(t)
#     ti = t[i]
#     plot(rvec, Cl[i,:], xaxis="Radial Location", yaxis="Coefficient of Lift", lab="t=$ti", ylims=(0.0, 0.9))
# end
# gif(anim, pathname*"anim_fps15.gif", fps = 5)

nothing