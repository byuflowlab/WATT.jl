using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations, FLOWMath, NLsolve

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")
include("../../src/static.jl")
include("../../src/extra.jl")

function getfieldnames(x)
    return fieldnames(typeof(x))
end

#= Test whether or not my wrapper can create a beam and solve a simple spinning beam case, with twist.

Cantilevered beam with a constant distributed load, a constant cross sectional distribution, and a constant twist.... the only difference from the last test is that this includes an angular velocity. 
=#


## Base inputs
nelem = 10 
rhub = 0
rtip = 60

L = 60 #Length of beam
load = 100*4.48 # Newtons 
E = 6.83e10 #Young's Modulus
h = 0.25 # Thickness (meters)
w = 3.4 # Width (meters)
twist = 0.0

precone = 0.0
tsr = 7.55
vinf = 10.0
rotorR = rtip*cos(precone)
omega = vinf*tsr/rotorR


### Create GXBeam Inputs to create function for DifferentialEquations


## DifferentialEquations inputs
tspan = (0.0, 6.0)




## Calculated inputs
r = range(0, L, length=nelem+2)
rvec = r[2:end-1] 
chordvec = ones(length(rvec)).*w
twistvec = ones(length(rvec)).*twist
thickvec = ones(length(rvec)).*h

## Create parameters
p, xp, xe = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

## Create models
gxmodel = gxbeam(xp, xe)
env = environment(0.0, 0.0, 0.0, 0.0, omega, 0.0, 0.0) #I'm not using any of these inputs in this test. 

## Create distributed load
function dsl(t) 
    # f = rotate_x(pi/2-twist)'*[0.0, -load, 0.0]
    f = [0.0, -load, 0.0]
    return SVector(f[1], f[2], f[3])
end











### Create GXBeam structures to analyze with GXBeam. 
assembly = create_gxbeam_assembly(gxmodel, p)

prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Root section is fixed. 

fload = dsl(0)

distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> fload[2], fz= (s) -> fload[3]) for ielem in 1:nelem)

# distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> -load) for ielem in 1:nelem) #The previous line acheives the same output. 
Omega = SVector(0.0, 0.0, omega)

g = 0.0
grav(t) = SVector(-g*cos(pi/2 - env.Omega(t)), -g*sin(pi/2 - env.Omega(t)), 0.0)

system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega, gravity=grav)

state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)


## Post process GXBeam solution
x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]


force_scaling = GXBeam.default_force_scaling(assembly)
xs = convert_assemblystate(state, assembly) #./force_scaling








### Test GXBeam DAE formulation
## run initial condition analysis to get consistent set of initial conditions
system2, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions, distributed_loads, angular_velocity = Omega, gravity=grav)

## construct a DAEProblem
prob = GXBeam.DAEProblem(system2, assembly, tspan; prescribed_conditions, distributed_loads, angular_velocity= Omega, gravity=grav) 

## solve DAEProblem
sol_gxbeam = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true, dtmin=0.01) 



## Initialize
# x0 = initialize_gxbeam2(gxmodel, p, dsl)
x0 = sol_gxbeam[1]
dx0 = zeros(length(x0))



###### Test if GXBeam's dynamic residual returns zero with dx = 0
gxfun = prob.f

fakegxouts = zero(x0)

gxfun(fakegxouts, dx0, xs, prob.p, 0.0) #GXBeam's dynamic residual also doesn't return steady values. 

# @show fakegxouts


#### Try a static solve of GXBeam's dynamic residual. 
nn = length(x0)
lb = fill(-Inf, nn)
ub = fill(-Inf, nn)
solved = static_solve(gxfun, xs, prob.p, 0.0, lb, ub)
xgxs = solved.zero

stategx = AssemblyState(system, assembly, xgxs; prescribed_conditions)

x_static = [assembly.points[ipoint][1] + stategx.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_static = [stategx.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]






### Post Process my data
## Create gxbeam function. 
fun = create_gxbeamfun(gxmodel, env, dsl; g=0.0, b=0.0)
diffvars = differentialvars(gxmodel)


probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars)


## Solve
sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01) #Currently taking 4.6-4.8 seconds for 0.5 time domain simulation. 

history_me = [AssemblyState(system2, assembly, sol[it]; prescribed_conditions) for it in eachindex(sol)]

tback = 325
nt_me = length(sol.t)
tidx_me = nt_me-tback

x_me = [assembly.points[ipoint][1] + history_me[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_me = [history_me[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]






####### Test if the steady state produces a zero residual
fakeouts = zero(x0)

fun(fakeouts, dx0, xgxs, p, sol.t[tidx_me]) #using the scaled states is much better. There are still some discrepancies, but much much better. I wonder if I use the statically solved states. -> Would you look at that. It's almost magical. Okay, so I need to add some scaling to GXBeam states. I'll add that in the initialize functions and functions that read states. 

# @show fakeouts


######## Test if any of my states produce a zero residual
tvec = 0:0.01:tspan[2]
xme = Array(sol(tvec))'
dxme = derivative_me(sol, tvec)

fakeouts2 = zero(x0)

fun(fakeouts2, dxme[end, :], xme[end,:], p, tvec[end])

# @show fakeouts2 #This checks out. The residual is zero. Which makes sense. The solver can solve my system. 


### Analytical Solution
boxx = [w/2, w/2, -w/2, -w/2]
boxy = [-h/2, h/2, h/2, -h/2]

points = zeros(4,3)
points[:,1] = boxx
points[:,2] = boxy

points[1,:] = rotate_z(pi/2 - twistvec[1])*points[1,:]
points[2,:] = rotate_z(pi/2 - twistvec[2])*points[2,:]
points[3,:] = rotate_z(pi/2 - twistvec[3])*points[3,:]
points[4,:] = rotate_z(pi/2 - twistvec[4])*points[4,:]

# Aplt = plot(points[:,1], points[:,2], aspectratio=:equal)
# display(Aplt)

Iz = w*(h^3)/12 #Second moment of area
Iy = h*(w^3)/12

Iz1, Iy1, Izy = secondmomentofarea(points[:,1], points[:,2])

Iz2, Iy2, Izy2 = rotate_smoa(Iz, Iy, 0.0, pi/2 - twist)

EIz = assembly.elements[1].compliance[6,6]
EIy = assembly.elements[1].compliance[5,5]

Iz3 = 1/(EIz*E) #Todo: Correct values, but swapped with each other. 
Iy3 = 1/(EIy*E)

# println("Outside functions (script)")
# @show Iz, Iy
# @show Iz1, Iy1
# @show Iz2, Iy2
# @show Iz3, Iy3

load2 = load + 2.69e3*g*w*h

yfun(x) = load*x^2*(4*L*x -x^2 -6*L^2)/(24*E*Iy1) #Shigley's
# Using Iz because the force in the -Y causes a negative moment about the Z axis. 

xvec = range(0, L, length=nelem+1)
yvec = yfun.(xvec)






# ### Plot
plt = plot(leg=:bottomleft, xaxis="Beam length (m)", yaxis="Beam Position (m)")
# plot!(xvec, yvec, lab="Theory", markershape=:cross)
scatter!(x, deflection, lab="Steady")
scatter!(x_me, deflection_me, lab="Current", markershape=:x, markersize=8)
scatter!(x_static, deflection_static, lab="Static")
display(plt)

nothing