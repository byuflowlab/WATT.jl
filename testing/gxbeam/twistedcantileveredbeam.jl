using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")

function getfieldnames(x)
    return fieldnames(typeof(x))
end

#= Test whether or not my wrapper can create a beam and solve a simple beam case, with twist 

Cantilevered beam with a constant distributed load, a constant cross sectional distribution, and a constant twist.

Note that this should be just another case of the simple beam, but now we have a new Izz and Iyy. 
=#

L = 60 #Length of beam
load = 100*4.48 # Newtons 
E = 6.83e10 #Young's Modulus
h = 0.25 # Thickness (meters)
w = 3.4 # Width (meters)
twist = pi/8

Iz = w*(h^3)/12 #Second moment of area
Iy = h*(w^3)/12

### Create GXBeam Inputs to create function for DifferentialEquations
## Base inputs
nelem = 10 
rhub = 0
rtip = L

## DifferentialEquations inputs
tspan = (0.0, 0.5) #Note: I feel like this is much slower. Like I know it is slower, and I don't have a solid number to give. But without twist it could solve the steady state, both time simulations in just a matter of moments, for 10 seconds. For this one, it didn't solve 10 seconds in quite a few minutes. Like.... Yeah... It doesn't seem to be oscillating a lot. 

#Also... the deflection isn't as much as theory... but it does oscillate about the steady state GXBeam solution. I wonder why it is taking so much longer. #Todo: I wonder if I personally rotate the compliance matrix, then compare... if I get what theory predicts.... another way to ask the question, is when you rotate the beam in GXBeam, how does it change the properties? 

## Calculated inputs
r = range(0, L, length=nelem+2)
rvec = r[2:end-1] 
chordvec = ones(length(rvec)).*w
twistvec = ones(length(rvec)).*twist
thickvec = ones(length(rvec)).*h

## Create parameters
n, p = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

## Create models
gxmodel = gxbeam(n)
env = environment(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) #I'm not using any of these inputs in this test. 

## Create distributed load
function dsl(t) 
    # f = rotate_x(pi/2 - twist)'*[0.0, -load, 0.0]
    f = [0.0, -load, 0.0]
    return SVector(f[1], f[2], f[3])
end

## Create gxbeam function. 
fun = create_gxbeamfun(gxmodel, env, dsl, g=0.0, damping=false)
diffvars = differentialvars(gxmodel)


## Initialize
x0 = initialize_gxbeam2(gxmodel, p, dsl)
dx0 = zeros(length(x0))

probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars)


## Solve
sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01) #Currently taking 4.6-4.8 seconds for 0.5 time domain simulation. 









### Create GXBeam structures to analyze with GXBeam. 
assembly = create_gxbeam_assembly(gxmodel, p)

prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Root section is fixed. 

fload = dsl(0)

distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> fload[2], fz= (s) -> fload[3]) for ielem in 1:nelem)

# distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> -load) for ielem in 1:nelem) #The previous line acheives the same output. 

system, converged = static_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false)

state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)


## Post process GXBeam solution
x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]






### Test GXBeam DAE formulation
## run initial condition analysis to get consistent set of initial conditions
system2, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions, distributed_loads)

## construct a DAEProblem
prob = GXBeam.DAEProblem(system2, assembly, tspan; prescribed_conditions, distributed_loads) 

## solve DAEProblem
sol_gxbeam = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true, dtmin=0.01) #Currently taking 0.08 seconds for a 0.5 second simulation.... so I have some speeding up to do on my implementation. -> Although... I think it might just be slower. Because we're putting information for the assembly in P... so it might inherently be slower. I need to talk to Dr. Ning about this. 



history = [AssemblyState(system2, assembly, sol_gxbeam[it]; prescribed_conditions) for it in eachindex(sol_gxbeam)]

nt = length(sol_gxbeam.t)
tidx = nt-15

xd = [assembly.points[ipoint][1] + history[tidx].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflectiond = [history[tidx].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]


t = sol_gxbeam.t
def_tip = [history[it].points[end].u[2] for it = 1:nt]

tplt = plot(t, def_tip, xaxis="Time (s)", yaxis="Tip Deflection (m)")
display(tplt)







### Post Process my data
history_me = [AssemblyState(system2, assembly, sol[it]; prescribed_conditions) for it in eachindex(sol)]

nt_me = length(sol.t)
tidx_me = nt_me-15

x_me = [assembly.points[ipoint][1] + history_me[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_me = [history_me[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]









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

yfun(x) = load*x^2*(4*L*x -x^2 -6*L^2)/(24*E*Iz1) #Shigley's
# Using Iz because the force in the -Y causes a negative moment about the Z axis. 

xvec = range(0, L, length=nelem+1)
yvec = yfun.(xvec)






### Plot
plt = plot(leg=:bottomleft, xaxis="Beam length (m)", yaxis="Beam Position (m)")
plot!(xvec, yvec, lab="Theory", markershape=:cross)
scatter!(x, deflection, lab="GXBeam")
scatter!(xd, deflectiond, lab="Differential Equations") # The Differential Equations solution passes through the steady state solution. Who knows if it converges to the steady state solution... cause... I don't want to look through to see if that works. At least, not right now. 
scatter!(x_me, deflection_me, lab="Current", markershape=:x, markersize=8)
display(plt)