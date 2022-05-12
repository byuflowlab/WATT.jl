using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")

function getfieldnames(x)
    return fieldnames(typeof(x))
end

#= Test whether or not my wrapper can create a beam and solve a simple beam case. 

Cantilevered beam with a constant distributed load. 
=#

L = 60 #Length of beam
load = 100*4.48 # Newtons 
E = 6.83e10 #Young's Modulus
h = 0.25 # Thickness (meters)
w = 3.4 # Width (meters)
Iz = w*(h^3)/12 #Second moment of area
Iy = h*(w^3)/12

### Create GXBeam Inputs to create function for DifferentialEquations
## Base inputs
nelem = 10 
rhub = 0
rtip = L

## DifferentialEquations inputs
tspan = (0.0, 10.0)

## Calculated inputs
r = range(0, L, length=nelem+2)
rvec = r[2:end-1] 
chordvec = ones(length(rvec)).*w
twistvec = zeros(length(rvec))
thickvec = ones(length(rvec)).*h

## Create parameters
n, p = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

## Create models
gxmodel = gxbeam(n)
env = environment(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) #I'm not using any of these inputs in this test. 

## Create distributed load
function dsl(t) 
    return SVector(0.0, -load, 0.0)
end

## Create gxbeam function. 
fun = create_gxbeamfun(gxmodel, env, dsl, g=0.0)
diffvars = differentialvars(gxmodel)


## Initialize
x0 = initialize_gxbeam2(gxmodel, p, dsl)
dx0 = zeros(length(x0))

probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars)


## Solve
sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01) 
### #Instability detected... 








### Create GXBeam structures to analyze with GXBeam. 
assembly = create_gxbeam_assembly(gxmodel, p)

prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Root section is fixed. 

distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> -load) for ielem in 1:nelem)

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
sol_gxbeam = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true, dtmin=0.01) #It bottoms out and hits dtmin... 

# x_gx = Array(sol_gxbeam)' #Well... ain't that interesting. His DAE formulation doesn't use as many states. -> length(system.x) = 132, length(system2.x) = 192, 6*(nelem + 1) + 18*(nelem) = 246
#= 
From the element.jl file, I see that there are a different number of state variables for a static element and a dynamic element. Looks like static has 12 state variables, and dynamic has 18. ... So why does that second system have only 192.... It should be as many as my equations. 

Aight, so I asked Taylor, GXBeam drops points states if there are only two beams attached to a point and there are no forces or constraints on the point. So in this instance, all of the points except for the first one won't have state variables. Since I'm constraining u and theta for the first point, there are six states that are required. For the 10 elements, I have 18 states per element, so 180... which means.... there are another 6 states somewhere that I don't know about. 

=#

history = [AssemblyState(system2, assembly, sol_gxbeam[it]; prescribed_conditions) for it in eachindex(sol_gxbeam)]

nt = length(sol_gxbeam.t)
tidx = nt-556

xd = [assembly.points[ipoint][1] + history[tidx].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflectiond = [history[tidx].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]


t = sol_gxbeam.t
def_tip = [history[it].points[end].u[2] for it = 1:nt]

tplt = plot(t, def_tip, xaxis="Time (s)", yaxis="Tip Deflection (m)")
# display(tplt)







### Post Process my data
history_me = [AssemblyState(system2, assembly, sol[it]; prescribed_conditions) for it in eachindex(sol)]

nt_me = length(sol.t)
tidx_me = nt_me-556

x_me = [assembly.points[ipoint][1] + history[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_me = [history[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]









### Analytical Solution
Izz, Iy, Izy = rotate_smoa(Iz, Iy, 0.0, pi/2 - twist)

yfun(x) = load*x^2*(4*L*x -x^2 -6*L^2)/(24*E*Izz) #Shigley's
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