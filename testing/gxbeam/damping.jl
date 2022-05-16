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
twist = 0.0

## DifferentialEquations inputs
tspan = (0.0, 0.5)
# tspan = (0.0, 5.0)

dt = 0.01
# dt = 0.1

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


outs = zeros(length(x0))

# fun(outs, dx0, x0, p, 0.0)


## Solve
sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=dt) 
### 

states = Array(sol)'

t1plt = plot(sol.t, states[:,8])
display(t1plt)



### Create GXBeam structures to analyze with GXBeam. 
assembly = create_gxbeam_assembly(gxmodel, p)

prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Root section is fixed. 

distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> -load) for ielem in 1:nelem)

system, converged = static_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false)

state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)


## Post process GXBeam solution
x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]






### Post Process my data
## run initial condition analysis to get consistent set of initial conditions
system2, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions, distributed_loads)

history_me = [AssemblyState(system2, assembly, sol[it]; prescribed_conditions) for it in eachindex(sol)] #Todo: I'm not sure that this is doing what I want it to. The solution for the tip is quite different from the states. (The solution across time.) I wonder if system2 screws with what I want. 

nt_me = length(sol.t)
tidx_me = nt_me

x_me = [assembly.points[ipoint][1] + history_me[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_me = [history_me[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]

ytip = [history_me[it].points[end].u[2] for it = 1:length(sol.t)]

tplt = plot(sol.t, ytip, leg=false, xaxis="Time (s)", yaxis="Tip Deflection (m)")
display(tplt)







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
scatter!(x_me, deflection_me, lab="Current", markershape=:x, markersize=8)
# display(plt)
nothing 