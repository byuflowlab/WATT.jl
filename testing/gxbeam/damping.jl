using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")

function getfieldnames(x)
    return fieldnames(typeof(x))
end

#= Test whether or not my wrapper can create a beam and solve a simple beam case, but let's include some damping!!! 

Cantilevered beam with a constant distributed load. 
=#

L = 60 #Length of beam
load = 100*4.48 # Newtons 
E = 6.83e10 #Young's Modulus
h = 0.25 # Thickness (meters)
w = 3.4 # Width (meters)

### Create GXBeam Inputs to create function for DifferentialEquations
## Base inputs
nelem = 10 
rhub = 0
rtip = L
twist = 0.0

## DifferentialEquations inputs
tspan = (0.0, 30.0)
# tspan = (0.0, 5.0)

dt = 0.01
# dt = 0.1

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
env = environment(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) #I'm not using any of these inputs in this test. 

## Create distributed load
function dsl(t) 
    return SVector(0.0, -load, 0.0)
end

## Create gxbeam function. 
fun = create_gxbeamfun(gxmodel, env, dsl, g=0.0, b=5000.0) #Note: For b=5000.0, that takes an awfully long time to go to steady state. We're talking 25 ish seconds. I don't have enough experience to say whether or not that would be at steady state. And I feel like if the value is super huge, then it could cause some problems. 
diffvars = differentialvars(gxmodel)








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
## Initialize
x0 = initialize_gxbeam2(gxmodel, p, dsl)
# x0 = sol_gxbeam[1]

dx0 = zeros(length(x0))

probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars)


outs = zeros(length(x0))



## Solve
sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=dt) 



## run initial condition analysis to get consistent set of initial conditions
system2, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions, distributed_loads)

history_me = [AssemblyState(system2, assembly, sol[it]; prescribed_conditions) for it in eachindex(sol)] 

nt_me = length(sol.t)
tidx_me = nt_me

x_me = [assembly.points[ipoint][1] + history_me[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_me = [history_me[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]

ytip = [history_me[it].points[end].u[2] for it = 1:length(sol.t)]

tplt = plot(sol.t, ytip, leg=:topright, xaxis="Time (s)", yaxis="Tip Deflection (m)", lab="Current")
hline!([deflection[end]], lab="Static GXBeam")
display(tplt)







### Analytical Solution
Iz = w*(h^3)/12 #Second moment of area
Iy = h*(w^3)/12
Izz, Iy, Izy = rotate_smoa(Iz, Iy, 0.0, twist)

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