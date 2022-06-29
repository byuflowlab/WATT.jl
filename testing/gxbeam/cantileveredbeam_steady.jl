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
twist = 0.0

### Create GXBeam Inputs to create function for DifferentialEquations
## Base inputs
nelem = 10 
rhub = 0
rtip = L

## DifferentialEquations inputs
tspan = (0.0, 8.0) #Todo. Why is it taking so much longer (in simulation time) to oscillate? -> The beam is really freaking heavy. With no gravity, there isn't really anything to help this thing accelerate. It's a beast. 

## Calculated inputs
r = range(0, L, length=nelem+2)
rvec = r[2:end-1] 
chordvec = ones(length(rvec)).*w
twistvec = zeros(length(rvec))
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
# fun = create_gxbeamfun(gxmodel, env, dsl; g=0.0, damping=false)
# diffvars = differentialvars(gxmodel) 


## Initialize
# x0 = initialize_gxbeam2(gxmodel, p, dsl)  

# ns = 18*gxmodel.ne + 12
# x0 = zeros(ns) #It appears that gxbeam's DAE starts from zeros. 
# dx0 = zeros(ns) 

# probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars)

# ## Solve
# sol_me = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01, reltol=1e-10) 









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
gxdaesystem, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions=prescribed_conditions, distributed_loads=distributed_loads)

pgxdae = (; prescribed_conditions=prescribed_conditions, distributed_loads=distributed_loads)

## construct a DAEProblem 
gxdae_prob = GXBeam.DAEProblem(gxdaesystem, assembly, tspan, pgxdae)


## solve DAEProblem
sol_gxbeam = DifferentialEquations.solve(gxdae_prob, DABDF2(), force_dtmin=true, dtmin=0.01) #Currently taking 18.78 seconds to run for a 8 second simulation... not bueno. 



history = [AssemblyState(gxdaesystem, assembly, sol_gxbeam[it]; prescribed_conditions) for it in eachindex(sol_gxbeam)]

nt = length(sol_gxbeam.t)
tback = 360
tidx = nt - tback

xd = [assembly.points[ipoint][1] + history[tidx].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflectiond = [history[tidx].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]


t = sol_gxbeam.t
def_tip = [history[it].points[end].u[2] for it = 1:nt]

tplt = plot(t, def_tip, xaxis="Time (s)", yaxis="Tip Deflection (m)")
display(tplt)







# ### Post Process my data
# history_me = [AssemblyState(system2, assembly, sol_me[it]; prescribed_conditions) for it in eachindex(sol_me)]

# nt_me = length(sol_me.t)
# tidx_me = nt_me - tback

# x_me = [assembly.points[ipoint][1] + history_me[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

# deflection_me = [history_me[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]






### Analytical Solution
Iz = w*(h^3)/12 #Second moment of area 
Iy = h*(w^3)/12

Izz, Iyy, Izy = rotate_smoa(Iz, Iy, 0.0, twist)

yfun(x) = load*x^2*(4*L*x -x^2 -6*L^2)/(24*E*Izz) #Shigley's
# Using Iz because the force in the -Y causes a negative moment about the Z axis. 

xvec = range(0, L, length=nelem+1)
yvec = yfun.(xvec)




### Plot
plt = plot(leg=:bottomleft, xaxis="Beam length (m)", yaxis="Beam Position (m)")
plot!(xvec, yvec, lab="Theory", markershape=:cross)
scatter!(x, deflection, lab="GXBeam")
# scatter!(xd, deflectiond, lab="Differential Equations") 
# scatter!(x_me, deflection_me, lab="Current", markershape=:x, markersize=8)
display(plt)