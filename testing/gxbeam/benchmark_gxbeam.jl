using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations, FLOWMath, BenchmarkTools, Traceur #, Profile #, PProf

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")

function getfieldnames(x)
    return fieldnames(typeof(x))
end

#= Test whether or not my wrapper can create a beam and solve a simple spinning beam case, with twist.

Cantilevered beam with a constant distributed load, a constant cross sectional distribution, and a constant twist.... the only difference from the last test is that this includes an angular velocity. 
=#


## Base inputs
nelem = 10 
rhub = 0.0
rtip = 60.0

L = 60.0 #Length of beam
load = 100*4.48 # Newtons 
E = 6.83e10 #Young's Modulus
h = 0.25 # Thickness (meters)
w = 3.4 # Width (meters)
twist = pi/8

precone = 0.0
tsr = 7.55
vinf = 10.0
rotorR = rtip*cos(precone)
omega = vinf*tsr/rotorR


### Create GXBeam Inputs to create function for DifferentialEquations


## DifferentialEquations inputs
tspan = (0.0, 0.1)




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
env = environment(0.0, 0.0, 0.0, 0.0, omega, 0.0, 0.0) #I'm not using any of these inputs in this test. 

## Create distributed load
function dsl(t) 
    # f = rotate_x(pi/2-twist)'*[0.0, -load, 0.0]
    # f = [0.0, -load, 0.0]
    # return SVector(f[1], f[2], f[3])
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
if !@isdefined(sol)
    sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01) #Currently taking 4.6-4.8 seconds for 0.5 time domain simulation. 
end


xf = Array(sol)'[end,:]

function derivativeofstates(x, y, xi)  
    fit = Akima(x,y)
    return gradient(fit, xi)
end

dxf = [derivativeofstates(sol.t, Array(sol)'[:,i], sol.t[end]) for i in 1:length(xf)]

outs = zeros(length(dxf))

# @btime fun(outs, dxf, xf, p, 0.11)

xfp = SVector{length(xf)}(xf)
dxfp = SVector{length(xf)}(dxf)
pp = SVector{length(p)}(p) #That increased the number of allocations. 

@time fun(outs, dxfp, xfp, p, 0.11)  #0.140568 seconds (285.55 k allocations: 18.899 MiB, 98.18% compilation time) #Todo: Why is 98% compliation time? What is there to compile?.... Wait... Is this part of the global variable that I'm using? 
@time fun(outs, dxfp, xfp, p, 0.11)


# @profile (fun(outs, dxfp, xfp, p, 0.11))
# Profile.print() #Nothing useful is coming out. Why isn't it profiling? 

# @profile fun(outs, dxfp, xfp, p, 0.11)
# pprof()

# @pprof fun(outs, dxfp, xfp, p, 0.11) #Interesting, but didn't really show me anything spectacular. Something about size of line indicating memory usage. It seemed to be spending a lot of time in the function itself, and not in the places it calls. Which goes along with the idea that I've looked at the individual parts, and they are all fast. I wonder if the fact that this function is dependent on other functions makes it slower. Like.... If it was just one function, not a function within a function, would it be faster? 


@time sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01)

# @time sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01) #Runs faster on the second run when it doesn't have to compile. 


### Create GXBeam structures to analyze with GXBeam. 
# assembly = create_gxbeam_assembly(gxmodel, p)

# prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Root section is fixed. 

# fload = dsl(0)

# distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> fload[2], fz= (s) -> fload[3]) for ielem in 1:nelem)

# Omega = SVector(0.0, 0.0, omega)

# system2, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions, distributed_loads, angular_velocity = Omega)

# ## construct a DAEProblem
# prob = GXBeam.DAEProblem(system2, assembly, tspan; prescribed_conditions, distributed_loads, angular_velocity= Omega) 

## solve DAEProblem
# println("     GXBeam's DAE solved by DifferentialEquations:")
# @time sol_gxbeam = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true, dtmin=0.01) 



nothing
