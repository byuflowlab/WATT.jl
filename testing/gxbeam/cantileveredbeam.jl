using Revise, DifferentialEquations, StaticArrays, FLOWMath, GXBeam, Plots, CurveFit, BenchmarkTools, LinearAlgebra, DelimitedFiles, Statistics

#=
Test that my coupling can handle a simple cantilevered beam with a constant distributed load. (That means no angular velocity, no gravity, and nothing is really changing with time (well the simulation will be through time, but my inputs will not change)).
=#

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../../src/gxbeam.jl")

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end

### Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
B = 3.0
pitch = 0.0
precone = 0.0 
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

thickvec = chordvec.*0.08 #Assume 8 percent thickness

tspan = (0.0, 5.0) #It's taking longer for half the amount of time. 





n = length(rvec)
twistvec = zeros(n) #See if getting rid of the twist will get rid of the twisting in the modell. -> It did, led me to find that my compliance matrix sucked. 
chordvec = ones(n).*mean(chordvec)
thickvec = ones(n).*mean(thickvec)

### Create models
n, p = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)
gxmodel = gxbeam(n)

# env = environment(rho, mu, a, vinf, omega, 0.0, 0.0)
env = environment(rho, mu, a, vinf, 0.0, 0.0, 0.0)

### Create distributed load
function dsl(t) 
    # f = zeros(3)
    f1 = 44.8*10 # approximately 10 pounds of thrust distributed
    return SVector(f1, 0.0, 0.0)
end

### Create gxbeam function. 
fun = create_gxbeamfun(gxmodel, env, dsl, g=0.0)
diffvars = differentialvars(gxmodel)
# nn = length(diffvars)

# x0 = initialize_gxbeam(gxmodel, p, dsl)
x0 = initialize_gxbeam2(gxmodel, p, dsl)

# x0 = ones(nn)  #It might not be solving because of the initial conditions. 
dx0 = zeros(length(x0))


probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars)

# sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=0.01)

# t = sol.t

# x = Array(sol)'

# endofpoints = 6*(n+1)

# uxtip = x[:,endofpoints-5]

# uxplt = plot(t, uxtip, xaxis="Time (s)", yaxis="Tip Deflection (m)", legend=:bottomright, lab="dt=0.01")
# display(uxplt)




### Solve the problem using GXBeam's differential equations solver
assembly = create_gxbeam_assembly(gxmodel, p)
elements = view(assembly.elements, :)

prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed

f = dsl(0.0)

# ## follower distributed loads
# f1 = f2 = @SVector zeros(3) 
# m1 = m2 = @SVector zeros(3) 
# m = @SVector zeros(3)

# distributed_loads = Dict() #Todo. Distributed loads are defined in the local reference frame, so if I want to apply forces in the overall reference frame, then I need to rotate back. 
# for ielem in 1:gxmodel.n
#     f1_f = f2_f = rotate_z(-twistvec[ielem])*elements[ielem].L*f/2
#     m1_f = m2_f = rotate_z(-twistvec[ielem])*elements[ielem].L*m/2
# 
#     f1_follower = f2_follower = SVector(f1_f[1], f1_f[2], f1_f[3])
#     m1_follower = m2_follower = SVector(m1_f[1], m1_f[2], m1_f[3])

#     distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
# end

distributed_loads = Dict(
    ielem => DistributedLoads(assembly, ielem; fx = (s) -> dsl(0)[1]) for ielem in
    1:n
)

# angular_velocity = SVector(omega, 0.0, 0.0)

#run initial condition analysis to get consistent set of initial conditions
system, converged = GXBeam.initial_condition_analysis(assembly, tspan[1]; prescribed_conditions, distributed_loads) #, angular_velocity


### construct a DAEProblem
prob = GXBeam.DAEProblem(system, assembly, tspan; prescribed_conditions, distributed_loads) #, angular_velocity

#### solve DAEProblem
# sol = DifferentialEquations.solve(prob) #IDAS Error
dtm = 0.0001

sol = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true, dtmin=dtm) #No change from 0.0001 to 0.00001

# sol = DifferentialEquations.solve(prob, alg_hints=[:stiff], force_dtmin=true) #Stalls out after 9.1e-5 seconds (IDAS error)

# sol = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true) #Still stalls out... I let it run 10+ minutes, and we didn't get very far. 

#This is taking a long time, I think because I used dead loads. Actually, it's probably from oscillations... which could come from the fact that I used a dead load and it's rotating.  -> Or from the fact that the beam is stiff and there isn't any damping. 

# sol2 = DifferentialEquations.solve(prob, DABDF2(), force_dtmin=true, dtmin=0.005) #Can artifically filter out oscillations by arbitrarilly increasing the time step. The larger the time step, the quicker it filters out. Converged to the same value. Did not converge to the same value that steady state solver did. 



t = sol.t
# t2 = sol2.t

x = Array(sol)'
# x2 = Array(sol2)'

endofpoints = 6*(n+1)

uxtip = x[:,endofpoints-5]
## uxtip2 = x2[:,endofpoints-5]

uxplt = plot(t, uxtip, xaxis="Time (s)", yaxis="Tip Deflection (m)", legend=:bottomright, lab="dt=$dtm", title="DiffEQ")
## plot!(t2, uxtip2, lab="dt=0.005")
display(uxplt)


############ Try solving the problem with just GXBeam.  

###### Steady analysis
system, converged = static_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false)

state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

xsteady = convert_assemblystate(state)

deltax = zeros(n+1)
for i = 2:n+1
    local idx = 6*(i-1)

    deltax[i] = xsteady[idx+1] #Some angular displacement in z (on points)
end

pointz = zeros(n+1) #Todo: I'm not sure that pointz should go from zero to the tip. Shouldn't it go from rhub to rtip. 
for i = 2:n+1
    local idx = 3*(i-1)

    pointz[i] = p[idx+3]
end

load = dsl(0.0)[1]
E = 6.83e10
Iyy = thickvec[1]*(chordvec[1]^3)/12 #Todo: isn't this the wrong I? Shouldn't it be Ixx
Ixx = chordvec[1]*(thickvec[1]^3)/12

L = rtip-rhub
yfun(x) = load*x^2*(4*L*x -x^2 -6*L^2)/(24*E*Iyy)
# yfun(x) = load*x^2*(4*L*x -x^2 -6*L^2)/(24*E*Ixx)

xvec = range(0, L, length=n)
yvec = yfun.(xvec)

steadyzplt = scatter(pointz, deltax, legend=:topleft, xaxis="Beam Length (m)", yaxis="X Displacement (m)", title="Steady GXBeam", lab="GXBeam")
plot!(xvec, abs.(yvec), lab="Theory")
display(steadyzplt)



#### Time simulation - GXBEAM SOLVED

t = range(tspan[1], tspan[2], length=2000)

system, history, converged = time_domain_analysis(assembly, t; prescribed_conditions=prescribed_conditions, distributed_loads = distributed_loads) 

idxvec = zeros(1)
for i=1:length(history)
    try
        local a = history[i]
    catch
        idxvec[1] = i-1
        break
    end
    idxvec[1] = i
end

idx = Int(idxvec[1])
newh = history[1:idx]
newt = t[1:idx] #The time domain analysis is diverging after 18 steps. I need to see how to get this not to diverge. 

# cyl = readdlm("../../data/airfoils/circle35.cor")
# af = readdlm("../../data/airfoils/naca2412.cor")

# sections = zeros(3, size(af,1), n+1)
# for ip = 1:n+1
#     if ip<=2
#         airfoil = cyl
#     else
#         airfoil = af
#     end
#     if ip==1
#         chord = chordvec[1]
#     else
#         chord = chordvec[ip-1]
#     end

#     sections[1, :, ip] .= chord .* (airfoil[:,1] .- 0.5)
#     sections[2, :, ip] .= chord .* airfoil[:,2]
#     sections[3, :, ip] .= 0.0
# end

### Twist the sections - It looks a little janky, but it's in there. 
# for i = 1:n+1
#     for j = 1:n+1
#         if j==1
#             twist = twistvec[j]
#         else
#             twist = twistvec[j-1]
#         end

#         sections[:,i,j] = rotate_z(-twist)*sections[:,i,j]
#     end
# end


# pathname = "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/bladeopt/coupling/coupling/mycoupling/vtk/"

# write_vtk(pathname*"shape_time_twist", assembly, newh, newt; sections=sections)

x = convert_history(newh)

endofpoints = 6*(n+1)

deltax_tip = x[:,endofpoints-5]

dtgx = t[3]-t[2]

tipdisplacementplt = plot(newt, deltax_tip, legend=true, xaxis="Time (s)", yaxis="Tip X Displacement (m)", title="GXBeam solved", lab="dt=$dtgx")
display(tipdisplacementplt)

nothing