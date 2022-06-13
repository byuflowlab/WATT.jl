using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations, CCBlade, FLOWMath, CurveFit, NLsolve, DelimitedFiles, ForwardDiff

include("../src/extra.jl")
include("../src/blades.jl")
include("../src/environments.jl")
include("../src/bem.jl")
include("../src/riso.jl")
include("../src/bem-riso.jl")
include("../src/gxbeam.jl")
include("../src/static.jl")
include("../src/bem-gxbeam.jl")


function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end



#### Define variables. 
## Define simplified NREL 5MW Turbine constants and other info. 
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
thickvec = chordvec.*0.40495
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
tspan = (0.0, 1.0)
dt = 0.01

#### Look at an airfoil cross section for a reference thickness
# du40_a17 = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/Rotors/data/airfoils/DU40_A17.cor")


# airfoilplt = scatter(du40_a17[:,1], du40_a17[:,2], aspectratio=:equal)
# display(airfoilplt)

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
p_a = zeros(7*n)

for i = 0:n-1

    localpolar = hcat(airfoils[i+1].alpha, airfoils[i+1].cl, airfoils[i+1].cd)

    afs[i+1] = complexairfoil(localpolar)

    p_ccblade = [rvec[i+1], chordvec[i+1], twistvec[i+1], pitch, rhub, rtip, hubht]

    p_a[1+(7*i):7+(7*i)] = p_ccblade

end

blade = Blade(afs)

bemmodel = bem(;shearexp=shearexp)
dvb = differentialvars(bemmodel, n)

env = environment(rho, mu, a, vinf, omega, 0.0, 0.0)
gxmodel = gxbeam(n)
dvs = differentialvars(gxmodel)

_, p_s = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

p = vcat(p_a, p_s)

diffvars = vcat(dvb, dvs)

fun = create_bemgxbeamfun(rvec, chordvec, twistvec, rhub, rtip, bemmodel, blade, env)
fun2 = create_explicitbemgxfun(rhub, rtip, bemmodel, blade, env)


#### Initialize states 
println("Initialize states")
println("")
x0, dx0 = initialize_bemgx_states(bemmodel, gxmodel, env, blade, p)

x02 = x0[n+1:end]
dx02 = dx0[n+1:end]

fakeouts = zero(x0)
fo2 = zero(x02)

# fun(fakeouts, dx0, x0, p, 0.0) #Note: I expect this to be equal to zero... since the initialize_bemgx_states function uses the fixed point iteration static solve. So we should be at steady state... and shouldn't need any further convergence... unless my read states function is garbage (for GXBeam). But.... I don't think it is... since I used that function to validate my GXBeam wrapper. 

# @show fakeouts
# println("Call explicit function")
# println("")
# fun2(fo2, dx02, x02, p, 0.0)

# lb = fill(-Inf, length(x02))
# ub = fill(Inf, length(x02))
# solved = static_solve(fun2, x02, p, 0.0, lb, ub) #Solves, I don't know if it's a good solve. 

# @show fakeouts
# @show fo2


println("Create the DAE")
# probdae = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars=diffvars) #This runs... forever. And looking at phi, it looks like it is struggling to find a consistent set of states. -> After normalizing the GXBeam states like they should be (during initialization), this still didn't solve after 20+ minutes. 


probdae = DifferentialEquations.DAEProblem(fun2, dx02, x02, tspan, p, differential_vars=dvs) #phi is also getting shoved past 1... which means that the twist that I'm adding is... a lot? I'm not sure. Something is off with the twist. It must be. -> Well it is erroring with the tip correction... which is not handy. Is the solution just happen to through the tip correction off? What if I don't do a tip correction. -> I changed how the states were initialized (I normalized the GXBeam states like they should be), and this actually solved. 

savefile = "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/Rotors/testing/savingfile.csv"

prepareextra(tspan, savefile, n)


## Solve
println("Solve the DAE")
sol = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=dt, initializealg=NoInit())

savedmat = readextra(savefile) #-> Some of these inflow angles are ridiculous. Like 172 degrees. That's not realistic... right? That'd be psycho. So what the duece happened? 



## Step through the simulation
# integrator = init(probdae, DABDF2(); force_dtmin=true, dtmin=dt) #Well.... actually, we might be stalling out at this step. I'm not sure why though. The function goes pretty fast. -> Fixed the velocities and now this works.... but solve doesn't. It goes to a negative number. 

# integrator = init(probdae, DABDF2(); force_dtmin=true, dtmin=dt, initializealg=NoInit()) #This didn't work either.
# integrator = init(probdae, DABDF2(); force_dtmin=true, dtmin=dt, initializealg=ShampineCollocationInit()) #Changes both x0 and dx0 to find a consistent set of states. 

# step!(integrator) #This stalls out even (waited for a half hour). So it can't even take a single time step. 

# iterations_me = 1
# while integrator.t<tspan[2]
#     try
#         step!(integrator)
#     catch
#         println("Failed after $iter iterations")
#     end
#     global iterations_me += 1
# end

# lb_riso = fill(0.01, 4*n) #zeros(4*n)
# for i = 1:n
#     idx = 4*(i-1)
#     lb_riso[idx+1:idx+3] .= -Inf #I don't think states 1-3 of Riso matter if they go negative
# end
# lb_bem = zeros(n)
# lb_struc = fill(-Inf, 12+(18*n))

# lowbounds = vcat(lb_riso, lb_bem, lb_struc)
# upbounds = fill(Inf, length(x0))

# result = static_solve(fun, x0, p, 0.0, lowbounds, upbounds; iterations=10000) 
#Todo: It runs but doesn't converge. I think that since I'm forcing dx to be equal to zero, It can't solve because some of the states would have to be changing... although... moving doesn't mean accelerating. But moving does mean position is changing... so if there is any sort of change of position... then we'd have a problem. Although, it might be moving but not deflecting. 

############### Extract the data. 


### Create GXBeam structures to analyze with GXBeam. 
assembly = create_gxbeam_assembly(gxmodel, p_s)

prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) #Root section is fixed. 

fload = SVector(0.0, 0.0, 0.0)

distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> fload[2], fz= (s) -> fload[3]) for ielem in 1:n)


system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false)


history = [AssemblyState(system, assembly, sol[it]; prescribed_conditions) for it in eachindex(sol)]

tback = 0
nt_me = length(sol.t)
tidx_me = nt_me-tback

x_me = [assembly.points[ipoint][1] + history[tidx_me].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]

deflection_me = [history[tidx_me].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]


### Plot
# plt = plot(leg=:topleft, xaxis="Beam length (m)", yaxis="Beam Position (m)")
# scatter!(x_me, deflection_me, lab="Current", markershape=:x, markersize=8)
# display(plt)



# ################ Animate
# println("Beginning Animation. ")
# numframes = length(sol.t)
# tvec = range(tspan[1], tspan[2], length=numframes)
# solt = sol(tvec)
# historyt = [AssemblyState(system, assembly, solt[it]; prescribed_conditions) for it in eachindex(solt)]

# anim = @animate for i in 1:numframes
#     x_m = [assembly.points[ipoint][1] + historyt[i].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
#     deflection_m = [historyt[i].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]

#     scatter(x_m, deflection_m, leg=false, xaxis="Beam length (m)", yaxis="Beam Position (m)", markershape=:circle, markersize=8, xlims = (0, 65), ylims = (0, 30))
# end
# gif(anim, "BeamDeflection_3.0s.gif", fps=30)


#####Take a look at the forces. 


loadsT, loadsN = extractBEMloads(savefile, Array(sol)', bemmodel, blade, env, p)


# ###### Plot the forces
# x_elem = [assembly.elements[ielem].x[1] for ielem = 1:length(assembly.elements)]

# anim = @animate for i in 1:length(sol.t)
    
#     plt1 = plot(x_elem, loadsT[i,:], leg=:topright, xaxis="Beam length (m)", yaxis="Distributed Load (N)", lab="Tangential", ylims = (0,1500))
#     plt2 = plot(x_elem, loadsN[i,:], lab="Normal", ylims = (0,3e4))
#     plot(plt1, plt2, layout=(2,1))
# end
# gif(anim, "Beamloads_1.0s.gif", fps=15)



####### Compare with the fixed point solution. 
outs_fp, state_fp, system_fp, _, _, _, _, _ = fixedpoint(bemmodel, gxmodel, env, blade, p)


## Animate
x_elem = [assembly.elements[ielem].x[1] for ielem = 1:length(assembly.elements)]

anim = @animate for i in 1:length(sol.t)
    
    plt1 = plot(x_elem, loadsT[i,:], leg=:topright, xaxis="Beam length (m)", yaxis="Tangential Load (N)", lab="Dynamic", ylims = (0,1500))
    plot!(x_elem, outs_fp.Tp, lab="Static")

    plt2 = plot(x_elem, loadsN[i,:], lab="Dynamic", ylims = (0,3e4), yaxis="Normal Load (N)", legend=false)
    plot!(x_elem, outs_fp.Np, lab="Static")

    plot(plt1, plt2, layout=(2,1))
end
gif(anim, "dynamicbeamloads_wfixedpoint_1.0s.gif", fps=15)

### Todo: x_elem != rvec..... what's up with that? 

#### Visualize the deflections through Paraview

# root_chord = chordvec[1]
# tip_chord =  chordvec[end]
# airfoil = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/Rotors/data/airfoils/DU40_A17.cor")

# x = [assembly.points[ip][1] for ip in 1:length(assembly.points)]
# L = rtip-rhub

# sections = zeros(3, size(airfoil, 1), length(assembly.points))
# for ip = 1:length(assembly.points)
#     chord = root_chord * (1 - x[ip]/L) + tip_chord * x[ip]/L
#     sections[1, :, ip] .= 0
#     sections[2, :, ip] .= chord .* (airfoil[:,1] .- 0.5)
#     sections[3, :, ip] .= chord .* airfoil[:,2]
# end

# write_vtk("wind-turbine-blade-simulation", assembly, history, sol.t; sections = sections)








######### Test to see if a state from the explicit solve will work in the implicit solver. 
totalstates = hcat(savedmat[:,2:end], Array(sol)')
xs2 = vcat(savedmat[end,2:end], Array(sol)'[end,:])
dxs2 = mat_derivative(totalstates, savedmat[:,1])

# fakeouts = zero(x0)

# fun(fakeouts, dxs2[end,:], xs2, p, tspan[2]) 

# @show fakeouts  #It's fairly well converged. The aero residauls aren't great, but they aren't awful. 


###### Use initial conditions from explicit solve for an implicit solve. 
# tspan2 = (tspan[2], tspan[2]+1.0)

# println("Create the new DAE")
# probdae = DifferentialEquations.DAEProblem(fun, dxs2[end,:], xs2, tspan2, p, differential_vars=diffvars) 
 

# ## Solve
# println("Solve the new DAE")
# sol2 = DifferentialEquations.solve(probdae, DABDF2(), force_dtmin=true, dtmin=dt, initializealg=NoInit())








nothing