#=
Testing to see if I can get derivatives through GXBeam. 

Adam Cardoza 
=#
println("Precomiling...")
using GXBeam, FiniteDiff, ForwardDiff, LinearAlgebra, Plots, StaticArrays, Statistics, BenchmarkTools


#### First run a simple simulation (static)
println("Generating assembly...")

L = 10
EI = 1e6

rhoA = 10

load = 1000


# create points
nelem = 8
x = range(0, L, length=nelem+1)
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

# index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

# compliance matrix for each beam element
compliance = fill(Diagonal([1/EI, 1/EI, 1/EI, 1/EI, 1/EI, 1/EI]), nelem)
mass = fill(Diagonal([rhoA, rhoA, rhoA, 1e-6, 1e-6, 1e-6]), nelem)



# create assembly of interconnected nonlinear beams
dvec = [0.001*(@SVector ones(6)) for i in 1:nelem]
assembly = Assembly(points, start, stop, compliance=compliance, damping=dvec, mass=mass)

# pre-initialize system storage
system = StaticSystem(assembly)



# create dictionary of prescribed conditions
prescribed_conditions = Dict(1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    nelem+1 => PrescribedConditions(Fz = load))

println("Static solution...")
# perform a static analysis
system, state, converged = static_analysis(assembly;
    prescribed_conditions = prescribed_conditions,
    linear = false)


def_z = [state.points[i].u[3] for i in eachindex(points)]

r = [points[i][1] for i in eachindex(points)]

defplt = plot(r, def_z, leg=false, xaxis="Radius (m)", yaxis="Deflection (m)")
# display(defplt)


###### Simple time domain simulation
tvec = 0:0.001:.25

println("Time domain solution...")
system, history, converged = time_domain_analysis(assembly, tvec;
    prescribed_conditions = prescribed_conditions,
    structural_damping = true, steady_state=false)

tipdef_z = [history[i].points[end].u[3] for i in eachindex(tvec)]

tipdefplt = plot(tvec, tipdef_z, leg=false, xaxis="Time (s)", yaxis="Tip Deflection (m)")
display(tipdefplt)

# cii = zeros(ForwardDiff.Dual{ForwardDiff.Tag{typeof(objective), Float64}, Float64, 8}, 6, 6)
# cinternal = fill(cii, nelem)
# cinternal = fill(Diagonal([1/EI, 1/EI, 1/EI, 1/EI, 1/EI, 1/EI]), nelem)
tvec_internal = 0:0.001:.25

# function update_compliance!(comp, x, k)
#     for i = 1:6
#         comp[k][i,i] = compliance[k][i,i]*x[k]
#     end
# end

# function objective(x)

#     # @show x[1]

#     # @show cinternal[1]
#     # @show cinternal[end]
#     for i = 1:nelem
#         # @show i
#         # for j = 1:6
#         #     cinternal[i][j,j] = compliance[i][j,j]*x[i]
#         # end

#         # cinternal[i] = compliance[i].*x[i] #Question: Why does this work and not what's above? -> Hopefully there aren't other inplace modifications that aren't broken. 
#         update_compliance!(cinternal, x, i)
#     end
#     # @show cinternal[1]

#     # ci = compliance.*x

#     # create assembly of interconnected nonlinear beams
#     assembly_internal = Assembly(points, start, stop, compliance=cinternal, damping=dvec, mass=mass)

#     _, history_internal, _ = time_domain_analysis(assembly_internal, tvec_internal;
#     prescribed_conditions = prescribed_conditions,
#     structural_damping = true, steady_state=false)

#     # @show history_internal[end].points[end].V

#     return history_internal[end].points[end].V[3]
#     # return history_internal[end].points[end].u[3]
# end

# construct pfunc to overwrite prescribed conditions
pfunc = (p, t) -> begin


    # compliancefactors = p[2:nelem+1] #Todo: We might not need to do a compliance factor, the problem was that it wasn't updating, not that it wasn't scaled properly. 

    cinternal = [Diagonal([1/EI, 1/EI, 1/EI, 1/EI, 1/EI, 1/EI].*p[i]) for i in 1:nelem] #Todo: I need a more efficient way of allocating this. 

    assembly_internal = Assembly(points, start, stop, compliance=cinternal, damping=dvec, mass=mass)

    # return named tuple with new arguments
    return (; assembly=assembly_internal)
end

# construct objective function
objfun = (p) -> begin

    _, history_internal, _ = time_domain_analysis(assembly, tvec_internal;
    prescribed_conditions = prescribed_conditions,
    structural_damping = true, steady_state=false, pfunc, p)

    return [history_internal[end].points[end].V[3]]
end

# x0 = [compliance[i][1,1] for i in 1:nelem]
x0 = ones(nelem)

xfake = [0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1]
# @show objective(xfake)

println("Calculating gradients...") 
println("   Finite diff")
dVdx_fd = FiniteDiff.finite_difference_jacobian(objfun, x0)
@btime FiniteDiff.finite_difference_jacobian(objfun, x0)
println("   Forward diff")
dVdx_fad = ForwardDiff.jacobian(objfun, x0)
@btime ForwardDiff.jacobian(objfun, x0)

@show dVdx_fd
@show dVdx_fad

#=
7/28/23 For whatever reason, the compliances of the first 15 elements seem to have no effect on the final position of the final element. That seems... wrong. 

8/16/23 I can't seem to directly deduce what's going. Tyler sugggested that maybe the influence of the preceeding elements was being dwarfed by the size of the step. So I tried normalizing. Perhaps the way I normalized didn't get the outcome I wanted... but the fact that the valeus are zero point zero, it makes me think that something else is going on. 

8/17/23 The reason why the derivative wasn't getting calculated correctly was because the internal compliance matrix wasn't updating for whatever reason. Now the derivatives make sense. 
=#

nothing