using GXBeam, LinearAlgebra, StaticArrays
import ForwardDiff # for forward-mode automatic differentiation
import ReverseDiff # for reverse-mode automatic differentiation
using BenchmarkTools # for benchmarking function performance

L = 12 # inches
EI = 1e6
rhoA = 10

h = w = 1 # inches
# E = 30e6 # lb/in^4 Young's Modulus

# A = h*w
# Iyy = w*h^3/12
# Izz = w^3*h/12

# create points
nelem = 16
x = range(0, L, length=nelem+1)
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

# index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

# compliance matrix for each beam element
# compliance = fill(Diagonal([1/(E*A), 0, 0, 0, 1/(E*Iyy), 1/(E*Izz)]), nelem)

compliance = fill(Diagonal([1/EI, 1/EI, 1/EI, 1/EI, 1/EI, 1/EI]), nelem)
mass = fill(Diagonal([rhoA, rhoA, rhoA, 1e-6, 1e-6, 1e-6]), nelem)
dvec = [0.001*(@SVector ones(6)) for i in 1:nelem]

# create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop, compliance=compliance)

# cinternal = fill(Diagonal([1/EI, 1/EI, 1/EI, 1/EI, 1/EI, 1/EI]), nelem)
# cii = zeros(ForwardDiff.Dual{Float64}, 6, 6)
# cinternal = fill(cii, nelem)

# function update_compliance!(comp, x, k)
#     for i = 1:6
#         comp[k][i,i] = compliance[k][i,i]*x[k]
#     end
# end


# construct pfunc to overwrite prescribed conditions
pfunc = (p, t) -> begin

    # non-dimensional tip moment
    λ = p[1]

    # dimensionalized tip moment
    m = pi*E*Iyy/L
    M = λ*m

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # moment on right side
        nelem+1 => PrescribedConditions(Mz = M)
    )

    # compliancefactors = p[2:nelem+1] #Todo: We might not need to do a compliance factor, the problem was that it wasn't updating, not that it wasn't scaled properly. 

    # for i = 1:nelem
    #     update_compliance!(cinternal, compliancefactors, i)
    # end

    cinternal = [Diagonal([1/EI, 1/EI, 1/EI, 1/EI, 1/EI, 1/EI].*p[i+1]) for i in 1:nelem] #Todo: I need a more efficient way of allocating this. 

    @show typeof(cinternal)

    # for i = 1:nelem
    #     cinternal[i] = cinternal[i].*p[i+1]
    # end

    assembly_internal = Assembly(points, start, stop, compliance=cinternal, damping=dvec, mass=mass)

    # return named tuple with new arguments
    return (; assembly=assembly_internal, prescribed_conditions=prescribed_conditions)
end

# construct objective function
objfun = (p) -> begin

    # perform static analysis
    system, state, converged = static_analysis(assembly; pfunc, p)

    # return the desired outputs
    return [state.points[end].u[1], state.points[end].u[2]]
end

pvec = vcat(1.0, ones(nelem))

# compute sensitivities using ForwardDiff with λ = 1.0
ForwardDiff.jacobian(objfun, pvec)
#Todo: Compare against finite diff

#=
I'm not sure that I need to do the unsteady derivatives because they should behave the same as the steady ones... although... it's probably a good idea to try it. To make sure that I don't get a stack overflow error. 
=#