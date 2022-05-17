using Plots, GXBeam, StaticArrays, LinearAlgebra, DifferentialEquations

include("../../src/blades.jl")
include("../../src/environments.jl")
include("../src/bem.jl")
include("../src/riso.jl")
include("../src/bem-riso.jl")
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