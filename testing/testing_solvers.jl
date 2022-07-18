using DifferentialEquations, Plots

include("../src/solvers.jl")

function testode(y, p, t)
    dy = zeros(2)
    dy[1] = y[1]
    dy[2] = y[1] - y[2]
    return dy
end

function analytical(x)
    y1 = exp(x)
    y2 = exp(x)/2 + 3*exp(-x)/2
    return [y1, y2]
end

y0 = [1, 2]

solver = RK4()


x = 0:0.001:1 # I don't think it has anything to do with the time step. 
p = nothing
nx = length(x)

yrk4 = zeros(nx, 2)
yana = zeros(nx, 2)

yrk4[1,:] = y0
yana[1,:] = y0

for i=2:nx
    dx = x[i]-x[i-1]
    yrk4[i,:] = solver(testode, yrk4[i-1,:], p, x[i-1], dx)
    yana[i,:] = analytical(x[i])
end


prob = ODEProblem(testode, y0, (x[1], x[end]))
sol = solve(prob)
yDE = Array(sol)'

plt = plot(x, yrk4[:,1], lab="y1 - RK4", leg=:topleft)
plot!(x, yrk4[:,2], lab="y2 - RK4")
plot!(x, yana[:,1], lab="y1 - ana")
plot!(x, yana[:,2], lab="y2 - ana")
plot!(sol.t, yDE[:,1], lab="y1 - DE")
plot!(sol.t, yDE[:,2], lab="y2 - DE")
display(plt)


nothing