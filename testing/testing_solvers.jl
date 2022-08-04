using DifferentialEquations, Plots, NLsolve, Sundials

include("../src/solvers.jl")

#=
Adam Cardoza 8/3/22

Test the ODE and DAE solvers. 

=#

########### Test the ODE Solvers ###########

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
bdfsolve = BDF1()


x = 0:0.001:1 # I don't think it has anything to do with the time step. 
p = nothing
nx = length(x)

yrk4 = zeros(nx, 2)
ybdf1 = zeros(nx, 2)
yana = zeros(nx, 2)

yrk4[1,:] = y0
ybdf1[1,:] = y0
yana[1,:] = y0

for i=2:nx
    dx = x[i]-x[i-1]
    yrk4[i,:] = solver(testode, yrk4[i-1,:], p, x[i-1], dx)
    ybdf1[i,:] = bdfsolve(testode, ybdf1[i-1,:], p, x[i-1], dx)
    yana[i,:] = analytical(x[i])
end


prob = ODEProblem(testode, y0, (x[1], x[end]))
sol = solve(prob)
yDE = Array(sol)'

plt = plot(x, yrk4[:,1], lab="y1 - RK4", leg=:topleft)
plot!(x, yrk4[:,2], lab="y2 - RK4")
plot!(x, ybdf1[:,1], lab="y1 - BDF1")
plot!(x, ybdf1[:,2], lab="y2 - BDF1")
plot!(x, yana[:,1], lab="y1 - ana")
plot!(x, yana[:,2], lab="y2 - ana")
plot!(sol.t, yDE[:,1], lab="y1 - DE")
plot!(sol.t, yDE[:,2], lab="y2 - DE")
# display(plt)


################# Test the DAE Solvers

function fdae(out,du,u,p,t)
    out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
    out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end

u0 = [1.0, 0, 0]
du0 = [-0.04, 0.04, 0.0]
tspan = (0.0,10.0)

differential_vars = [true,true,false]
prob = DAEProblem(fdae, du0, u0, tspan, differential_vars=differential_vars)


sol = solve(prob,IDA())

tvecdae = sol.t 
nt = length(tvecdae)
dsolver1 = DBDF1!(differential_vars)

u1 = zeros(nt, dsolver1.n)
du1 = zeros(nt, dsolver1.nd)

u1[1,:] = u0
du1[1,:] = du0[1:2]
outs = zero(u0)
p = nothing

for i = 2:nt
    t = tvecdae[i-1]
    dt = tvecdae[i]-tvecdae[i-1]
    u1[i,:], du1[i,:] = dsolver1(fdae, outs, u1[i-1,:], du1[i-1,:], p, t, dt)
end

dsolver2 = DBDF2!(differential_vars)
u2 = zeros(nt, dsolver2.n)
du2 = zeros(nt, dsolver2.nd)

u2[1:2,:] = u1[1:2,:]
du2[1:2,:] = du1[1:2,:]

for i = 3:nt
    t = tvecdae[i-1]
    dt = tvecdae[i]-tvecdae[i-1]
    u2[i,:], du2[i,:] = dsolver2(fdae, outs, u2[i-2:i-1,:], du2[i-1,:], p, t, dt)
end

dplt = plot(xaxis="Time (s)", yaxis="State")
plot!(sol)
plot!(tvecdae, u1; linestyle=:dash, lab=["u1 BDF1" "u2 BDF1" "u3 BDF1"])
plot!(tvecdae, u2; linestyle=:dashdot, lab=["u1 BDF2" "u2 BDF2" "u3 BDF2"])
display(dplt)

nothing