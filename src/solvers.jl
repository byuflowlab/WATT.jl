#=
The purpose of this file is to provide some simple ODE solution methods. I used my own custom solvers over DifferentialEquations so that I'd be fairly confident in what was happening. I wanted to be able to control taking a single step a time, and be able to see what was happening during the integration. This was much easier in my own functions rather than looking inside DifferentialEquations. 

The plan for the future is either figure out how to use DifferentialEquations and replace or combine it with these methods. If not, then I might add some better integration methods. Maybe a Tsit method or BDF? I dunno. Kind of depends on how much time I have. Honestly, this is probably where the solvers will stay.
- Adam Cardoza, 7/20/22

=#


"""
    Solver

The abstract type of all the solvers. 
"""
abstract type Solver end


"""
    RK4()

A struct to leverage multiple dispatch on what solver the user would like to use. There are no fields. 
"""
struct RK4 <: Solver 
end

"""
    solver(fun, x, p, t, dt)

A method on the RK4 struct. This method takes a single timestep dt from the function evaluated at x, p, t using the Runge-Kutta fourth order method. 

### Inputs
- fun::Function - a function to return the state rates dx for a given set of states x. The arguments of this function should be fun(x, p, t)
- x::Union{TF, Array{TF,1}} - either a state or a vector of states, depending on fun. 
- p::Any - a vector to pass parameters into the ode state rates function. 
- t::TF - the current time value. 
- dt::TF - the desired time step. 

### Outputs
- x_new::Union{TF, Array{TF,1}} - the updated state(s). 
"""
function (s::RK4)(fun, x, p, t, dt)
    k1 = fun(x, p, t)
    k2 = fun(x .+ (dt/2).*k1, p, t + (dt/2))
    k3 = fun(x .+ (dt/2).*k2, p, t + (dt/2))
    k4 = fun(x .+ dt.*k3, p, t + dt)

    return @. x + (k1 + k2*2 + k3*2 + k4)*dt/6
end

struct BDF1 <: Solver
end

function (s::BDF1)(fun, x, p, t, dt) #Takes about 2x the time of the RK4
    tn = t + dt
    
    fun! = function(residual, xn)
        residual .= x .+ dt.*fun(xn, p, tn) - xn
    end

    xn = nlsolve(fun!, x) #TODO: Previously I had a initial guess that I initialized with an RK4 for the nonlinear solve. I'm using the current state right now, but I might try using the RK4 later. 
    return xn.zero
end

struct DiffEQ <: Solver #TODO: This could probably be implemented a lot better. 
end

function (s::DiffEQ)(fun, x, p, t, dt) #Takes about the same amount of time as the BDF1. 
    prob = ODEProblem(fun, x, (t, t+dt), p)
    sol = DifferentialEquations.solve(prob)
    return sol(t+dt)
end