#=
The purpose of this file is to provide some simple ODE solution methods. I used my own custom solvers over DifferentialEquations so that I'd be fairly confident in what was happening. I wanted to be able to control taking a single step a time, and be able to see what was happening during the integration. This was much easier in my own functions rather than looking inside DifferentialEquations. 

The plan for the future is either figure out how to use DifferentialEquations and replace or combine it with these methods. If not, then I might add some better integration methods. Maybe a Tsit method or BDF? I dunno. Kind of depends on how much time I have. Honestly, this is probably where the solvers will stay.
- Adam Cardoza, 7/20/22

=#
#Todo: I might consider moving all of these methods over to DynamicStallModels. 

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
    (::RK4)(fun, x, p, t, dt) -> x_new

Take a single explicit step of size `dt` using the classical fourth-order
Runge-Kutta method.

**Arguments**
- `fun::Function`: Returns the state rates `dx` for a given state. Called as `fun(x, p, t)`.
- `x::Union{TF, AbstractVector{TF}}`: The current state, or vector of states.
- `p::Any`: Parameters forwarded to `fun`.
- `t::TF`: The current time.
- `dt::TF`: The step size.

**Returns**
- `x_new::Union{TF, AbstractVector{TF}}`: The state advanced to `t + dt`.

**Notes**
Written with broadcast arithmetic and no in-place mutation of `x`, so
ForwardDiff and ReverseDiff propagate through it cleanly. This AD transparency
is the reason the package carries its own integrators instead of calling
DifferentialEquations.
"""
function (s::RK4)(fun, x, p, t, dt)
    k1 = fun(x, p, t)
    k2 = fun(x .+ (dt/2).*k1, p, t + (dt/2))
    k3 = fun(x .+ (dt/2).*k2, p, t + (dt/2))
    k4 = fun(x .+ dt.*k3, p, t + dt)

    return @. x + (k1 + k2*2 + k3*2 + k4)*dt/6
end

"""
    BDF1()

A struct to leverage multiple dispatch on what solver the user would like to use.
There are no fields.

Selects the first-order backward differentiation formula (implicit / backward
Euler). Use it in place of [`RK4`](@ref) when the state rates are stiff enough
that the explicit step goes unstable at the desired `dt`.
"""
struct BDF1 <: Solver
end

"""
    (::BDF1)(fun, x, p, t, dt) -> x_new

Take a single implicit step of size `dt` using the first-order backward
differentiation formula, solving `x_new = x + dt*fun(x_new, p, t+dt)` for
`x_new` with `NLsolve`.

**Arguments**
- `fun::Function`: Returns the state rates `dx` for a given state. Called as `fun(x, p, t)`.
- `x::Union{TF, AbstractVector{TF}}`: The current state, or vector of states.
- `p::Any`: Parameters forwarded to `fun`.
- `t::TF`: The current time.
- `dt::TF`: The step size.

**Returns**
- `x_new::Union{TF, AbstractVector{TF}}`: The state advanced to `t + dt`.

**Notes**
Roughly 2x the cost of an [`RK4`](@ref) step, in exchange for the stability of
an implicit method. The nonlinear solve is seeded with the current state `x`;
seeding it with an explicit RK4 predictor instead would likely cut the iteration
count and has not been tried.

Because the root-find is not wrapped in `ImplicitAD.implicit`, differentiating
through this solver propagates duals through every `NLsolve` iteration rather
than applying the implicit function theorem at the solution. That works, but it
is more expensive than it needs to be.
"""
function (s::BDF1)(fun, x, p, t, dt) #Takes about 2x the time of the RK4
    tn = t + dt
    
    fun! = function(residual, xn)
        residual .= x .+ dt.*fun(xn, p, tn) - xn
    end

    xn = nlsolve(fun!, x) #TODO: Previously I had a initial guess that I initialized with an RK4 for the nonlinear solve. I'm using the current state right now, but I might try using the RK4 later. 
    return xn.zero
end

# struct DiffEQ <: Solver
#     algorithm
# end

# #The idea behind this solver is to take advantage of the integrator/step approach and use a callback to update p. 
# function (s::DiffEQ)(fun, x, p, t, dt) #I think the integrator keeps track of x, t and what not. 
#     fun.p .= p
#     DifferentialEquations.step!(fun, dt, true)
#     return fun.u
# end

# struct DiffEQinit <: Solver #TODO: This could probably be implemented a lot better. 
# end

# function (s::DiffEQinit)(fun, x, p, t, dt) #Takes about the same amount of time as the BDF1. 
#     prob = ODEProblem(fun, x, (t, t+dt), p)
#     sol = DifferentialEquations.solve(prob)
#     return sol(t+dt)
# end

# import DS.Indicial

# function (s::DS.Indicial)(fun, x, p, t, dt)
# end









