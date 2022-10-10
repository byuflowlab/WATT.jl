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

struct DiffEQ <: Solver
    algorithm
end

#The idea behind this solver is to take advantage of the integrator/step approach and use a callback to update p. 
function (s::DiffEQ)(fun, x, p, t, dt) #I think the integrator keeps track of x, t and what not. 
    fun.p .= p
    DifferentialEquations.step!(fun, dt, true)
    return fun.u
end

struct DiffEQinit <: Solver #TODO: This could probably be implemented a lot better. 
end

function (s::DiffEQinit)(fun, x, p, t, dt) #Takes about the same amount of time as the BDF1. 
    prob = ODEProblem(fun, x, (t, t+dt), p)
    sol = DifferentialEquations.solve(prob)
    return sol(t+dt)
end

# import DS.Indicial

# function (s::DS.Indicial)(fun, x, p, t, dt)
# end









struct DBDF1!{TI} <: Solver
    diffvars::Array{TI, 1} #which variables are states
    constvars::Array{TI, 1} #which variables are constraints
    n::TI #Number of states and constaints
    nd::TI #Number of states
end

function DBDF1!(diffvars::Array{Bool, 1})
    n = length(diffvars)
    nd = 0
    diffs = Int[] #TODO: There's probably a better way to do this. -> I could probably iterate through and get nd, then create diffs and consts
    consts = Int[]
    for i = 1:n
        if diffvars[i]
            nd += 1
            push!(diffs, i)
        else
            push!(consts, i)
        end
    end
    return DBDF1!(diffs, consts, n, nd)
end


function (s::DBDF1!)(fun!, outs, x, dx, p, t, dt)
    tn = t + dt

    residfun! = function(resids, rn)  
        
        xn = view(rn, 1:s.n)
        dxn = view(rn, s.n+1:s.n+s.nd) #Should be length s.nd (length of the diffvars)

        fun!(outs, dxn, xn, p, tn)

        resids[1:s.n] = outs #Converge the derivatives and constraints

        for i in 1:s.nd #Converge the derivative/state relationships
            diffidx = s.diffvars[i]
            resids[diffidx] = x[diffidx] + dt*dxn[i] - xn[diffidx] #Insure that the new states satisfy backwards compatibility. 
        end
    end

    residuals = nlsolve(residfun!, vcat(x, dx)) 

    # @show residuals.zero[s.n+1:end] #The state rates are zero... which is problematic. 

    return residuals.zero[1:s.n], residuals.zero[s.n+1:end]
end

struct DBDF2!{TI} <: Solver
    diffvars::Array{TI, 1} #which variables are states
    constvars::Array{TI, 1} #which variables are constraints
    n::TI #Number of states and constaints
    nd::TI #Number of states
end

function DBDF2!(diffvars::Array{Bool, 1})
    n = length(diffvars)
    nd = 0
    diffs = Int[] #TODO: There's probably a better way to do this. -> I could probably iterate through and get nd, then create diffs and consts
    consts = Int[]
    for i = 1:n
        if diffvars[i]
            nd += 1
            push!(diffs, i)
        else
            push!(consts, i)
        end
    end
    
    return DBDF2!(diffs, consts, n, nd)
end

function (s::DBDF2!)(fun!, outs, x, dx, p, t, dt)
    tn = t + dt

    #Coefficients for the method
    a = 4/3
    b = 2*dt/3
    c = 1/3
    
    residfun! = function(resids, rn)  
        
        xn = view(rn, 1:s.n)
        dxn = view(rn, s.n+1:s.n+s.nd) #Should be length s.nd (length of the diffvars)

        fun!(outs, dxn, xn, p, tn)

        resids[1:s.n] = outs #Converge the derivatives and constraints

        for i in 1:s.nd #Converge the derivative/state relationships
            diffidx = s.diffvars[i]
            # resids[diffidx] = x[diffidx] + dt*dxn[i] - xn[diffidx] #Insure that the new states satisfy backwards compatibility. 
            resids[diffidx] = a*x[2, diffidx] + b*dxn[i] - c*x[1, diffidx] - xn[diffidx]
        end
    end

    residuals = nlsolve(residfun!, vcat(x[2,:], dx)) 

    return residuals.zero[1:s.n], residuals.zero[s.n+1:end]
end


function fixedpointbem(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, phi0; atol=1e-8, maxiter::Int=500)

    phistar = phi0 
    converged = false
    iter = 0
    while !converged
        # println(phistar)
        r, out = CCBlade.residual(phistar, rotor, section, op) 
        # println(r)

        if rotor.turbine
            phistar = atan(op.Vx*(1-out.a), op.Vy*(1+out.ap))
        else
            phistar = atan(op.Vx*(1+out.a), op.Vy*(1-out.ap))
        end
        # println(phistar)
        iter += 1

        if r<=atol
            converged = true
        elseif iter>=maxiter #I'm not convinced that this is converging. Can it not hit the level of tolerance? 
            # @warn("Maximum number of iterations hit.") #It's hitting 500 iterations... like alot. 
            break
        end

        # phistar = phinew

    end
    _, out = CCBlade.residual(phistar, rotor, section, op)
    return out
end


function secant(fun, x0, x1; atol=1e-8, rtol=4*eps(), maxiter=100)

    # if isapprox(x0, x1; atol) #Todo: The method will break if the inputs are the same... which... I wanted to put the two previous inflow angles in... which would break this. 
    #     error("secant() can't have x0 and x1 be the same argument.")
    # end

    r0 = fun(x0)
    r1 = fun(x1)
    converged = false
    iter = 1
    while !converged

        top = r1*(x0-x1)
        bot = r0 - r1
        x2 = x1 - top/bot

        r2 = fun(x2)
        if r2<atol
            converged = true
        elseif abs(r2-r1)<rtol
            converged = true
        end

        r0 = r1
        r1 = r2

        x0 = x1
        x1 = x2
        iter += 1
    end

    return x2
end