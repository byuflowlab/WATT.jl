abstract type Solver end


struct RK4 <: Solver 
end

function (s::RK4)(fun, x, p, t, dt)
    k1 = fun(x, p, t)
    k2 = fun(x .+ (dt/2).*k1, p, t + (dt/2))
    k3 = fun(x .+ (dt/2).*k2, p, t + (dt/2))
    k4 = fun(x .+ dt.*k3, p, t + dt)

    return @. x + (k1 + k2*2 + k3*2 + k4)*dt/6
end