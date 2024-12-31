#=

The blade element momentum theory (BEMT) code. Most of this wraps the package CCBlade.jl

TODO: I should probably rename this file to BEMT. It appears that most people refer to panel codes as BEMs, and refer to BEMT as BEMT. 

=#

import CCBlade.afeval

function afeval(af::DS.Airfoil, alpha, Re, Mach)
    return af.cl(alpha), af.cd(alpha)
end

function update_BEM_variables!(xv, blade, airfoil, env, r, twist, Vx, Vy, pitch)
    # [r, airfoil.c, twist, blade.rhub, blade.rtip, Vx, Vy, env.rho, pitch, env.mu, env.a]
    # @show typeof(twist) #Shows as a tracked real when it needs to.
    xv[1] = r
    xv[2] = airfoil.c
    xv[3] = twist
    xv[4] = blade.rhub
    xv[5] = blade.rtip
    xv[6] = Vx
    xv[7] = Vy
    xv[8] = env.rho
    xv[9] = pitch
    xv[10] = env.mu
    xv[11] = env.a
end



"""
    solve_BEM(rotor, blade, env, idx, Vx, Vy, pitch, op; npts=10)
Solve the BEM equations for given rotor geometry and operating point. This function is based on CCBlade's solve function. 

### Inputs
- `rotor::Rotor` - rotor properties
- `blade::Blade` - section properties
- `env::Environment` - operating point
- `idx::Int` - The 
- `npts::Int` - number of discretization points to find bracket in residual solve

### Outputs
- `outputs::Outputs`: BEM output data including loads, induction factors, etc.
"""
function solve_BEM(rotor::Rotor, blade::Blade, env::Environment, phi0, idx, Vx, Vy, pitch; npts::Int=10, newbounds::Bool=true)
    xv = zeros(11)

    return solve_BEM!(rotor, blade, env, phi0, idx, Vx, Vy, pitch, xv; npts, newbounds)
end

function solve_BEM!(rotor::Rotor, blade::Blade, env::Environment, phi0, idx, Vx, Vy, pitch, xv; twist=blade.twist[idx], npts::Int=10, newbounds::Bool=true)
#TODO: This shouldn't have an exclamation point in the name because it doesn't actually appear to be inplace. And it returns a struct, not a vector, so it's not an allocation... Or shouldn't be. 


    airfoil = blade.airfoils[idx]
    # rR = blade.rR[idx]
    r = blade.r[idx]
    

    # check if we are at hub/tip
    # if isapprox(rR, 0.0, atol=1e-6) || isapprox(rR, 1.0, atol=1e-6)
    #     return Outputs()  # no loads at hub/tip
    # end

    ### unpack
    theta = twist + pitch

    # package up variables and parameters for residual 
    update_BEM_variables!(xv, blade, airfoil, env, r, twist, Vx, Vy, pitch) #TODO: I actually might be able to do this with a tuple, because I don't think a tuple allocates anything. 
    pv = (airfoil, rotor.B, rotor.turbine, rotor.re, rotor.mach, rotor.rotation, rotor.tip)

    if newbounds
        res0, outputs0 = CCBlade.residual_and_outputs(phi0, xv, pv)

        if isapprox(res0, 0.0, atol=1e-6) #TODO: I might consider raising that tolerance... so I'm not solving the BEMT every time. (FLOWMath has a tolerance of 2e-12)
            # println("BEMT: Using previous solution...")
            return outputs0
        end
    end

    # ---- determine quadrants based on case -----
    Vx_is_zero = isapprox(Vx, 0.0, atol=1e-6)
    Vy_is_zero = isapprox(Vy, 0.0, atol=1e-6)

    ### quadrants
    epsilon = 1e-6
    q1 = [epsilon, pi/2]
    q2 = [-pi/2, -epsilon]
    q3 = [pi/2, pi-epsilon]
    q4 = [-pi+epsilon, -pi/2]

    if Vx_is_zero && Vy_is_zero
        return CCBlade.Outputs()

    elseif Vx_is_zero

        startfrom90 = false  # start bracket at 0 deg.

        if Vy > 0 && theta > 0
            order = (q1, q2)
        elseif Vy > 0 && theta < 0
            order = (q2, q1)
        elseif Vy < 0 && theta > 0
            order = (q3, q4)
        else  # Vy < 0 && theta < 0
            order = (q4, q3)
        end

    elseif Vy_is_zero

        startfrom90 = true  # start bracket search from 90 deg

        if Vx > 0 && abs(theta) < pi/2
            order = (q1, q3)
        elseif Vx < 0 && abs(theta) < pi/2
            order = (q2, q4)
        elseif Vx > 0 && abs(theta) > pi/2
            order = (q3, q1)
        else  # Vx < 0 && abs(theta) > pi/2
            order = (q4, q2)
        end

    else  # normal case

        startfrom90 = false

        if Vx > 0 && Vy > 0
            order = (q1, q2, q3, q4)
        elseif Vx < 0 && Vy > 0
            order = (q2, q1, q4, q3)
        elseif Vx > 0 && Vy < 0
            order = (q3, q4, q1, q2)
        else  # Vx[i] < 0 && Vy[i] < 0
            order = (q4, q3, q2, q1)
        end

    end

    # ----- solve residual function ------
    # pull out first argument
    residual(phi, x, p) = CCBlade.residual_and_outputs(phi, x, p)[1]

    

    success = false
    for j = 1:length(order)  # quadrant orders.  In most cases it should find root in first quadrant searched.
        phimin, phimax = order[j]

        # check to see if it would be faster to reverse the bracket search direction
        backwardsearch = false
        if !startfrom90
            if phimin == -pi/2 || phimax == -pi/2  # q2 or q4
                backwardsearch = true
            end
        else
            if phimax == pi/2  # q1
                backwardsearch = true
            end
        end

        # find bracket
        success, phiL, phiU = CCBlade.firstbracket(phi -> residual(phi, xv, pv), phimin, phimax, npts, backwardsearch)

        

        # once bracket is found, solve root finding problem and compute loads
        if success
            if newbounds
                resL = residual(phiL, xv, pv)
                resU = residual(phiU, xv, pv)
    
                if phiL<=phi0<=phiU
                    if resL*res0<0
                        phiU = phi0
                        # println("replacing upper bound")
                    elseif resU*res0<0
                        phiL = phi0
                        # println("replacing lower bound")
                    end
                else
                    # println("Not using bounds...")
                end
            end
    
            function solve(x, p) #TODO: Is there a more efficient way to do this instead of a closure? -> Is this a problem? If phiL and phiU were vectors, I think it might be a problem, but since they are scalars, I don't think it is.
                phistar, _ = FLOWMath.brent(phi -> residual(phi, x, p), phiL, phiU)
                return phistar
            end

            phistar = IAD.implicit(solve, residual, xv, pv)
            @show typeof(phistar)
            _, outputs = CCBlade.residual_and_outputs(phistar, xv, pv) #TODO: Instead of creating a new function... I could just check if x0 is inbetween a and b outside of the loop. If it is, then I can check if it is a zero. If not, then replace one of the bounds into Brent's method. -> Is this even a problem? 
            return outputs
        end    
    end    

    # it shouldn't get to this point.  if it does it means no solution was found
    # it will return empty outputs
    # alternatively, one could increase npts and try again
    
    @warn "Invalid data (likely) for this section.  Zero loading assumed."
    return CCBlade.Outputs()
end

function show_dual_vec(x::AbstractArray{<:ForwardDiff.Dual})
    for i in eachindex(x)
        println(x[i].value)
    end
end

"""
    solve_BEM!(rotor, blade, env, idx, Vx, Vy, pitch, xv; twist=blade.twist[idx], npts::Int=10, epsilon=1.1e-3)

Solve the BEMT equations for given rotor geometry and operating point. This function is based on CCBlade's solve function.

**Arguments**
- `rotor::Rotor`: rotor object
- `blade::Blade`: blade object
- `env::Environment`: environment object
- `idx::Int`: index of the blade section
- `Vx::Float64`: x-component of the wind velocity (Freestream)
- `Vy::Float64`: y-component of the wind velocity (Rotational)
- `pitch::Float64`: pitch angle (radians)
- `xv::Vector{Float64}`: vector of variables for the residual function
- `twist::Float64`: twist angle (radians)
- `npts::Int`: number of discretization points to find bracket in residual solve
- `epsilon::Float64`: small value to avoid division by zero

**Returns**
- `outputs::Outputs`: BEMT output data including inflow angle, angle of attack, loads, induction factors, etc.
"""
function solve_BEM!(rotor::Rotor, blade::Blade, env::Environment, idx, Vx, Vy, pitch, xv; twist=blade.twist[idx], npts::Int=10, epsilon=1e-6) #epsilon=1.1e-3)
    #TODO: This shouldn't have an exclamation point in the name because it doesn't actually appear to be inplace. And it returns a struct, not a vector, so it's not an allocation... Or shouldn't be. 
    
    
    airfoil = blade.airfoils[idx]
    # rR = blade.rR[idx]
    r = blade.r[idx]
    

    # check if we are at hub/tip
    # if isapprox(rR, 0.0, atol=1e-6) || isapprox(rR, 1.0, atol=1e-6)
    #     return Outputs()  # no loads at hub/tip
    # end

    ### unpack
    theta = twist + pitch



    # ---- determine quadrants based on case -----
    Vx_is_zero = isapprox(Vx, 0.0, atol=1e-6)
    Vy_is_zero = isapprox(Vy, 0.0, atol=1e-6)

    phi_geo = atan(Vx, Vy)
    if abs(phi_geo)<=epsilon
        warn("The geometric inflow angle is approximately the same as the solution ϵ, which may cause the BEMT to fail.")
    end

    ### quadrants
    # epsilon = 1e-6
    q1 = [epsilon, pi/2]
    q2 = [-pi/2, -epsilon]
    q3 = [pi/2, pi-epsilon]
    q4 = [-pi+epsilon, -pi/2]

    if Vx_is_zero && Vy_is_zero
        println("Vx and Vy is zero.")
        return CCBlade.Outputs()

    elseif Vx_is_zero

        startfrom90 = false  # start bracket at 0 deg.

        if Vy > 0 && theta > 0
            order = (q1, q2)
        elseif Vy > 0 && theta < 0
            order = (q2, q1)
        elseif Vy < 0 && theta > 0
            order = (q3, q4)
        else  # Vy < 0 && theta < 0
            order = (q4, q3)
        end

    elseif Vy_is_zero

        startfrom90 = true  # start bracket search from 90 deg

        if Vx > 0 && abs(theta) < pi/2
            order = (q1, q3)
        elseif Vx < 0 && abs(theta) < pi/2
            order = (q2, q4)
        elseif Vx > 0 && abs(theta) > pi/2
            order = (q3, q1)
        else  # Vx < 0 && abs(theta) > pi/2
            order = (q4, q2)
        end

    else  # normal case

        startfrom90 = false

        if Vx > 0 && Vy > 0
            order = (q1, q2, q3, q4)
        elseif Vx < 0 && Vy > 0
            order = (q2, q1, q4, q3)
        elseif Vx > 0 && Vy < 0
            order = (q3, q4, q1, q2)
        else  # Vx[i] < 0 && Vy[i] < 0
            order = (q4, q3, q2, q1)
        end

    end

    # ----- solve residual function ------
    # pull out first argument
    residual(phi, x, p) = CCBlade.residual_and_outputs(phi, x, p)[1]

    # package up variables and parameters for residual 
    update_BEM_variables!(xv, blade, airfoil, env, r, twist, Vx, Vy, pitch) #TODO: I actually might be able to do this with a tuple, because I don't think a tuple allocates anything. 
    pv = (airfoil, rotor.B, rotor.turbine, rotor.re, rotor.mach, rotor.rotation, rotor.tip)

    success = false
    for j = 1:length(order)  # quadrant orders.  In most cases it should find root in first quadrant searched.
        phimin, phimax = order[j]

        # check to see if it would be faster to reverse the bracket search direction
        backwardsearch = false
        if !startfrom90
            if phimin == -pi/2 || phimax == -pi/2  # q2 or q4
                backwardsearch = true
            end
        else
            if phimax == pi/2  # q1
                backwardsearch = true
            end
        end

        # find bracket
        # @show xv
        if isa(xv[1], ReverseDiff.TrackedReal)
            @show airfoil.c.value
            println([xv[i].value for i in eachindex(xv)])
        end
        success, phiL, phiU = CCBlade.firstbracket(phi -> residual(phi, xv, pv), phimin, phimax, npts, backwardsearch)

        # once bracket is found, solve root finding problem and compute loads
        if success
            function solve(x, p) #todo: Is there a more efficient way to do this instead of a closure? 
                phistar, _ = FLOWMath.brent(phi -> residual(phi, x, p), phiL, phiU)
                
                return phistar
            end
            
            phistar = IAD.implicit(solve, residual, xv, pv)
            
            _, outputs = CCBlade.residual_and_outputs(phistar, xv, pv) #TODO: Instead of creating a new function... I could just check if x0 is inbetween a and b outside of the loop. If it is, then I can check if it is a zero. If not, then replace one of the bounds into Brent's method. 
                return outputs
            end    
        end    

    # it shouldn't get to this point.  if it does it means no solution was found
    # it will return empty outputs
    # alternatively, one could increase npts and try again
    
    @warn "Invalid data (likely) for this section.  Zero loading assumed."
    return CCBlade.Outputs()
end


"""
    thrusttorque(rotor, sections, outputs::AbstractVector{TO}) where TO

integrate the thrust/torque across the blade, 
including 0 loads at hub/tip, using a trapezoidal rule.

**Arguments**
- `rotor::Rotor`: rotor object
- `sections::Vector{Section}`: rotor object
- `outputs::Vector{Outputs}`: output data along blade

**Returns**
- `T::Float64`: thrust (along x-dir see Documentation).
- `Q::Float64`: torque (along x-dir see Documentation).
"""
function thrusttorque(rvec, N, T, B=1; precone=0.0)
    #Todo: I think I have this function implemented several times across the package. 

    # integrate Thrust and Torque (trapezoidal)
    thrust = N*cos(precone)
    torque = T.*rvec*cos(precone)

    T = B * FLOWMath.trapz(rvec, thrust)
    Q = B * FLOWMath.trapz(rvec, torque)

    return T, Q
end

function ccthrusttorque(rvec, Rhub, Rtip, Np, Tp, B; precone=0.0)

    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    # rvec = [s.r for s in sections]
    rfull = [Rhub; rvec; Rtip]
    Npfull = [0.0; Np; 0.0]
    Tpfull = [0.0; Tp; 0.0]

    # integrate Thrust and Torque (trapezoidal)
    thrust = Npfull*cos(precone)
    torque = Tpfull.*rfull*cos(precone)

    T = B * FLOWMath.trapz(rfull, thrust)
    Q = B * FLOWMath.trapz(rfull, torque)

    return T, Q
end