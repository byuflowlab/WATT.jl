#=

=#

import CCBlade.afeval

function afeval(af::DS.Airfoil, alpha, Re, Mach)
    return af.cl(alpha), af.cd(alpha)
end

function update_BEM_variables!(xv, blade, airfoil, env, r, twist, Vx, Vy, pitch)
    # [r, airfoil.c, twist, blade.rhub, blade.rtip, Vx, Vy, env.rho, pitch, env.mu, env.a]
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
function solve_BEM(rotor::Rotor, blade::Blade, env::Environment, idx, Vx, Vy, pitch; npts::Int=10)
    xv = zeros(11)

    return solve_BEM!(rotor, blade, env, idx, Vx, Vy, pitch, xv; npts)
end

function solve_BEM!(rotor::Rotor, blade::Blade, env::Environment, idx, Vx, Vy, pitch, xv; npts::Int=10)

    airfoil = blade.airfoils[idx]
    rR = blade.rR[idx]
    # r = rR*(blade.rtip-blade.rhub)
    # r = sqrt(blade.rx[idx]^2 + blade.ry[idx]^2 + blade.rz[idx]^2)
    r = blade.r[idx]
    twist = blade.twist[idx]

    # check if we are at hub/tip
    if isapprox(rR, 0.0, atol=1e-6) || isapprox(rR, 1.0, atol=1e-6)
        return Outputs()  # no loads at hub/tip
    end

    ### unpack
    # Vx = op.Vx
    # Vy = op.Vy
    # theta = section.theta + op.pitch
    theta = twist + pitch

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

    # package up variables and parameters for residual
    # xv = [r, airfoil.c, twist, blade.rhub, blade.rtip, Vx, Vy, env.rho, pitch, env.mu, env.a] #Todo: I wonder how I can make this vector inplace, or limit allocations here. -> Just pass in a dummy vector that get's replaced every function call.  
    update_BEM_variables!(xv, blade, airfoil, env, r, twist, Vx, Vy, pitch)
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
        success, phiL, phiU = CCBlade.firstbracket(phi -> residual(phi, xv, pv), phimin, phimax, npts, backwardsearch)

        function solve(x, p)
            phistar, _ = FLOWMath.brent(phi -> residual(phi, x, p), phiL, phiU)
            return phistar
        end

        # once bracket is found, solve root finding problem and compute loads
        if success
            phistar = IAD.implicit(solve, residual, xv, pv)
            _, outputs = CCBlade.residual_and_outputs(phistar, xv, pv)
            return outputs
        end    
    end    

    # it shouldn't get to this point.  if it does it means no solution was found
    # it will return empty outputs
    # alternatively, one could increase npts and try again
    
    @warn "Invalid data (likely) for this section.  Zero loading assumed."
    return CCBlade.Outputs()
end




const InductionLimit = 1000000.0


# function solve_inflowangle_aerodyn(U, Omega, r, rtip, rhub, B, c, beta, clfit, cdfit, maxiters, tolerance)

    
#     lambda = Omega*r/U #Local tip speed ratio
#     Vx = U
#     Vy = Omega*r
#     sigma = B*c/(2*pi*r) #local solidity #x

#     # println("sigma: ", sigma)


#     # Estimate the axial induction factor
#     sqrt_term = sqrt(4 - 4*pi*lambda*sigma + pi*(lambda^2)*sigma*(8*beta + pi*sigma))
#     a = (2 + pi*lambda*sigma - sqrt_term)/4 #Estimated axial induction factor, EQ 20
#     ap = 0.0 #Assume zero an initial tangential induction (just befor EQ 21)
#     phi = 0.0

#     # phi = 0.0001304168554257867 #It converges to the same solution, even when it's given the correct solution. 
#     # a = 0.9973918231563605
#     # ap = -0.09984575201529174

#     # println("a initial: ", a)

#     aold = a
#     apold = ap
#     phiold = 0.0
#     R = 1.0
#     converged = false

#     for i = 1:maxiters
#         # Estimate the inflow angle
#         phi = atan(Vx*(1-a), Vy*(1+ap)) #Inflow angle, more or less EQ 21
#         # println("phi: ", phi)

#         cphi = cos(phi)
#         sphi = sin(phi)

#         #Calculate the thrust based on inflow angle
#         ## Look up the lift and drag
#         alpha = phi - beta #Angle of attack, in the paragraph before EQ 20
#         cl = clfit(alpha) #Static section coefficient of lift
#         cd = cdfit(alpha) #Static section coefficient of drag

#         # Ct_top = sigma*((1-a)^2)*(cl*cphi + cd*sphi) # A portion of Equation 22
#         # Ct = 1 + Ct_top/(sphi^2) #The thrust coefficient, EQ 22



#         #Calculate the correction factors
#         # ftip = B*(rtip-r)/(2*r*sin(phi)) #Part of equation 23
#         # Ftip = 2*acos(exp(-ftip)) #Tip correction factor, EQ 23
#         # fhub = B*(r-rhub)/(2*r*sin(phi)) #Part of equation 24
#         # Fhub = 2*acos(exp(-fhub)) #Prandtl hub correction factor, EQ 24
#         # F = Ftip*Fhub #Correction factor, EQ 25
#         F = 1.0

#         ## Calculate the induction factors
#         # if Ct>0.96*F #Use the modified Glauert correction
#         #     atop = 18*F - 20 - 3*sqrt(Ct*(50-36*F) + 12*F*(3*F - 4))
#         #     abot = 36*F - 50
#         #     a = atop/abot #Axial induction factor with the modified Glauert correction, EQ 26
#         # else #Standard BEM theory
#         #     atop = 4*F*(sphi^2)
#         #     abot = sigma*(cl*cphi + cd*sphi)
#         #     a = 1/(1 + atop/abot)
#         # end

#         # aptop = 4*F*sphi*cphi
#         # apbot = sigma*(cl*sphi - cd*cphi)
#         # ap = 1/(aptop/apbot - 1)

#         # if skew
#         #     a = a*(1 + 15*pi*r *tan(chi/2)*cos(psi)/(32*rtip)) #apply the skewed wake correction factor, EQ 29
#         # end

#         ### What AeroDyn (BEMTUncoupled.f90) has in their code. 
#         cx = cl*cphi + cd*sphi
#         cy = cl*sphi - cd*cphi

#         a, ap, R = inductionfactors(phi, U, c, r, B, cx, cy, F, lambda)
#         # @show R
        
#         #Check convergence
#         aerr = abs(a-aold)
#         aperr = abs(ap-apold)
#         phierr = abs(phi-phiold)

#         # @show aperr

#         # if (aerr<tolerance)&&(aperr<tolerance)&&(phierr<tolerance)&&(R<tolerance)
#         if (aerr<tolerance)&&(aperr<tolerance)&&(phierr<tolerance)
#             converged = true
#             println("Converged after $i iterations")
#             # @show R
#             # @show ap
#             break
#         end

#         # @show phi, a, ap, Ct
#         println("")
        

#         aold = a
#         apold = ap
#         phiold = phi
#     end
#     return phi, a, ap, converged
# end

"""
signf(x, y)

Fortran's sign transfer function. 
"""
function signf(x, y)
    if y>=0
        return abs(x)
    end
    return -abs(x)
end

function inductionfactors(phi, U, omega, c, r, B, cn, ct, F)
    sphi = sin(phi)
    cphi = cos(phi)


    sigma = B*c/(2*pi*r) #x
    k = sigma*cn/(4*F*(sphi^2)) #x

    lambda = omega*r/U

    momentum_region = (phi>0 && U>=0)||(phi<0 && U<0)

    if momentum_region #Momentum Region
        if k < 2/3
            if isapprox(k, -1)
                a = -signf(InductionLimit, 1+k)
            else
                a = k/(1+k)
            end
            
            if k<-1
                # error("Invalid solution")
                @warn("Invalid solution")
                a = 0.0
            end
        else #Glauert (Buhl) Correction
            temp = 2*F*k
            g1 = temp - (10/9 - F)
            g2 = temp - (4/3 - F)*F
            g3 = temp - (25/9 -2*F)

            if abs(g3)<1e-6
                a = 1 - sqrt(g2)/2
            else
                a = (g1 - sqrt(g2))/g3
            end
        end
    else #Propeller Brake
        if isapprox(k, 1)
            @warn("Invalid Solution")
            a = InductionLimit
        else
            a = k/(k-1)
        end

        if k<=1
            # error("Invalid Solution: Propeller Brake: k<=1")
            @warn("Invalid Solution: Propeller Brake: k<=1")
            # a = InductionLimit
            a = 0.0
        end
    end

    if isapprox(cphi, 0)
        ap = -1
        kp = signf(InductionLimit, ct*sphi)*signf(1, U)
    else
        kp = sigma*ct/(4*F*sphi*cphi)
        if U<0
            kp *= -1
        end
        if isapprox(kp, 1)
            ap = signf(InductionLimit, 1-kp)
        else
            ap = kp/(1-kp)
        end
    end

    if momentum_region
        if isapprox(a, 1)
            R = -cphi*(1-kp)/lambda
        else
            R = sphi/(1-a) -cphi*(1-kp)/lambda
        end
    else
        R = sphi*(1-k) - cphi*(1-kp)/lambda
    end

    return a, ap, R
end