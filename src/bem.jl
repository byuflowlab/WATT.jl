const InductionLimit = 1000000.0


function solve_inflowangle_aerodyn(U, Omega, r, rtip, rhub, B, c, beta, clfit, cdfit, maxiters, tolerance)

    
    lambda = Omega*r/U #Local tip speed ratio
    Vx = U
    Vy = Omega*r
    sigma = B*c/(2*pi*r) #local solidity #x

    # println("sigma: ", sigma)


    # Estimate the axial induction factor
    sqrt_term = sqrt(4 - 4*pi*lambda*sigma + pi*(lambda^2)*sigma*(8*beta + pi*sigma))
    a = (2 + pi*lambda*sigma - sqrt_term)/4 #Estimated axial induction factor, EQ 20
    ap = 0.0 #Assume zero an initial tangential induction (just befor EQ 21)
    phi = 0.0

    # phi = 0.0001304168554257867 #It converges to the same solution, even when it's given the correct solution. 
    # a = 0.9973918231563605
    # ap = -0.09984575201529174

    # println("a initial: ", a)

    aold = a
    apold = ap
    phiold = 0.0
    R = 1.0
    converged = false

    for i = 1:maxiters
        # Estimate the inflow angle
        phi = atan(Vx*(1-a), Vy*(1+ap)) #Inflow angle, more or less EQ 21
        # println("phi: ", phi)

        cphi = cos(phi)
        sphi = sin(phi)

        #Calculate the thrust based on inflow angle
        ## Look up the lift and drag
        alpha = phi - beta #Angle of attack, in the paragraph before EQ 20
        cl = clfit(alpha) #Static section coefficient of lift
        cd = cdfit(alpha) #Static section coefficient of drag

        # Ct_top = sigma*((1-a)^2)*(cl*cphi + cd*sphi) # A portion of Equation 22
        # Ct = 1 + Ct_top/(sphi^2) #The thrust coefficient, EQ 22



        #Calculate the correction factors
        # ftip = B*(rtip-r)/(2*r*sin(phi)) #Part of equation 23
        # Ftip = 2*acos(exp(-ftip)) #Tip correction factor, EQ 23
        # fhub = B*(r-rhub)/(2*r*sin(phi)) #Part of equation 24
        # Fhub = 2*acos(exp(-fhub)) #Prandtl hub correction factor, EQ 24
        # F = Ftip*Fhub #Correction factor, EQ 25
        F = 1.0

        ## Calculate the induction factors
        # if Ct>0.96*F #Use the modified Glauert correction
        #     atop = 18*F - 20 - 3*sqrt(Ct*(50-36*F) + 12*F*(3*F - 4))
        #     abot = 36*F - 50
        #     a = atop/abot #Axial induction factor with the modified Glauert correction, EQ 26
        # else #Standard BEM theory
        #     atop = 4*F*(sphi^2)
        #     abot = sigma*(cl*cphi + cd*sphi)
        #     a = 1/(1 + atop/abot)
        # end

        # aptop = 4*F*sphi*cphi
        # apbot = sigma*(cl*sphi - cd*cphi)
        # ap = 1/(aptop/apbot - 1)

        # if skew
        #     a = a*(1 + 15*pi*r *tan(chi/2)*cos(psi)/(32*rtip)) #apply the skewed wake correction factor, EQ 29
        # end

        ### What AeroDyn (BEMTUncoupled.f90) has in their code. 
        cx = cl*cphi + cd*sphi
        cy = cl*sphi - cd*cphi

        a, ap, R = inductionfactors(phi, U, c, r, B, cx, cy, F, lambda)
        # @show R
        
        #Check convergence
        aerr = abs(a-aold)
        aperr = abs(ap-apold)
        phierr = abs(phi-phiold)

        # @show aperr

        # if (aerr<tolerance)&&(aperr<tolerance)&&(phierr<tolerance)&&(R<tolerance)
        if (aerr<tolerance)&&(aperr<tolerance)&&(phierr<tolerance)
            converged = true
            println("Converged after $i iterations")
            # @show R
            # @show ap
            break
        end

        # @show phi, a, ap, Ct
        println("")
        

        aold = a
        apold = ap
        phiold = phi
    end
    return phi, a, ap, converged
end

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

    # cphi:   0.999999991536066     
    # sphi:   1.301071420651440E-004

    # println("phi: ", phi)
    # println("sphi: ", sphi)
    # println("cphi: ", cphi)

    sigma = B*c/(2*pi*r) #x
    k = sigma*cn/(4*F*(sphi^2)) #x #73595.0332412394

    lambda = omega*r/U


    # println("cn: ", cn)
    # println("F: ", F)
    # println("sigma: ", sigma)
    # println("k: ", k) 

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