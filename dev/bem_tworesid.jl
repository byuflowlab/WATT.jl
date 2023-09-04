#=
Writing a two residual BEM to converge for comparison with AeroDyn. 

Adam Cardoza 3/17/23

The following equations are from the AeroDyn Theory documentation

1. Estimate the axial induction factor.
2. Determine the thrust coefficient, Ct. 
3. Calculate the corrections, F.
4. Calculate the new axial induction factor, and the tangential induction factor.
5. Iterate between steps 2, 3, and 4 until a, a', and phi converge. 

a. There must be a calculate the flow angle step in there because the equations for the induction factors are functions of phi. 
-> tan(phi) = U*(1-a)/(Omega*r*(1+a'))
=#

function calcPhi(a, ap, U, lambda)
    top = U*(1-a)
    bot = lambda*(1+ap)
    return atan(top, bot)
end

function 

function solveBEM(U, Omega, r, maxiters, tolerance, clfit, cdfit)

    # Estimate the axial induction factor
    lambda = Omega*r
    sqrt_term = sqrt(4 - 4*pi*lambda*sigma + pi*(lambda^2)*sigma*(8*beta + pi*sigma))
    a = (2 + pi*lambda*sigma - sqrt_term)/4 #Estimated axial induction factor, EQ 20
    ap = 0.0 #Assume zero an initial tangential induction (just befor EQ 21)

    aold = a
    apold = ap
    phiold = 0.0
    

    for i = 1:maxiters
        # Estimate the inflow angle
        phi = calcPhi(a, ap, U, lambda) #Inflow angle, more or less EQ 21

        #Calculate the thrust based on inflow angle
        ## Look up the lift and drag
        alpha = phi - beta #Angle of attack, in the paragraph before EQ 20
        cl = clfit(alpha) #Static section coefficient of lift
        cd = cdfit(alpha) #Static section coefficient of drag

        Ct_top = sigma*((1-a)^2)*(cl*cos(phi) + cd*sin(phi))
        Ct = 1 + Ct_top/(sin(phi)^2)

        #Calculate the correction factors

        # Calculate the induction factors

        #Check convergence
        aerr = abs(a-aold)
        aperr = abs(ap-apold)
        phierr = abs(phi-phiold)

        if (aerr<tolerance)&&(aperr<tolerance)&&(phierr<tolerance)
            break
        end
    end
    return phi, a, ap
end


function bemresiduals()


end
