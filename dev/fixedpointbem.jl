abstract type StallModel end

struct IndicialRiso
end

struct StateSpaceRiso
end

struct NoStall
end


function afeval(stallmodel::NoStall, alpha, alphadot, w, wdot, state, chord, airfoil) 

    return airfoil.cl(alpha), airfoil.cd(alpha)
end

function afeval(stallmodel::IndicialRiso, alpha, alphadot, w, wdot, state, chord, airfoil) 

    y = [w, wdot, 0.0, 0.0, alpha, alphadot] #TODO: I'm using the quasi-inflow velocity W (the norm of Vx and Vy) rather than the actual inflow angle, or the two velocities seperate... because that's what I've always done... but now I wonder if I'm double counting angles and such.  

    return riso_coefs(state, y, chord, airfoil)
end

function afeval(stallmodel::StateSpaceRiso, alpha, alphadot, w, wdot, state, chord, airfoil) 

    y = [w, wdot, 0.0, 0.0, alpha, alphadot] #TODO: I'm using the quasi-inflow velocity W (the norm of Vx and Vy) rather than the actual inflow angle, or the two velocities seperate... because that's what I've always done... but now I wonder if I'm double counting angles and such.  

    return riso_coefs(state, y, chord, airfoil)
end





abstract type Tipcorrection end #TODO: CamelCase this when I've made this into a package for good style. (Currently avoiding conflict with CCBlade.)

struct Prandtltiphub <: Tipcorrection #TODO: Are there things I can move into this struct? I don't think Dr. Ning did. 
end

function calc_F(tipcorrection::Prandtltiphub, r, rhub, rtip, B, phi) 

    #Hack fix
    # if phi == 0
    #     return 1
    # end


    # Prandtl's tip and hub loss factor
    asphi = abs(sin(phi))

    factortip = B/2.0*(rtip/r - 1)/asphi
    Ftip = 2.0/pi*acos(exp(-factortip))

    factorhub = B/2.0*(r/rhub - 1)/asphi
    Fhub = 2.0/pi*acos(exp(-factorhub))

    F = Ftip * Fhub

    # if isnan(F) #phi is coming in as a NaN
    #     @show r, rhub, rtip, B, phi
    # end
    return F
end







#Todo. Validate this implementation. For not dynamic stall model, with Prandtl Tip/Hub correction, this converges to the CCBlade solution. Boom (after 54 iterations with a relaxation factor of 0.7)
function update_inflowangle(phi, Vx, Vy, state, r, chord, twist, airfoil, pitch, rhub, rtip, B, turbine, tipcorrection, stallmodel)

    ### Constants
    sigma = B*chord/(2*pi*r)
    sphi = sin(phi)
    cphi = cos(phi)

    ### Calculate angle of attack
    alpha = (twist + pitch) - phi

    # if turbine
    #     alpha = -((twist + pitch) - phi)
    # else
    #     error("Propeller Analysis not yet setup.")
    # end
    # @show alpha

    ### CCBlade's lookup. 
        # airfoil cl/cd
        # if rotor.turbine
        #     cl, cd = afeval(af, -alpha, Re, Mach)
        #     cl *= -1
        # else
        #     cl, cd = afeval(af, alpha, Re, Mach)
        # end

    ### Airfoil look up
    if turbine #TODO: Should I also do a negative alphadot
        alpha *= -1
    end

    alphadot = 0.0 #Todo: I need something here. 
    w = sqrt(Vx^2 + Vy^2)
    wdot = 0.0
    cl, cd = afeval(stallmodel, alpha, alphadot, w, wdot, state, chord, airfoil)

    if turbine
        cl *= -1
    end

    ### Rotate loads
    cn = cl*cphi - cd*sphi
    ct = cl*sphi + cd*cphi

    # @show cn, ct

    ### Calculate correction factor
    F = calc_F(tipcorrection, r, rhub, rtip, B, phi) 
    # @show F  #The initial thing went to NaN because of phi? -> Initialize with the twist vector. .. Didn't completely fix it.  

    ### Calculate induction factor
    k = cn*sigma/(4.0*F*sphi*sphi)

    # @show sigma, k #These are both the same as CCBlade. s
    # @show #sphi  #This is the correct value. 

    if phi < 0
        k *= -1
    end

    if k >= -2.0/3  # momentum region #TODO: Dr. Ning updated this in a newer version of CCBlade. It makes a difference in calculating a. 
        a = k/(1 - k)
        # println("Momentum Region")

    else  # empirical region
        # println("empirical region")
        g1 = F*(2*k - 1) + 10.0/9
        g2 = F*(F - 2*k - 4.0/3)
        g3 = 2*F*(1 - k) - 25.0/9

        if isapprox(g3, 0.0, atol=1e-6)  # avoid singularity
            a = 1.0/(2.0*sqrt(g2)) - 1
        else
            a = (g1 + sqrt(g2)) / g3
        end
    end

    # u = a * Vx

    ## Tangential induction  #I'm not sure that I even need to calculate the tangential angle. I don't think that I pass the induced velocities to the DS model. -> But I might want the loads.... which Dr. Ning has so graciously provided. :) -> Also... equation 6.43 of the notes

    kp = ct*sigma/(4.0*F*sphi*cphi) #Copied and pasted from CCBlade. 

    if Vx < 0
        kp *= -1
    end

    # if isapprox(kp, -1.0, atol=1e-6)  # state corresopnds to Vy=0, return any nonzero residual
    #     return 1.0, Outputs()
    # end

    ap = kp/(1 + kp)
    # v = ap * Vy

    # @show a, ap

    # # ------- loads ---------
    # W = sqrt((Vx + u)^2 + (Vy - v)^2)
    # Np = cn*0.5*rho*W^2*chord
    # Tp = ct*0.5*rho*W^2*chord

    # The BEM methodology applies hub/tip losses to the loads rather than to the velocities.  
    # This is the most common way to implement a BEM, but it means that the raw velocities are misleading 
    # as they do not contain any hub/tip loss corrections.
    # To fix this we compute the effective hub/tip losses that would produce the same thrust/torque.
    # In other words:
    # CT = 4 a (1 + a) F = 4 a G (1 + a G)\n
    # This is solved for G, then multiplied against the wake velocities.
    
    # if isapprox(Vx, 0.0, atol=1e-6)
    #     G = sqrt(F)
    # elseif isapprox(Vy, 0.0, atol=1e-6)
    #     G = F
    # else
    #     G = (-1.0 + sqrt(1.0 + 4*a*(1.0 + a)*F))/(2*a)
    # end
    # u *= G
    # v *= G


    ### Calculate Inflow angle. 
    phi = atan(Vx*(1+a), Vy*(1-ap))   

    return phi, F

    # if rotor.turbine ####### FROM CCBLADE. 
    #     return R, Outputs(-Np, -Tp, -a, -ap, -u, -v, phi, -alpha, W, -cl, cd, -cn, -ct, F, G)
    # else
    #     return R, Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)
    # end

    #Todo: I need to correct this depending if it is for a turbine or a propeller.
    # if turbine #Todo: I'm not sure this is working correctly
    #     return -phi #TODO: Pass out loadings
    # else
    #     return phi
    # end
end




function update_inflowangle_blade!(phivec, Vxvec, Vyvec, states, rvec, chordvec, twistvec, blade, pitch, rhub, rtip, B, turbine, tipcorrection, stallmodel, Fvec)
    na = length(rvec)

    for i = 1:na
        idx = 4*(i-1)
        phivec[i], Fvec[i] = update_inflowangle(phivec[i], Vxvec[i], Vyvec[i], states[idx + 1:idx+4], rvec[i], chordvec[i], twistvec[i], blade.airfoils[i], pitch, rhub, rtip, B, turbine, tipcorrection, stallmodel)
    end
end

function fixedpointbem(rvec, chordvec, twistvec, blade::Blade, pitch, rhub, rtip, B, env::Environment; turbine=true, tipcorrection=Prandtltiphub(), stallmodel=NoStall(), tolerance=1e-5, maxiterations=100, alpha = 0.7, verbose=false, speakiter=100)
    na = length(rvec)

    Vxvec = [env.Vinf(0) for i = 1:na]
    Vyvec = [env.RS(0)*rvec[i] for i = 1:na]

    phivec = deepcopy(twistvec)
    oldphivec = deepcopy(twistvec.*2)
    residuals = zeros(na)
    converged = false
    iter = 0
    states = Array{Float64}(undef, 4*na)

    Fvec = zeros(na)

    while !converged
        update_inflowangle_blade!(phivec, Vxvec, Vyvec, states, rvec, chordvec, twistvec, blade, pitch, rhub, rtip, B, turbine, tipcorrection, stallmodel, Fvec)

       
        residuals .= abs.(phivec .- oldphivec)

        for i = 1:na #relaxation
            phivec[i] = oldphivec[i] + (phivec[i] - oldphivec[i])*alpha
        end

        for i = 1:na #Copy by value, not by reference. #TODO: There's probably a more elegant way to do this. 
            oldphivec[i] = phivec[i]
        end

        iter += 1

        if maximum(residuals) < tolerance #& iter>1
            converged = true
            if verbose
                println("Converged at $iter iterations")
            end
            break
        end

        if iter >= maxiterations
            if verbose
                println("Did not converge")
            end
            break
        end

        if verbose && mod(iter, speakiter)==0
            println("Iteration $iter  - Max residual: ", maximum(residuals))
        end
    end

    return phivec, Fvec, converged
end