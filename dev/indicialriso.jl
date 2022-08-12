#=


=#


function seperation_point(alpha, afm, afp, clfit, dcldalpha, alpha0)
    #TODO: I'm not really sure that using the minimum of these two is really the way to avoid the problem of this blowing up to infinity. (When alpha=alpha0) (This check happens in the if statement.)

    #TODO: I'm not sure that using the absolute value here is the correct way to dodge the problem of crossing the x axis at different times.


    if !(afm < alpha < afp) #Check if alpha is in the bounds. 
        # println("f was set to zero. ")
        return typeof(alpha)(0)
    end
    
    cl_static = clfit(alpha)
    cl_linear = dcldalpha*(alpha-alpha0)
    f = (2*sqrt(abs(cl_static/cl_linear))-1)^2

    if f>1 #Question: What if I don't return this? I might get Inf.... or possibly NaN... but I will less likely get 1.0... which is my problem child in the seperated coefficient of lift function. -> I fixed the fully seperated coefficient of lift function... I just plugged this function inside the other and simplified. 
        return typeof(alpha)(1)
    elseif isnan(f)
        # println("f return NaN")
        return typeof(alpha)(1)
    end

    return f
end

function update_states(xj, xjm, deltat, uj, ujm, udotj, udotjm, alphaj, alphajm, alphadotj, alphadotjm, c, dcldalpha, alpha0, afm, afp, A1, A2, b1, b2, Tp, Tf, clfit)
    Pi = b1*((uj + ujm)/c) + (udotj + udotjm)/(uj + ujm)
    Qi = b1*A1*(ujm*alphajm + uj*alphaj)/c 
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[1] = Ci*xjm[1] + Ii

    Pi = b2*((uj + ujm)/c) + (udotj + udotjm)/(uj + ujm)
    Qi = b2*A2*(ujm*alphajm + uj*alphaj)/c
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[2] = Ci*xjm[2] + Ii

    alphaEjm = alphajm*(1-A1-A2) + xjm[1] + xjm[2]
    alphaEj = alphaj*(1-A1-A2) +xj[1] + xj[2]
    t1 = dcldalpha*(alphaEjm + alphaEj - 2*alpha0)
    t2 = pi*c*(alphadotjm/ujm + alphadotj/uj)/2

    Pi = 1/Tp
    Qi = (t1 + t2)/(2*Tp)
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[3] = Ci*xjm[3] + Ii

    alphafj = xj[3]/dcldalpha + alpha0 
    alphafjm = xjm[3]/dcldalpha + alpha0 
    fj = seperation_point(alphafj, afm, afp, clfit, dcldalpha, alpha0)
    fjm = seperation_point(alphafjm, afm, afp, clfit, dcldalpha, alpha0)

    Pi = 1/Tf
    Qi = (fj + fjm)/(2*Tf)
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[4] = Ci*xjm[4] + Ii
    return xj
end

function update_states_blade()
end

function riso_coefs(X, y, c, airfoil)

    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    liftfit = airfoil.cl
    dragfit = airfoil.cd
    A1 = airfoil.A[1]
    A2 = airfoil.A[2]
    b1 = airfoil.b[1]
    b2 = airfoil.b[2]
    Tp = airfoil.T[1]
    Tf = airfoil.T[2]

    u, udot, v, vdot, theta, thetadot = y
    # @show theta #Okay, it's passing in the correct value here. 

    U = sqrt(u^2 + v^2)
    Udot = (2*u*udot + 2*v*vdot)/(2*sqrt(u^2 + v^2))

    psi = atan(v, u)
    alpha = psi + theta #Turbine formulation
    # @show alpha #This is the correct value

    psidot = (u*vdot - v*udot)/((u^2)*(((v/u)^2)+1))
    alphadot = psidot + thetadot

    ae = alpha*(1-A1-A2) + X[1] + X[2] #With the current alpha values, the angle of attack makes no difference, regardless of the state. 
    Tu = c/(2*U)

    # @show ae
    clfs = Clfs(ae, liftfit, dcldalpha, alpha0)
    # @show X

    Cl = dcldalpha*(ae-alpha0)*X[4] + clfs*(1-X[4]) + pi*Tu*alphadot

    fae = fst(ae, liftfit, dcldalpha, alpha0) 
    
    fpp = X[4]
    if X[4] < 0
        fpp = 0
    end
    fterm = (sqrt(fae)-sqrt(fpp))/2 - (fae-fpp)/4 #Todo: Is there a way to re-write this line so that it takes a negative argument and gives us what we want? 
    Cd = dragfit(ae) + (alpha-ae)*Cl + (dragfit(ae)-dragfit(alpha0))*fterm
    return SVector(Cl, Cd)
end