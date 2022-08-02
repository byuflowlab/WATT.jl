"""
Riso{TF} <: AbstractModel

Aerodynamic model based on the Beddoes-Leishman (1989,1990) papers with state variables ``x_1,
x_2, x_3, x_4, x_5, x_6, x_7, x_8, C'_N, f'', \\tau_v, C^v_N, f_m''``, inputs ``u, \\dot{u}, v, \\dot{v} \\omega, ,``, and parameters `` m, T_p, T_f, T_vl, T_v, eta``
"""
struct Riso
end




# --- Internal Methods --- #

function fst(alpha, liftfit, dcdalpha, alpha0)
    
    f =  (2*sqrt(abs(liftfit(alpha)/(dcdalpha*(alpha-alpha0)))) - 1)^2
    if f>= 1 || isnan(f)
        return 1.0
    elseif f<0
        return 0.0
    else
        return f
    end
    
end

function find_seperation_alpha(liftfit, dcldalpha, alpha0; n=10)

    ### Create a residual function to solve. 
    residual(alpha) = abs(dcldalpha*(alpha-alpha0)/4) - abs(liftfit(alpha))

    ### Find a bracketing range for the positive fully stalled angle. 
    aoa = range(0,pi/2; length = n)
    rng = zeros(2)
    for i = 1:n-1
        if residual(aoa[i])*residual(aoa[i+1])<0
            rng[:] = aoa[i:i+1]
            break
        end
    end
    
    ### Find the positive fully stalled angle
    alpha_positive = find_zero(residual, rng, Bisection())

    ### Bracket the negative fully stalled. 
    aoa = range(0,-pi/2; length = n)
    for i = 1:n-1
        if residual(aoa[i])*residual(aoa[i+1])<0
            rng[:] = aoa[i:i+1]
            break
        end
    end

    ### Find the negative fully stalled angle
    alpha_negative = find_zero(residual, rng, Bisection())
    return alpha_positive, alpha_negative
end

function seperationpoint(alpha, afm, afp, clfit, dcldalpha, alpha0)

    if !(afm < alpha < afp) #Check if alpha is in the bounds of the fully seperated angles of attack. 
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

    #Todo. Hansen must have some sort of switch that stops this function from reattaching when the aoa gets really high. -> like the one where you automatically set f=0 when you're outside the bounds of afm, afp
    return f
end

function Clfs(alpha, liftfit, dcldalpha, alpha0)
    f = fst(alpha, liftfit, dcldalpha, alpha0)
    Cl = 0.0
    if f>= 1.0 
        Cl =  liftfit(alpha)/2
    else
        Cl = (liftfit(alpha) - (dcldalpha*(alpha-alpha0))*f)/(1-f)  
    end

    if isnan(Cl)
        Cl = liftfit(alpha)/2
    end

    return Cl
end


function riso_states(X, u, udot, v, vdot, theta, thetadot, c, airfoil)
    #X[1] - 
    #X[2] - 
    #X[3] - 
    #X[4] - 
    # u - 
    # theta - Angle of attack

    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    liftfit = airfoil.cl
    A1 = airfoil.A[1]
    A2 = airfoil.A[2]
    b1 = airfoil.b[1]
    b2 = airfoil.b[2]
    Tp = airfoil.T[1]
    Tf = airfoil.T[2]

    ### Convert from u, v, omega to U, alpha and q 
    U = sqrt(u^2 + v^2)
    Udot = (2*u*udot + 2*v*vdot)/(2*sqrt(u^2 + v^2))

    psi = atan(v, u) #Note: psi \dne phi
    alpha = psi + theta #Turbine formulation

    # println(psi)

    psidot = (u*vdot - v*udot)/((u^2)*(((v/u)^2)+1))
    alphadot = psidot + thetadot
    # println(psidot)
    
    ### Calculate constants
    Tu = c/(2*U) #This will go to NaN if U=0
    # println(typeof(X[1]))
    # @show U, Tu, c
    # @show X[1]
    # @show Tp
    
    ### Calculate the state rates
    dx1 = (b1*A1*alpha/Tu) - X[1]*(b1+ (c*Udot/(2*(U)^2)))/Tu 

    dx2 = (b2*A2*alpha/Tu) - X[2]*(b2+ (c*Udot/(2*(U)^2)))/Tu

    ae = alpha*(1-A1-A2) + X[1] + X[2] #Effective Angle of Attack
    # @show dcldalpha, ae, alpha0, Tu, alphadot, Tp, X[3]
    dx3 = (dcldalpha*(ae-alpha0) + pi*Tu*alphadot)/Tp - X[3]/Tp #TODO: Should this be able to go negative? 

    alphaf = (X[3]/dcldalpha)+alpha0 #Seperation Angle of Attack
    fp = fst(alphaf, liftfit, dcldalpha, alpha0)
    dx4 = fp/Tf - X[4]/Tf 

    return SVector(dx1, dx2, dx3, dx4)
end

function riso_residual(dx, x, y, p, t, airfoil)  

    ### Extract inputs
    u, udot, v, vdot, theta, thetadot = y

    # println(y)

    ### Extract parameters
    c = p 

    ### Calculate Aero state rates and factors for BEM 
    d_X = riso_states(x, u, udot, v, vdot, theta, thetadot, c, airfoil)

    ### Dynamic Stall Model State Rate Residuals
    r1 = dx[1] - d_X[1]
    r2 = dx[2] - d_X[2]
    r3 = dx[3] - d_X[3]
    r4 = dx[4] - d_X[4] 
    
    if isnan(r1)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    elseif isnan(r2)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    elseif isnan(r3)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    elseif isnan(r4)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    end
    # println("riso residual: ", [r1, r2, r3, r4])
    return SVector(r1, r2, r3, r4)
end

function get_riso_y(twist, env, frequency, amplitude, t) 

    u = env.Vinf(t)
    udot = env.Vinfdot(t)

    v = 0.0
    vdot = 0.0

    theta = twist + amplitude*cos(frequency*t)
    thetadot = -amplitude*frequency*sin(frequency*t)

    return [u, udot, v, vdot, theta, thetadot]
end

function create_risofun(twistvec, blade::Blade, env::Environment, frequency, amplitude)
    function risofun(outs, dx, x, p, t)
        n = length(blade.airfoils)
        for i = 1:n
            idx = 4*(i-1)
            dxs = dx[1+idx:idx+4]
            xs = x[1+idx:idx+4]
            ps = p[i]
            ys = get_riso_y(twistvec[i], env, frequency, amplitude, t) #Steady

            outs[1+idx:idx+4] = riso_residual(dxs, xs, ys, ps, t, blade.airfoils[i])
        end
    end
    return risofun
end

function differentialvars(model::Riso, n)
    return fill(true, 4*n)
end

function create_risoODE(twistvec, blade::Blade, env::Environment, frequency, amplitude)
    function risofun(x, p, t)
        n = length(blade.airfoils)
        dx = zeros(4*n)
        for i = 1:n
            idx = 4*(i-1)
            xs = x[1+idx:idx+4]
            ps = p[i]
            ys = get_riso_y(twistvec[i], env, frequency, amplitude, t) 
            u, udot, v, vdot, theta, thetadot = ys
            c = ps

            dx[1+idx:4+idx] = riso_states(xs, u, udot, v, vdot, theta, thetadot, c, blade.airfoils[i])
        end
        return dx
    end
    return risofun
end

function parsesolution(model::Riso, blade::Blade, env::Environment, p, sol, twistvec, frequency, amplitude)
    u = Array(sol)'
    t = sol.t
    m = length(t)
    n = length(blade.airfoils)
    Cl = zeros(m, n)
    Cd = zeros(m, n)

    for j = 1:m
        ut = u[j,:]
        for i = 1:n
            idx = 4*(i-1)
            xs = ut[1+idx:idx+4]
            ps = p[i]
            ys = get_riso_y(twistvec[i], env, frequency, amplitude, t[i]) 
    
            Cl[j, i], Cd[j, i] = riso_coefs(xs, ys, ps, blade.airfoils[i])
        end
    end
    return t, Cl, Cd
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

function riso_coefs_funs(X, U, c, alphadot, airfoil)
    dcldalpha = airfoil.dcldalpha 
    alpha0 = airfoil.alpha0
    liftfit = airfoil.cl
    dragfit = airfoil.cd
    A1 = airfoil.A[1]
    A2 = airfoil.A[2]
    

    Tu = c/(2*U)

    function clfun(alpha)
        ae = alpha*(1-A1-A2) + X[1] + X[2]

        clfs = Clfs(ae, liftfit, dcldalpha, alpha0)

        return dcldalpha*(ae-alpha0)*X[4] + clfs*(1-X[4]) + pi*Tu*alphadot
    end

    function cdfun(alpha)
        ae = alpha*(1-A1-A2) + X[1] + X[2]

        clfs = Clfs(ae, liftfit, dcldalpha, alpha0)

        Cl = dcldalpha*(ae-alpha0)*X[4] + clfs*(1-X[4]) + pi*Tu*alphadot

        fae = fst(ae, liftfit, dcldalpha, alpha0) 
        
        fterm = (sqrt(fae)-sqrt(X[4]))/2 - (fae-X[4])/4 
        return dragfit(ae) + (alpha-ae)*Cl + (dragfit(ae)-dragfit(alpha0))*fterm
    end
    return clfun, cdfun
end


