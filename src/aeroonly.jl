#=
An aerodynamics only coupling. So just the dynamic stall model and CCBlade in a time-marching approach. The bem and ds model are tightly coupled, but in a format to solve only one time step at a time. 


Adam Cardoza 8/1/22
=#


function getbemrisoy(phi, x_riso, chord, twist, pitch, Vx, Vy, Vxdot, Vydot, airfoil::Airfoil)

    ### Find velocities in dynamic stall frame. (Not accounting for induced velocities in Dynamic stall model.)
    u = sqrt(Vy^2 + Vx^2)  
    v = 0 

    udot = sqrt(Vxdot^2 + Vydot^2)
    vdot = 0.0

    alpha = -((twist + pitch) - phi)
    alphadot = 0.0 

    yriso = SVector(u, udot, v, vdot, alpha, alphadot)

    ### Create Airfoil object 
    Cl, Cd = riso_coefs(x_riso, yriso, chord, airfoil)

    ybem = SVector(Cl, Cd, Vx, Vy)

    return ybem, yriso
end


"""

p = radius_vec, chord_vec, twist_vec
x_section = [risostates, bemcontraints]
"""
function createbemrisofun(B, rvec, chordvec, twistvec, rhub, rtip, hubHt, pitch, blade::Blade, env::Environment; turbine=true, tipcorrection=CCBlade.PrandtlTipHub(), shearexp=0.0)
    n = length(blade.airfoils)
    
    afarray = Array{CCBlade.AFType}(undef,1)
    oparray = Array{CCBlade.OperatingPoint}(undef, 1)
    
    alphavec = collect(-pi:.1:pi)
    rotor = CCBlade.Rotor(rhub, rtip, B, precone=0.0, turbine=turbine)


    function bemrisofun(outs, dx, x, p, t)
        clvec = ones(eltype(outs), length(alphavec))  
        cdvec = ones(eltype(outs), length(alphavec))
        # @show t #This is changing. 

        for i=1:n
            bem_idx = 4*n + i
            xbem = view(x, bem_idx) 

            riso_idx = 4*(i-1)
            dxs_riso = view(dx, 1+riso_idx:riso_idx+4) 
            xs_riso = view(x, 1+riso_idx:riso_idx+4)  

            Vx = env.Vinf(t)
            Vy = env.RS(t)*rvec[i]
            Vxdot = env.Vinfdot(t)
            Vydot = env.RSdot(t)*rvec[i]

            ys_bem, ys_riso = getbemrisoy(xbem[1], xs_riso, chordvec[i], twistvec[i], pitch, Vx, Vy, Vxdot, Vydot, blade.airfoils[i])

            outs[1+riso_idx:riso_idx+4] = riso_residual(dxs_riso, xs_riso, ys_riso, chordvec[i], t, blade.airfoils[i])
            # if i==n These aren't changing in time. Which implies that x and dx aren't changing. 
            #     @show outs[1+riso_idx:riso_idx+4]
            # end

            U = sqrt(ys_bem[3]^2 + ys_bem[4]^2)
            clvec[:] .= ys_bem[1]
            cdvec[:] .= ys_bem[2]

            afarray[1] = CCBlade.AlphaAF(alphavec, clvec, cdvec, "", env.rho*U*chordvec[i]/env.mu, U/env.a) #Note: We definitely want to leave the airfoils in the replacing array, they take a crap ton of space. #TODO: If this works, then we'll try using the dynamic stall states. #Question: Why am I including the Reynold's number? It doesn't make any difference in the actual call.  

            ### Create section object
            section = CCBlade.Section(rvec[i], chordvec[i], twistvec[i], afarray[1])

            ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. #TODO: I can just directly create the operating point. 
            oparray[1] = CCBlade.windturbine_op(ys_bem[3], ys_bem[4]/rvec[i], pitch, rvec[1], 0.0, 0.0, 0.0, 0.0, hubHt, shearexp, env.rho, env.mu, env.a)

            get_bem_residual!(outs, bem_idx, xbem[1], rotor, section, oparray[1])
        end
    end
    return bemrisofun
end


"""
    simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, tspan, dt; solver=RK4(), verbose=false, coupling=OneWay())

Simulates the described wind turbine across the given time span using a BEM (I am using CCBlade.) that is extended by the Riso dynamic stall method. A structural response is simulated by GXBeam. The BEM is assumed to be constant across a given time step. The BEM is influenced by the current state of the dynamic stall method. The dynamic stall and structural models are then integrated forward based on the BEM solution. The dynamic stall model is integrated using the provided solver (defaults to a Runge-Kutta 4th order) and the structural model is solved using it's internal dynamic solver. 

### Inputs
- rvec::Array{TF, 1} - A vector holding the radial distance from the center of rotation of all of the aerodynamic nodes (meters)
- chordvec::Array{TF, 1} - A vector holding the chord lengths of the blade at all of the aerodyanmic nodes (meters)
- twistvec::Array{TF, 1} - A vector holding the aerodynamic twist at every aerodynamic node (radians)
- rhub::TF - the radial location of the hub (meters) #TODO: Should be able to remove this since it's in the blade. 
- rtip::TF - the radial location of the tip (meters) # ""
- blade::Blade - a blade struct
- env::Environment
- tspan::Tuple{TF, TF}
- dt::TF
- solver::Solver - the type of integration method used for integrating the Riso states. Defaults to RK4()
- verbose::Bool - Whether or not to print out information at every time step. 
- coupling::Coupling - a struct to declare to the model which type of coupling used. OneWay() or TwoWay(). 

### Outpus


### Notes
- The structural model does not include damping. 
- The initial condition is assumed to be the loads from the undeflected solution (to add the option to start from the fixed point iteration solution later). Those loads are then used to find the steady state deflections and inflow velocities. The inflow velocities are used to calculate the steady state states for the Riso model.
"""
function simulate(rvec, chordvec, twistvec, rhub, rtip, hubHt, blade, env, tspan, dt, B; solver=DBDF1!, verbose=false, turbine=true, tipcorrection=CCBlade.PrandtlTipHub(), shearexp=0.0)

    t0 = tspan[1]
    tvec = tspan[1]:dt:tspan[2]
    nt = length(tvec)

    ### Create CCBlade inputs
    rotor = CCBlade.Rotor(rhub, rtip, B, precone=0.0, turbine=true)
    na = length(blade.airfoils)
    sections = [Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i in 1:na] 


    ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
    Vxvec = [env.Vinf(t0) for i in 1:na] #Freestream velocity
    Vyvec = [env.RS(t0)*rvec[i] for i in 1:na]
    pitch = 0.0 #TODO: I'm not sure if this pitch does something or not. And I don't know if I want to pitch yet. 
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na] 
    

     
    ### Solve BEM residual
    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints) 

    


    ### Prepare DAE for solution. 
    dae = createbemrisofun(B, rvec, chordvec, twistvec, rhub, rtip, hubHt, pitch, blade::Blade, env::Environment; turbine=turbine, tipcorrection=tipcorrection, shearexp=shearexp)

    p_aero = nothing

    x0_ds = zeros(4*na) #TODO: I think that one of these states might need to be initialized as phi. #TODO: There's room to initialize this with a steady state solution. 
    #Note: I wonder if the first three states are in an intersting spot, so if I artifically push them positive if the steady state solve will find something more reasonable. 
    # x0[1:4:4*na] .= 0.01 
    # x0[2:4:4*na] .= 0.01
    # x0[3:4:4*na] .= 0.5
    #Note: Yeah, that did nothing. The steady state solution was still the same. I wonder what happens if I don't initialize with the steady state solution. 

    x0_ds[4:4:4*na] .= 1.0 

    x0 = vcat(x0_ds, ccout.phi)
    n = length(x0) #Number of states and constraints
    
    

    differential_vars = [i<=4*na ? true : false for i in 1:n] 

    solverobj = solver(differential_vars)

    x = zeros(nt, solverobj.n)
    dx = zeros(nt, solverobj.nd)
    outs = zeros(solverobj.n)

    x[1,:] = x0
    # dx[1,:] = dx0 #TODO: I don't know what to put this as. Like... I could assume that it starts from steady state or something.

    if isa(solverobj, DBDF1!)
        start = 2
    elseif isa(solverobj, DBDF2!)
        start = 3
        initialsolver = DBDF1!(differential_vars)
        dt0 = tvec[2]-tvec[1]
        x[2,:], dx[2,:] = initialsolver(dae, outs, x[1,:], dx[1,:], p_aero, t0, dt0)
    end


    ### Iterate through time steps
    for i = start:nt
        ## Update timestep
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]
        if isa(solverobj, DBDF1!)
            indices = (i-1,:)
        elseif isa(solverobj, DBDF2!)
            indices = (i-2:i-1,:)
        end

        x[i,:], dx[i,:] = solverobj(dae, outs, x[indices...], dx[i-1,:], p_aero, t, dt) #Todo: There might be something wrong with the solver. I'm wondering if it's setting dxn equal to something that causes it not to move. I'm not sure. Because the dynamic stall can't be at steady state. -> I could wrap the DifferentialEquations solver and see what happens. -> It was the differential_vars... I think. -> Well it solves, but it doesn't look like it is solving right. Or yeah... It's also really slow. 


        if verbose
            println("t=", round(t, digits=4))
        end
    end

    if verbose
        println("t=", round(tvec[end], digits=4))
    end

    return x, dx, tvec
end

