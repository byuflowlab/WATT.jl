#=
An aerodynamics only coupling. So just the dynamic stall model and CCBlade in a time-marching approach.  


Adam Cardoza 8/1/22
=#




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
function simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, tspan, dt, B; solver=RK4(), verbose=false, coupling=OneWay())

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
    
    # @show Vxvec #Really small if I use the ramp up environment. 
    # @show Vyvec #Seems reasonable. 
     
    ### Solve BEM residual
    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints) 

    

    ### Get steady state Riso states
    Wvec = [sqrt(Vxvec[i]^2 + Vyvec[i]^2) for i = 1:na] #Total inflow velocity (assumed) -> Question: will the inflow velocity from CCBlade be different? 
    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i])^2) for i = 1:na] #Derivative of the total inflow velocity

    risoode = createrisoode(blade) 

    p_aero = vcat(chordvec, twistvec, ccout.phi, Wvec, Wdotvec, pitch) #Create p_aero 

    x0 = zeros(4*na) #TODO: I think that one of these states might need to be initialized as phi.
    #Note: I wonder if the first three states are in an intersting spot, so if I artifically push them positive if the steady state solve will find something more reasonable. 
    # x0[1:4:4*na] .= 0.01 
    # x0[2:4:4*na] .= 0.01
    # x0[3:4:4*na] .= 0.5
    #Note: Yeah, that did nothing. The steady state solution was still the same. I wonder what happens if I don't initialize with the steady state solution. 

    x0[4:4:4*na] .= 1.0 

    ## Steady state solution of dynamic stall odes at initial condition. 
    prob = SteadyStateProblem(risoode, x0, p = p_aero)
    sol = DifferentialEquations.solve(prob)
    xds = sol.u 
    # xds = x0


    dsstates = [RisoState(xds[4*(i-1)+1:4*(i-1)+4], Wvec[i], chordvec[i], blade.airfoils[i]) for i = 1:na]


    
    if isa(coupling, ThreeWay)
        sections = [Section(rvec[i], chordvec[i], twistvec[i], dsstates[i]) for i in 1:na]
    end 

    
    cchistory = [ccout for i = 1:nt]
    dshistory = [dsstates for i = 1:nt] #TODO: The way this is allocated allows the first entry to get changed as dsstates changes. 
    dshistory[1] = deepcopy(dsstates)


    ### Iterate through time steps
    for i = 1:nt-1
        ## Update timestep
        t = tvec[i]


        ## Update CCBlade inputs
        for j = 1:na
            Vxvec[j] = env.Vinf(t)
            Vyvec[j] = env.RS(t)*rvec[j]
            Wvec[j] = sqrt(Vxvec[j]^2 + Vyvec[j]^2)
            Wdotvec[j] = sqrt(env.Vinfdot(t)^2 + (env.RSdot(t)*rvec[j])^2) #TODO: I need to get the structural accelerations in here as well.

            operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
            sections[j] = Section(rvec[j], chordvec[j], twistvec[j], dsstates[j])
        end

        # @show Vyvec
        #WorkLocation: 
        #Todo: Where do I update the sections? Right, I have to update the sections because I have the dynamic stall states... I wonder how it's having any change whatsoever. 

        # Solve BEM residual
        cchistory[i+1] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)


        ## Solve Riso states based on inflow and deflection #Todo I need to determine if I should be using the ith or the ith +1 state. -> I think I should be using the i+1 ccblade outputs, because I calculate the ccblade outputs, then I use those outputs to calculate the inputs for the dynamic stall model. Then I take a time step, calculate the new CCBlade outputs based on the current ds states. So I should be storing things into the i+1 slot, and if I rely on something that I've already solved, then I'll use the i+1 spot. 
        update_p_aero!(p_aero, twistvec, cchistory[i+1].phi, Wvec, Wdotvec) 
        

        # Note: I could just pass the angle of attack and inflow velocity from CCBlade to the DS model. 
        xds[:] = solver(risoode, xds, p_aero, t, dt) #TODO: I need to convert this to work in place. 

        update_dsstates!(xds, dsstates, Wvec, chordvec, blade)  



        dshistory[i+1] = deepcopy(dsstates)


        if verbose
            println("t=", round(t, digits=4))
        end
    end

    if verbose
        println("t=", round(tvec[end], digits=4))
    end

    return cchistory, dshistory, tvec
end

