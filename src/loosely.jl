#=

Recreating what they do in AeroDyn, which is solve the BEM, feed the inflow angle to the dynamic stall model, then calculate the loading. 

Adam Cardoza 8/6/22
=# 

export simulate

#Question: I don't know if using this many structs will make my package take forever to compile... So it might be better to create a nest of if statements... but we'll see.



"""
    simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec; turbine::Bool=true, dsmodel::DS.DSModel=DS.riso(blade.airfoils), dsmodelinit::ModelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, azimuth0=0.0)

Simulate the rotor's response for a given rotor and environmental condition. 

### Inputs
- rvec::Array{TF, 1}
- chordvec::Array{TF, 1}

### Outputs


### Notes
"""
function simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec; turbine::Bool=true, tipcorrection=CCBlade.PrandtlTipHub(), dsmodel::DS.DSModel=DS.riso(blade.airfoils), dsmodelinit::ModelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, azimuth0=0.0)
    if isa(dsmodel.detype, DS.Functional)
        error("Rotors isn't set up to simulate with functional forms of the dsmodel yet, choose the iterative model.")
    elseif isa(dsmodel.detype, DS.Iterative)
        return simulate_iterative(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade, env, tvec, turbine, tipcorrection, dsmodel, dsmodelinit, solver, verbose, speakiter, azimuth0)
    elseif isa(dsmodel.detype, DS.Indicial)
        if isa(dsmodel, DS.BeddoesLeishman)&&isa(dsmodelinit, Hansen)
            if verbose
                @warn("The initialization for the BeddoesLeishman dynamic stall model was Hansen()... switched to BeddoesLeishman().")
            end
            dsmodelinit = BeddoesLeishman()
        end
        return simulate_indicial(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade, env, tvec, turbine, tipcorrection, dsmodel, dsmodelinit, solver, verbose, speakiter, azimuth0)
    end
      
end



function simulate_iterative(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec, turbine::Bool, tipcorrection, dsmodel::DS.DSModel, dsmodelinit::ModelInit, solver::Solver, verbose::Bool, speakiter, azimuth0)

    if verbose
        println("Rotors.jl initializing solution...")
    end

    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. 
    #TODO: I might put some checks in to see if the airfoils are in degrees or radians. 
    #TODO: I might put some checks in to see if twist or pitch is in degrees or radians. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). 
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types. 

    ### Initialization information
    na = length(rvec)
    nt = length(tvec)

    t0 = tvec[1]

    Vxvec = Array{eltype(chordvec)}(undef, na)
    Vyvec = Array{eltype(chordvec)}(undef, na)

    azimuth = Array{eltype(chordvec)}(undef, nt)




    ### Initialize BEM solution
    azimuth[1] = azimuth0

    rotor = Rotor(rhub, rtip, B; precone, turbine, tip=tipcorrection)
    sections = [CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i = 1:na]

    # Vxvec = [env.Vinf(t0) for i in 1:na]
    # Vyvec = [env.RS(t0)*rvec[i]*cos(precone) for i in 1:na]

    for i = 1:na
        Vxvec[i], Vyvec[i] = get_aero_velocities(env, t0, rvec[i], azimuth[1], precone, tilt, yaw, hubht)
    end
    
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na]

    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)




    ### Initialize DS solution
     

    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2) for i in 1:na] #Todo: I probably need to update if there is precone, tilt, yaw, etc. -> Maybe I'll make a function to do this. 

    # p_ds = vcat(chordvec, twistvec, ccout.phi, ccout.W, Wdotvec, pitch) #p = [chordvec, twistvec, phivec, Wvec, Wdotvec, pitch]

    # xds = initializeDSstates(dsmodelinit, nt, na, p_ds, risoode)
    ode, xds, p_ds = initializeDSmodel(dsmodel, dsmodelinit, solver, turbine, nt, na, tvec, ccout.W, Wdotvec, chordvec, twistvec, ccout.phi, pitch, env.a) 


    # @show size(xds)

    ### Prepare data storage
    cchistory = Array{Array{Outputs{eltype(rvec)}, 1}, 1}(undef, nt) #[ccout for i = 1:nt]
    Cl = Array{eltype(rvec)}(undef,(nt, na))
    Cd = Array{eltype(rvec)}(undef,(nt, na))
    Cn = Array{eltype(rvec)}(undef,(nt, na))
    Ct = Array{eltype(rvec)}(undef,(nt, na))
    N = Array{eltype(rvec)}(undef,(nt, na))
    T = Array{eltype(rvec)}(undef,(nt, na))

    cchistory[1] = ccout
    N[1,:], T[1,:], Cn[1,:], Ct[1,:], Cl[1,:], Cd[1,:] = extractloads(dsmodel, xds[1,:], cchistory[1], t0, rvec, chordvec, twistvec, pitch, blade, env)
    
    


    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. 

        ### Update BEM inputs
        for j = 1:na
            # Vxvec[j] = env.Vinf(t)
            # Vyvec[j] = env.RS(t)*rvec[j]*cos(precone)
            Vxvec[j], Vyvec[j] = get_aero_velocities(env, t, rvec[j], azimuth[i], precone, tilt, yaw, hubht)
            operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
        end

        ### Solve BEM
        cchistory[i] = CCBlade.solve.(Ref(rotor), sections, operatingpoints) #TODO: Write a solver that is initialized with the previous inflow angle. 


        ### Update Dynamic Stall model inputs
        update_aero_parameters!(dsmodel, turbine, p_ds, na, rvec, cchistory[i].W, cchistory[i].phi, twistvec, pitch, env, t)



        ### Integrate Dynamic Stall model
        xds[i,:] = solver(ode, xds[i-1,:], p_ds, t, dt)



        ### Extract loads
        N[i,:], T[i,:], Cn[i,:], Ct[i,:], Cl[i,:], Cd[i,:] = extractloads(dsmodel, xds[i,:], cchistory[i], t, rvec, chordvec, twistvec, pitch, blade, env) 
        if verbose & (mod(i, speakiter)==0)
            println("Simulation time: ", t)
        end
    end

    Fx = N.*cos(precone) #TODO: The N and T might need some rotating for precone. I'm not sure that I'm doing this right.
    Fy = T.*cos(precone) #I copied Dr. Ning's code to get the thrust and torque, but it seems odd that they are both multiplied by the cosine of precone. Which makes sense. 

    return (N=N, T=T, Fx=Fx, Fy=Fy), cchistory, xds  
end

function simulate_indicial(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec, turbine::Bool, tipcorrection, dsmodel::DS.DSModel, dsmodelinit::ModelInit, solver::Solver, verbose::Bool, speakiter, azimuth0)

    if verbose
        println("Rotors.jl initializing solution...")
    end
    
    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. 
    #TODO: I might put some checks in to see if the airfoils are in degrees or radians. 
    #TODO: I might put some checks in to see if twist or pitch is in degrees or radians. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). 
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types. 

    ### Initialization information
    na = length(rvec)
    nt = length(tvec)

    t0 = tvec[1]

    Vxvec = Array{eltype(chordvec)}(undef, na)
    Vyvec = Array{eltype(chordvec)}(undef, na)

    azimuth = Array{eltype(chordvec)}(undef, nt)




    ### Initialize BEM solution
    azimuth[1] = azimuth0

    rotor = Rotor(rhub, rtip, B; precone, turbine, tip=tipcorrection)
    sections = [CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i = 1:na]

    # Vxvec = [env.Vinf(t0) for i in 1:na]
    # Vyvec = [env.RS(t0)*rvec[i]*cos(precone) for i in 1:na]

    for i = 1:na
        Vxvec[i], Vyvec[i] = get_aero_velocities(env, t0, rvec[i], azimuth[1], precone, tilt, yaw, hubht)
    end
    
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na]

    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)




    ### Initialize DS solution
     

    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2) for i in 1:na] #Todo: I probably need to update if there is precone, tilt, yaw, etc. -> Maybe I'll make a function to do this. 

    # p_ds = vcat(chordvec, twistvec, ccout.phi, ccout.W, Wdotvec, pitch) #p = [chordvec, twistvec, phivec, Wvec, Wdotvec, pitch]

    # xds = initializeDSstates(dsmodelinit, nt, na, p_ds, risoode)
    ode, xds, p_ds = initializeDSmodel(dsmodel, dsmodelinit, solver, turbine, nt, na, tvec, ccout.W, Wdotvec, chordvec, twistvec, ccout.phi, pitch, env.a) 
    # @show length(p_ds)
    # @show p_ds
    # @show size(xds)

    ### Prepare data storage
    cchistory = Array{Array{Outputs{eltype(rvec)}, 1}, 1}(undef, nt) #[ccout for i = 1:nt]
    Cl = Array{eltype(rvec)}(undef,(nt, na))
    Cd = Array{eltype(rvec)}(undef,(nt, na))
    Cn = Array{eltype(rvec)}(undef,(nt, na))
    Ct = Array{eltype(rvec)}(undef,(nt, na))
    Cx = Array{eltype(rvec)}(undef,(nt, na))
    Cy = Array{eltype(rvec)}(undef,(nt, na))

    N = Array{eltype(rvec)}(undef,(nt, na))
    T = Array{eltype(rvec)}(undef,(nt, na))
    Fx = Array{eltype(rvec)}(undef,(nt, na))
    Fy = Array{eltype(rvec)}(undef,(nt, na))
    # Pcopy = Array{eltype(rvec)}(undef, (nt, length(p_ds)))
    # Pcopy[1,:] = deepcopy(p_ds)

    cchistory[1] = ccout
    Cx[1,:], Cy[1,:], Cn[1,:], Ct[1,:], Cl[1,:], Cd[1,:] = extractloads(dsmodel, xds[1,:], cchistory[1], chordvec, twistvec, pitch, blade, env)
    

    ### Dimensionalize
    


    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. 

        ### Update BEM inputs
        for j = 1:na
            # Vxvec[j] = env.Vinf(t)
            # Vyvec[j] = env.RS(t)*rvec[j]*cos(precone)
            Vxvec[j], Vyvec[j] = get_aero_velocities(env, t, rvec[j], azimuth[i], precone, tilt, yaw, hubht)
            operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
        end

        # @show Vyvec
        # @show @. sqrt(Vxvec^2 + Vyvec^2)

        ### Solve BEM
        cchistory[i] = CCBlade.solve.(Ref(rotor), sections, operatingpoints) #TODO: Write a solver that is initialized with the previous inflow angle. 
        # cchistory[i] = fixedpointbem.(Ref(rotor), sections, operatingpoints, cchistory[i-1].phi)

        
        ### Update Dynamic Stall model inputs
        # println("")
        # @show p_ds[end-1:end]
        update_aero_parameters!(dsmodel, turbine, p_ds, na, rvec, cchistory[i].W, cchistory[i].phi, twistvec, pitch, env, t)
        # @show p_ds[end-1:end] #Looks like it's updating every iteration. 
        # Pcopy[i, :] = deepcopy(p_ds)


        ### Integrate Dynamic Stall model

        # @show i, t, dt

        xds[i,:] = ode(xds[i-1,:], p_ds, t, dt) #TODO: Instead of doing an entirely different function, I could have an if statement around this. 
    

        # for k=1:na
        #     idxs = 21*(k-1)+1:21*k
        #     @show xds[i,idxs]
        # end
        # println("")




        ### Extract loads #TODO: I don't think that I need to keep all of this data. I can probably define a function that calculates all the different loads. 
        Cx[i,:], Cy[i,:], Cn[i,:], Ct[i,:], Cl[i,:], Cd[i,:] = extractloads(dsmodel, xds[i,:], cchistory[i], chordvec, twistvec, pitch, blade, env) 



        ### Dimensionalize
        for j = 1:na
            u_i = cchistory[i].W[j]
            N[i,j] = Cn[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
            T[i,j] = Ct[i,j]*0.5*env.rho*u_i^2*chordvec[j]
            Fx[i,j] = Cx[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
            Fy[i,j] = Cy[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
        end



        if verbose & (mod(i, speakiter)==0)
            println("Simulation time: ", t)
        end
    end


    # Fx = N.*cos(precone) #TODO. The N and T might need some rotating for precone. I'm not sure that I'm doing this right.
    # Fy = T #.*cos(precone) #I copied Dr. Ning's code to get the thrust and torque, but it seems odd that they are both multiplied by the cosine of precone. Which makes sense. -> What on earth am I saying here? It seems odd then I turn around and say it makes sense..... thinking more about it... I'm not sure that it makes any sense whatsoever. After thinking through it, a precone should create deflection in the XZ plane, which means that the tangential force in the Y direction shouldn't get deflected (displaced, but not deflected). 

    return (N=N, T=T, Fx=Fx, Fy=Fy), (Cx=Cx, Cy=Cy, Cn=Cn, Ct=Ct, Cl=Cl, Cd=Cd), cchistory, xds, azimuth #, Pcopy
end


export rotorloads

function rotorloads(loads, rhub, rtip, rvec, B)

    nt, _ = size(loads.N)

    thrust = Array{eltype(loads.N)}(undef, nt)
    torque = Array{eltype(loads.N)}(undef, nt)

    rfull = [rhub; rvec; rtip]
    
    for i = 1:nt
        Fxfull = [0.0; loads.Fx[i,:]; 0.0]
        Fyfull = [0.0; loads.Fy[i,:]; 0.0]

        thrust[i] = B*FLOWMath.trapz(rfull, Fxfull)
        torque[i] = B*FLOWMath.trapz(rfull, Fyfull.*rfull)
    end

    return thrust, torque
end


function rotorloads(rhub, rtip, rvec, loads...)

    nt, _ = size(loads[1].N)

    thrust = Array{eltype(loads[1].N)}(undef, nt)
    torque = Array{eltype(loads[1].N)}(undef, nt)

    rfull = [rhub; rvec; rtip]
    
    for i = 1:nt
        for j = 1:length(loads)
            Fxfull = [0.0; loads[j].Fx[i,:]; 0.0]
            Fyfull = [0.0; loads[j].Fy[i,:]; 0.0]

            thrust[i] += FLOWMath.trapz(rfull, Fxfull)
            torque[i] += FLOWMath.trapz(rfull, Fyfull.*rfull)
        end
    end

    return thrust, torque
end