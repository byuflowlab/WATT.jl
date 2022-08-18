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
function simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec; turbine::Bool=true, dsmodel::DS.DSModel=DS.riso(blade.airfoils), dsmodelinit::ModelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, azimuth0=0.0)

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

    rotor = Rotor(rhub, rtip, B; precone, turbine)
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
    ode, xds, p_ds = initializeDSmodel(dsmodel, dsmodelinit, solver, turbine, nt, na, tvec, ccout.W, Wdotvec, chordvec, twistvec, ccout.phi, pitch) 




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



