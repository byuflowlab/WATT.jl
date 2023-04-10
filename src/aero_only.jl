#=

Recreating what they do in AeroDyn, which is solves the BEM, feeds the inflow angle to the dynamic stall model, then calculates the loading. 

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
function simulate(rvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec; turbine::Bool=true, tipcorrection=CCBlade.PrandtlTipHub(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, azimuth0=0.0)


    if verbose
        println("Rotors.jl initializing solution...")
    end

    #Todo: This will break without dsmodel exposed. 
    # if isa(dsmodel.detype, DS.Functional)
    #     error("Rotors isn't set up to simulate with functional forms of the dsmodel yet, choose the iterative model.")
    # elseif isa(dsmodel.detype, DS.Indicial)
    #     if isa(dsmodel, DS.BeddoesLeishman)&&isa(dsmodelinit, Hansen)
    #         if verbose
    #             @warn("The initialization for the BeddoesLeishman dynamic stall model was Hansen()... switched to BeddoesLeishman().")
    #         end
    #         dsmodelinit = BeddoesLeishman()
    #     end
    # end
    
    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. 
    #TODO: I might put some checks in to see if the airfoils are in degrees or radians. 
    #TODO: I might put some checks in to see if twist or pitch is in degrees or radians. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). 
    #TODO. It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types. -> I'm going to have different functions for optimization. I'm going to break it up into functions that calculate by each step. 

    ### Initialization information
    na = length(rvec)
    nt = length(tvec)

    t0 = tvec[1]

    airfoils = blade.airfoils
    chordvec = airfoils.c


    ### Initialize BEM solution
    Vxvec = Array{eltype(chordvec)}(undef, na)
    Vyvec = Array{eltype(chordvec)}(undef, na)
    azimuth = Array{eltype(chordvec)}(undef, nt)
    azimuth[1] = azimuth0

    rotor = Rotor(rhub, rtip, B; precone, turbine, tip=tipcorrection)
    sections = [CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i = 1:na]



    for i = 1:na
        Vxvec[i], Vyvec[i] = get_aero_velocities(env, t0, rvec[i], azimuth[1], precone, tilt, yaw, hubht)
    end
    
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na]

    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)






    ### Initialize DS solution
    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2) for i in 1:na] #TODO: I probably need to update if there is precone, tilt, yaw, etc. -> Maybe I'll make a function to do this. -> Currently not used by BLADG. 
    alphadotvec = zero(Wdotvec) #TODO: 

    
    xds, xds_idxs, p_ds = initialize_DS_model(airfoils, turbine, nt, tvec, ccout, Wdotvec, alphadotvec, twistvec, pitch)





    ### Prepare data storage
    cchistory = Array{Array{Outputs{eltype(rvec)}, 1}, 1}(undef, nt) #[ccout for i = 1:nt]
    # Cl = Array{eltype(rvec)}(undef,(nt, na))
    # Cd = Array{eltype(rvec)}(undef,(nt, na))
    # Cn = Array{eltype(rvec)}(undef,(nt, na))
    # Ct = Array{eltype(rvec)}(undef,(nt, na))
    Cx = Array{eltype(rvec)}(undef,(nt, na))
    Cy = Array{eltype(rvec)}(undef,(nt, na))
    Cm = Array{eltype(rvec)}(undef,(nt, na))

    # N = Array{eltype(rvec)}(undef,(nt, na))
    # T = Array{eltype(rvec)}(undef,(nt, na))
    Fx = Array{eltype(rvec)}(undef,(nt, na))
    Fy = Array{eltype(rvec)}(undef,(nt, na))


    cchistory[1] = ccout
    # Cx[1,:], Cy[1,:], Cn[1,:], Ct[1,:], Cl[1,:], Cd[1,:], Cm[1,:] = extractloads(dsmodel, xds[1,:], cchistory[1], chordvec, twistvec, pitch, blade, env)
    extract_ds_loads!(airfoils, view(xds, 1, :), xds_idxs, cchistory[1], p_ds, view(Cx, 1, :), view(Cy, 1, :), view(Cm, 1, :))
    

    ### Dimensionalize
    for j = 1:na
        u_1 = cchistory[1].W[j]
        # N[1,j] = Cn[1,j]*0.5*env.rho*u_1^2*chordvec[j] 
        # T[1,j] = Ct[1,j]*0.5*env.rho*u_1^2*chordvec[j]
        Fx[1,j] = Cx[1,j]*0.5*env.rho*u_1^2*chordvec[j] 
        Fy[1,j] = Cy[1,j]*0.5*env.rho*u_1^2*chordvec[j] 
    end







    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. 

        ### Update BEM inputs
        for j = 1:na
            Vxvec[j], Vyvec[j] = get_aero_velocities(env, t, rvec[j], azimuth[i], precone, tilt, yaw, hubht)
            operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
        end



        ### Solve BEM
        cchistory[i] = CCBlade.solve.(Ref(rotor), sections, operatingpoints) #TODO: Write a solver that is initialized with the previous inflow angle. 
        # cchistory[i] = fixedpointbem.(Ref(rotor), sections, operatingpoints, cchistory[i-1].phi)

        
        ### Update Dynamic Stall model inputs
        # update_aero_parameters!(dsmodel, turbine, p_ds, na, rvec, cchistory[i].W, cchistory[i].phi, twistvec, pitch, env, t)
        update_ds_inputs!(airfoils, p_ds, cchistory[i].W, cchistory[i].phi, twistvec, pitch, dt, turbine) #Todo: There is probably a more efficient way of passing in W, and Phi. Maybe I keep a running vector of W and phi that I just update. That way I'm not allocating every iteration. Better yet. I wonder if I should just copy out code from CCBlade and get rid of the Outputs struct and just store it in a matrix. 


        ### Integrate Dynamic Stall model
        # if isa(dsmodel.detype, DS.Indicial)
        #     xds[i,:] = ode(xds[i-1,:], p_ds, t, dt) 
        # else 
        #     xds[i,:] = solver(ode, xds[i-1,:], p_ds, t, dt)
        # end
        update_ds_states!(solver, airfoils, view(xds, i-1, :), view(xds, i, :), xds_idxs, p_ds, t, dt)





        ### Extract loads #TODO: I don't think that I need to keep all of this data. I can probably define a function that calculates all the different loads. -> If my vectors are predefined, it shouldn't matter how many I keep. 
        # Cx[i,:], Cy[i,:], Cn[i,:], Ct[i,:], Cl[i,:], Cd[i,:], Cm[i,:] = extractloads(dsmodel, xds[i,:], cchistory[i], chordvec, twistvec, pitch, blade, env) 
        extract_ds_loads!(airfoils, view(xds, i, :), xds_idxs, cchistory[i], p_ds, view(Cx, i, :), view(Cy, i, :), view(Cm, i, :))



        ### Dimensionalize
        for j = 1:na
            u_i = cchistory[i].W[j] #TODO: Should I be normalizing by the actual nodal velocity, or the undistrubed nodal velocity. 
            # N[i,j] = Cn[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
            # T[i,j] = Ct[i,j]*0.5*env.rho*u_i^2*chordvec[j]
            Fx[i,j] = Cx[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
            Fy[i,j] = Cy[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
        end



        if verbose & (mod(i, speakiter)==0)
            println("Simulation time: ", t)
        end
    end


    # return (N=N, T=T, Fx=Fx, Fy=Fy), (Cx=Cx, Cy=Cy, Cn=Cn, Ct=Ct, Cl=Cl, Cd=Cd, Cm=Cm), cchistory, xds, azimuth 
    return (Fx=Fx, Fy=Fy), (Cx=Cx, Cy=Cy, Cm=Cm), cchistory, xds, azimuth 
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

    nt, _ = size(loads[1].Fx)

    thrust = Array{eltype(loads[1].Fx)}(undef, nt)
    torque = Array{eltype(loads[1].Fx)}(undef, nt)

    rfull = [rhub; rvec; rtip]
    
    for i = 1:nt
        for j = eachindex(loads)
            Fxfull = [0.0; loads[j].Fx[i,:]; 0.0]
            Fyfull = [0.0; loads[j].Fy[i,:]; 0.0]

            thrust[i] += FLOWMath.trapz(rfull, Fxfull)
            torque[i] += FLOWMath.trapz(rfull, Fyfull.*rfull)
        end
    end

    return thrust, torque
end