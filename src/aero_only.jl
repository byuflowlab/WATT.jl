#=

Recreating what they do in AeroDyn, which is solves the BEM, feeds the inflow angle to the dynamic stall model, then calculates the loading. 

Adam Cardoza 8/6/22
=# 

export simulate

#Question: I don't know if using this many structs will make my package take forever to compile... So it might be better to create a nest of if statements... but we'll see.

function dimensionalize!(Fx, Fy, Mx, Cx, Cy, Cm, blade::Blade, env::Environment, W)
    
    for j in eachindex(blade.r)
        qinf = 0.5*env.rho*W[j]^2
        
        Fx[j] = Cx[j]*qinf*blade.airfoils.c[j] 
        Fy[j] = Cy[j]*qinf*blade.airfoils.c[j] 
        Mx[j] = Cm[j]*qinf*blade.airfoils.c[j]^2 #The coefficient of moment is positive about the negative Z aero axis, so we need the negative of this to move it to the structural X axis. 
    end
end

function update_aerostates!(aerostates::AeroStates, ccout, i, j)

    # aerostates.phi[i,j] = mesh.cchistory[j].phi
    # aerostates.alpha[i,j] = mesh.cchistory[j].alpha
    # aerostates.W[i,j] = mesh.cchistory[j].W

    aerostates.phi[i,j] = ccout.phi
    aerostates.alpha[i,j] = ccout.alpha
    aerostates.W[i,j] = ccout.W
    
end

"""
    take_aero_step!()

Take an in-place step of the aerodynamic models. 

"""
function take_aero_step!(aerostates::AeroStates, mesh::Mesh, rotor::Rotor, blade::Blade, env::Environment, tvec, i, pitch; solver::Solver=RK4(), pfunc = (p,t) -> (;), prepp=nothing, p=nothing)
    na = length(blade.r)

    t = tvec[i]
    dt = tvec[i] - tvec[i-1]
    
    
    if dt<0
        error("Time step is negative")
    end

    #update azimuthal position
    aerostates.azimuth[i] = env.RS(t)*dt + aerostates.azimuth[i-1] #Euler step for azimuthal position. #TODO: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

    if aerostates.azimuth[i]<aerostates.azimuth[i-1]
        @warn("Blade moved backwards")
    end

    ### Update BEM inputs and solve
    for j = 1:na
        ### Update base inflow velocities
        Vx, Vy = Rotors.get_aerostructural_velocities(rotor, blade, env, t, j, aerostates.azimuth[i], mesh.delta[j], mesh.def_theta[j], mesh.aerov[j])
        
        #TODO: Write a solver that is initialized with the previous inflow angle.
        # mesh.cchistory[j] = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)
        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc)

        update_aerostates!(aerostates, ccout, i, j)
        # update_aerostates!(aerostates, mesh, i, j)
    end
    
    ### Update Dynamic Stall model inputs 
    update_ds_inputs!(blade.airfoils, view(mesh.p_ds, :), view(aerostates.W, i, :), view(aerostates.phi, i, :), blade.twist, pitch, dt, rotor.turbine)
    
    ### Integrate Dynamic Stall model
    update_ds_states!(solver, blade.airfoils, view(aerostates.xds, i-1, :), view(aerostates.xds, i, :), mesh.xds_idxs, mesh.p_ds, t, dt)

    ### Extract loads 
    extract_ds_loads!(blade.airfoils, view(aerostates.xds, i, :), mesh.xds_idxs, view(aerostates.phi, i, :), mesh.p_ds, view(aerostates.cx, i, :), view(aerostates.cy, i, :), view(aerostates.cm, i, :))
 
    
    dimensionalize!(view(aerostates.fx, i, :), view(aerostates.fy, i, :), view(aerostates.mx, i, :), view(aerostates.cx, i, :), view(aerostates.cy, i, :), view(aerostates.cm, i, :), blade::Blade, env::Environment, view(aerostates.W, i, :))
    #These loads do not need to be rotated because they will be applied in the deflected frame (a follower load). This should also be true for things like precone, tilt, and yaw if they are defined correctly in GXBeam. 
end

"""
    simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec; turbine::Bool=true, dsmodel::DS.DSModel=DS.riso(blade.airfoils), dsmodelinit::ModelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, azimuth0=0.0)

Simulate the rotor's response for a given rotor and environmental condition. 

### Inputs
- rvec::Array{TF, 1}
- chordvec::Array{TF, 1}

### Outputs


### Notes
"""
function simulate(rotor::Rotors.Rotor, blade::Blade, env::Environment, tvec;
     pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter=100,
     azimuth0=0.0)


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
    na = length(blade.r)
    nt = length(tvec)

    t0 = tvec[1]

    turbine = rotor.turbine

    # rvec = @. sqrt(blade.rx^2 + blade.ry^2 + blade.rz^2)
    twistvec = blade.twist
    
    airfoils = blade.airfoils
    chordvec = airfoils.c


    ### Initialize BEM solution
    Vxvec = Array{eltype(chordvec)}(undef, na)
    Vyvec = Array{eltype(chordvec)}(undef, na)
    azimuth = Array{eltype(chordvec)}(undef, nt)
    azimuth[1] = azimuth0

    cchistory = Array{CCBlade.Outputs{eltype(chordvec)}, 2}(undef, nt, na)
    xcc = zeros(11) #Todo: Typing


    for j = 1:na
        Vxvec[j], Vyvec[j] = get_aero_velocities(rotor, blade, env, t0, j, azimuth[1])

        cchistory[1,j] = solve_BEM!(rotor, blade, env, j, Vxvec[j], Vyvec[j], pitch, xcc)
    end



    ### Initialize DS solution #Todo: Create a get Wdot function (from the environment struct. )
    Wdotvec = zeros(na) #[sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2) for i in 1:na] #TODO: I probably need to update if there is precone, tilt, yaw, etc. -> Maybe I'll make a function to do this. -> Currently not used by BLADG. 
    alphadotvec = zero(Wdotvec) #TODO: 

    #Todo: This is going to need to be rearrange now. 
    xds, xds_idxs, p_ds = initialize_ds_model(airfoils, turbine, nt, tvec, cchistory[1,:], Wdotvec, alphadotvec, twistvec, pitch)





    ### Prepare data storage
    Cx = Array{eltype(chordvec)}(undef,(nt, na))
    Cy = Array{eltype(chordvec)}(undef,(nt, na))
    Cm = Array{eltype(chordvec)}(undef,(nt, na))

    Fx = Array{eltype(chordvec)}(undef,(nt, na))
    Fy = Array{eltype(chordvec)}(undef,(nt, na))


    extract_ds_loads!(airfoils, view(xds, 1, :), xds_idxs, cchistory[1,:], p_ds, view(Cx, 1, :), view(Cy, 1, :), view(Cm, 1, :))
    

    ### Dimensionalize
    for j = 1:na
        u_1 = cchistory[1, j].W
        Fx[1,j] = Cx[1,j]*0.5*env.rho*u_1^2*chordvec[j] 
        Fy[1,j] = Cy[1,j]*0.5*env.rho*u_1^2*chordvec[j] 
    end



    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. 

        ### Update BEM inputs & solve BEM
        for j = 1:na
            Vxvec[j], Vyvec[j] = get_aero_velocities(rotor, blade, env, t, j, azimuth[i])

            cchistory[i, j] = solve_BEM!(rotor, blade, env, j, Vxvec[j], Vyvec[j], pitch, xcc)
        end


        
        ### Update Dynamic Stall model inputs
        update_ds_inputs!(airfoils, p_ds, cchistory[i, :].W, cchistory[i, :].phi, twistvec, pitch, dt, turbine) #Todo: There is probably a more efficient way of passing in W, and Phi. Maybe I keep a running vector of W and phi that I just update. That way I'm not allocating every iteration. Better yet. I wonder if I should just copy out code from CCBlade and get rid of the Outputs struct and just store it in a matrix. 


        ### Integrate Dynamic Stall model
        update_ds_states!(solver, airfoils, view(xds, i-1, :), view(xds, i, :), xds_idxs, p_ds, t, dt)





        ### Extract loads 
        extract_ds_loads!(airfoils, view(xds, i, :), xds_idxs, cchistory[i, :], p_ds, view(Cx, i, :), view(Cy, i, :), view(Cm, i, :))



        ### Dimensionalize
        for j = 1:na
            u_i = cchistory[i, j].W 
            Fx[i,j] = Cx[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
            Fy[i,j] = Cy[i,j]*0.5*env.rho*u_i^2*chordvec[j] 
        end



        if verbose & (mod(i, speakiter)==0)
            println("Simulation time: ", t)
        end
    end #End iterating through time. 


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