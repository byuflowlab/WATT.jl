#=

Recreating what they do in AeroDyn, which is solves the BEM, feeds the inflow angle to the dynamic stall model, then calculates the loading. 

Adam Cardoza 8/6/22
=# 

export simulate

#Question: I don't know if using this many structs will make my package take forever to compile... So it might be better to create a nest of if statements... but we'll see.

function dimensionalize!(Fx, Fy, Mx, Cx, Cy, Cm, blade::Blade, env::Environment, W)
    
    for j in eachindex(blade.r)
        q_local = 0.5*env.rho*W[j]^2 #Local dynamic pressure
        
        Fx[j] = Cx[j]*q_local*blade.c[j] 
        Fy[j] = Cy[j]*q_local*blade.c[j] 
        Mx[j] = Cm[j]*q_local*blade.c[j]^2 #The coefficient of moment is positive about the negative Z aero axis, so we need the negative of this to move it to the structural X axis. 
    end
end

# function update_aerostates!(aerostates::AeroStates, ccout, i, j)

#     # aerostates.phi[i,j] = mesh.cchistory[j].phi
#     # aerostates.alpha[i,j] = mesh.cchistory[j].alpha
#     # aerostates.W[i,j] = mesh.cchistory[j].W

#     aerostates.phi[i,j] = ccout.phi
#     aerostates.alpha[i,j] = ccout.alpha
#     aerostates.W[i,j] = ccout.W
    
# end


"""
    initialize()

Prepare data structures for an aerodynamic only simulation. 

**Arguments**
- blade::Blade - The blade object that contains the airfoils, twist, and chord.
- tvec::Vector{TF} - The time vector that the simulation will be run over.

**Outputs**
- aerostates::AeroStates - The aerodynamic states that are calculated during the simulation.
- mesh::Mesh - The mesh that is used to store the structural deflections and velocities.
"""
function initialize(blade::Blade, tvec; verbose::Bool=false, inittype=nothing)
    #Todo: This still needs to be completed. And it'll need a initial condition function. But for now, it looks like the simulate function initializes, finds the initial condition, and then simulates. 

    if verbose
        println("Rotors.jl initializing solution...")
    end

    
    # if warnings
    #     checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    # end

    
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types.

    ### Initialization information
    na = length(blade.rR)
    nt = length(tvec)

    t0 = first(tvec)

    if isnothing(inittype)
        inittype = find_inittype(blade.airfoils[1].c, blade.twist[1])
    end #TODO: I should probably check if the passed in type is a valid type. 


    ### ----- Prepare data storage for aerodynamic models ----- ###
    azimuth = Array{inittype}(undef, nt)
    phi = Array{inittype}(undef,(nt, na))
    alpha = Array{inittype}(undef,(nt, na))
    W = Array{inittype}(undef,(nt, na))

    Cx = Array{inittype}(undef,(nt, na))
    Cy = Array{inittype}(undef,(nt, na))
    Cm = Array{inittype}(undef,(nt, na))

    Fx = Array{inittype}(undef,(nt, na))
    Fy = Array{inittype}(undef,(nt, na))
    Mx = Array{inittype}(undef,(nt, na))

    # A vector that CCBlade uses for solving. 
    xcc = Vector{inittype}(undef, 11)

    # Initialize DS solution
    xds, xds_idxs, y_ds = initialize_ds_model(blade.airfoils, nt; inittype)  

    # Store everything in the aerostates 
    aerostates = (;azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds)


    

    #Placeholders for strucutral deflections
    delta = Vector{SVector{3, inittype}}(undef, na) #todo. What is this used for? -> This is the deflection from the structural mesh projected onto the aero mesh. We keep it as zeros for aero-only simulations.
    def_theta = Vector{SVector{3, inittype}}(undef, na) #todo. What is this used for
    #The structural velocities interpolated to the aerodynamic nodes.
    aerov = Vector{SVector{3, inittype}}(undef, na)

    for i = 1:na
        delta[i] = SVector{3, inittype}(0.0, 0.0, 0.0)
        def_theta[i] = SVector{3, inittype}(0.0, 0.0, 0.0)
        aerov[i] = SVector{3, inittype}(0.0, 0.0, 0.0)
    end

    
    mesh = (; delta, def_theta, aerov, xcc, xds_idxs, y_ds)

    return aerostates, mesh
end

"""
    take_aero_step!()

Take an in-place step of the aerodynamic models. 

**Arguments**
- phi::Vector{TF} - The inflow angle at every aerodynamic node. Calculating inplace. 
- alpha::Vector{TF} - The angle of attack at every aerodynamic node. Calculating inplace.
- W::Vector{TF} - The inflow velocity at every aerodynamic node. Calculating inplace. 
- xds::Vector{TF} - The new dynamic stall states at every aerodynamic node. Calculating inplace. 
- cx::Vector{TF} - The coefficient of force in the x direction at every aerodynamic node. Calculating inplace. 
- cy
- cm
- fx
- fy
- mx
- xds_old::Vector{TF} - The old dynamic stall states (referenced to calculate new states). 

"""
function take_aero_step!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, xds_old, azimuth, t, dt, pitch, mesh, rotor::Rotor, blade::Blade, env::Environment; solver::Solver=RK4())
    #TODO: I don't think that I need pfunc, nor prepp, nor p here. 
    #TODO: I'm thinking that I won't have any optional arguments here. It doesn't seem needed... Well... maybe it'd be handy when I'm using the function outside of the package for optimization. 

    # println("Using this function")

    na = length(blade.r)
    
    
    if dt<0
        error("Time step is negative")
    end

    #update azimuthal position
    # azimuth = env.RS(t)*dt + azimuth0 #Euler step for azimuthal position. #TODO: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

    

    ### Update BEM inputs and solve
    for j = 1:na
        ### Update base inflow velocities
        # @show mesh.aerov[j][1].value, mesh.aerov[j][2].value, mesh.aerov[j][2].value
        # @show mesh.aerov[j]
        Vx, Vy = Rotors.get_aerostructural_velocities(rotor, blade, env, t, j, azimuth, mesh.delta[j], mesh.def_theta[j], mesh.aerov[j])
        #Todo: Angles aren't updated with angular deflection.... 
        
        #TODO: Write a solver that is initialized with the previous inflow angle.
        
        # ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc) 
        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch - mesh.def_theta[j][1], mesh.xcc) 
        # ccout = solve_BEM!(rotor, blade, env, phi_old[j], j, Vx, Vy, pitch, mesh.xcc) #Todo: Need to create some sort of fail safe for not converging. 

        phi[j] = ccout.phi
        alpha[j] = ccout.alpha
        W[j] = ccout.W
    end
    
    ### Update Dynamic Stall model inputs 
    update_ds_inputs!(blade.airfoils, view(mesh.y_ds, :), W, phi, blade.twist, pitch, dt, rotor.turbine, blade)
    
    ### Integrate Dynamic Stall model
    update_ds_states!(solver, blade.airfoils, xds_old, xds, mesh.xds_idxs, mesh.y_ds, mesh.p_ds, t, dt)

    ### Extract loads 
    extract_ds_loads!(blade.airfoils, xds, mesh.xds_idxs, phi, mesh.y_ds, mesh.p_ds, cx, cy, cm)
 
    
    dimensionalize!(fx, fy, mx, cx, cy, cm, blade::Blade, env::Environment, W)
    #These loads do not need to be rotated because they will be applied in the deflected frame (a follower load). This should also be true for things like precone, tilt, and yaw if they are defined correctly in GXBeam. 

end




"""
    simulate!(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade::Blade, env::Environment, tvec; turbine::Bool=true, dsmodel::DS.DSModel=DS.riso(blade.airfoils), dsmodelinit::ModelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, azimuth0=0.0)

Simulate the rotor's response for a given rotor and environmental condition. 

### Inputs
- rvec::Array{TF, 1}
- chordvec::Array{TF, 1}

### Outputs


### Notes
"""
function simulate!(aerostates, mesh, rotor::Rotors.Rotor, blade::Blade, env::Environment, tvec;
     pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter=100,
     azimuth0=0.0)


    if verbose
        println("Rotors.jl finding initial solution...")
    end

    # if isnothing(inittype)
    #     inittype = find_inittype(blade.airfoils[1].c, blade.twist[1])
    # end #TODO: I should probably check if the passed in type is a valid type. 

    #TODO: This will break without dsmodel exposed. 
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

    @unpack azimuth, phi, alpha, W, Cx, Cy, Cm, Fx, Fy, Mx, xds = aerostates

    @unpack y_ds, xds_idxs = mesh



    # turbine = rotor.turbine #Flag: Is this a turbine or a propeller?

    # rvec = @. sqrt(blade.rx^2 + blade.ry^2 + blade.rz^2)
    # twistvec = blade.twist
    
    airfoils = blade.airfoils
    # chordvec = airfoils.c


    ### Initial Condition
    azimuth[1] = azimuth0

    phi0 = view(phi, 1, :)
    alpha0 = view(alpha, 1, :)
    W0 = view(W, 1, :)

    ### Initialize BEM solution
    for j = 1:na #TODO. Write a solver that is initialized with the previous inflow angle. -> Doesn't drastically change much.
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth[1])

        ccout = solve_BEM!(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc) #TODO: Does mesh.xcc allocate a new vector? 
        phi0[j] = ccout.phi
        alpha0[j] = ccout.alpha
        W0[j] = ccout.W
    end


    ### Initialize DS solution 
    xds0 = view(xds, 1, :)
    dsmodel_initial_condition!(xds0, phi0, W0, mesh, blade, rotor.turbine, t0, pitch) #

    Cx0 = view(Cx, 1, :)
    Cy0 = view(Cy, 1, :)
    Cm0 = view(Cm, 1, :)
    extract_ds_loads!(airfoils, xds0, xds_idxs, phi0, y_ds, Cx0, Cy0, Cm0)

    Fx0 = view(Fx, 1, :)
    Fy0 = view(Fy, 1, :)
    Mx0 = view(Mx, 1, :)

    dimensionalize!(Fx0, Fy0, Mx0, Cx0, Cy0, Cm0, blade, env, W0) 


    if verbose
        println("Rotors.jl starting simulation...")
    end
    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        if dt<0
            error("Time step is negative")
        end

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. 

        if azimuth[i]<azimuth[i-1]
            @warn("Blade moved backwards")
        end

        ### Unpack
        phi_i = view(phi, i, :)
        alpha_i = view(alpha, i, :)
        W_i = view(W, i, :)
        cx_i = view(Cx, i, :)
        cy_i = view(Cy, i, :)
        cm_i = view(Cm, i, :)
        fx_i = view(Fx, i, :)
        fy_i = view(Fy, i, :)
        mx_i = view(Mx, i, :)
        xds_i = view(xds, i, :)
        xds_im1 = view(xds, i-1, :)

        take_aero_step!(phi_i, alpha_i, W_i, xds_i, cx_i, cy_i, cm_i, fx_i, fy_i, mx_i, xds_im1, azimuth[i], t, dt, pitch, mesh, rotor, blade, env; solver)



        if verbose & (mod(i, speakiter)==0)
            println("Simulation time: ", t)
        end
    end #End iterating through time. 

    if verbose
        println("Rotors.jl simulation complete.")
    end
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