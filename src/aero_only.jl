#=

Recreating what they do in AeroDyn, which is solves the BEM, feeds the inflow angle to the dynamic stall model, then calculates the loading. 

Adam Cardoza 8/6/22
=# 

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
    initialize_aero(blade, tvec; inittype=nothing, verbose=false) -> aerostates, mesh

Pre-allocate buffers for an aerodynamics-only transient simulation (no
structural feedback). The returned pair is the canonical input to
[`simulate!`](@ref).

**Arguments**
- `blade::Blade`
- `tvec::AbstractVector`: Time vector.

**Keyword Arguments**
- `inittype`: Override the inferred element type (default: from `blade.c[1]`,
  `blade.twist[1]`).
- `verbose::Bool = false`

**Returns**
- `aerostates::AeroStates`: Time-indexed aero state history.
- `mesh::AeroMesh`: BEM/DS scratch buffers. `delta`/`def_theta`/`aerov` are
  allocated for shape compatibility with the coupled solver but stay zero.
"""
function initialize_aero(blade::Blade, tvec; verbose::Bool=false, inittype=nothing)
    #Todo: This still needs to be completed. And it'll need a initial condition function. But for now, it looks like the simulate function initializes, finds the initial condition, and then simulates. 

    if verbose
        println("WATT.jl initializing solution...")
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
        inittype = find_inittype(blade.c[1], blade.twist[1])
    end #TODO: I should probably check if the passed in type is a valid type. 


    ### ----- Prepare data storage for aerodynamic models ----- ###

    # Initialize DS state vector first so we know its width (ns).
    xds, xds_idxs, y_ds, p_ds = initialize_ds_model(blade, nt, inittype)
    ns = size(xds, 2)

    aerostates = AeroStates{inittype}(undef, nt, na, ns)
    copyto!(aerostates.xds, xds)
    xds = aerostates.xds

    xcc = Vector{inittype}(undef, 11)

    # Placeholders for structural quantities; stay zero in aero-only mode.
    delta = Vector{SVector{3, inittype}}(undef, na)
    def_theta = Vector{SVector{3, inittype}}(undef, na)
    aerov = Vector{SVector{3, inittype}}(undef, na)
    for i = 1:na
        delta[i] = SVector{3, inittype}(0.0, 0.0, 0.0)
        def_theta[i] = SVector{3, inittype}(0.0, 0.0, 0.0)
        aerov[i] = SVector{3, inittype}(0.0, 0.0, 0.0)
    end

    mesh = AeroMesh(delta, def_theta, aerov, xcc, xds_idxs, y_ds, p_ds)

    return aerostates, mesh
end

"""
    take_aero_step!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, xds_old, azimuth, t, dt, pitch, mesh, rotor, blade, env; solver=RK4())

Advance the aerodynamic state (BEM + DS) by one time step `dt`. Writes the
new aero outputs into the supplied views — does not allocate.

**Arguments**
- `phi, alpha, W::AbstractVector`: BEM outputs at every node (written).
- `xds::AbstractVector`: New DS states (written).
- `cx, cy, cm`: Force/moment coefficients (written).
- `fx, fy, mx`: Dimensional sectional loads (written).
- `xds_old::AbstractVector`: DS states at the previous step.
- `azimuth, t, dt, pitch::Real`
- `mesh::AbstractSimMesh`, `rotor::Rotor`, `blade::Blade`, `env::Environment`

**Keyword Arguments**
- `solver::Solver = RK4()`: DS state integrator.
"""
function take_aero_step!(phi, alpha, W, xds, cx, cy, cm, fx, fy, mx, xds_old, azimuth, t, dt, pitch, mesh, rotor::Rotor, blade::Blade, env::Environment; solver::Solver=RK4())

    na = length(blade.r)
    
    if dt<0
        error("Time step is negative")
    end

    ### Update BEM inputs and solve
    for j = 1:na
        ### Update base inflow velocities
        Vx, Vy = get_aerostructural_velocities(rotor, blade, env, t, j, azimuth, mesh.delta[j], mesh.def_theta[j], mesh.aerov[j])
        
        ccout = solve_BEMT(rotor, blade, env, j, Vx, Vy, pitch - mesh.def_theta[j][1], mesh.xcc) #Correct twist based on the structural deformation. 

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
    simulate!(aerostates, mesh, rotor, blade, env, tvec; pitch=0.0, solver=RK4(), kwargs...)

Pre-allocated aerodynamics-only transient simulation. Mutates `aerostates`
and the mutable fields of `mesh` in place.

**Arguments**
- `aerostates::AeroStates`, `mesh::AeroMesh`: From [`initialize_aero`](@ref).
- `rotor::Rotor`, `blade::Blade`, `env::Environment`
- `tvec::AbstractVector`: Time vector.

**Keyword Arguments**
- `pitch::Real = 0.0`: Blade pitch (rad).
- `solver::Solver = RK4()`: DS state integrator.
- `azimuth0::Real = 0.0`: Initial azimuth.
- `verbose::Bool = false`, `speakiter::Int = 100`: Progress printing.
"""
function simulate!(aerostates, mesh, rotor::Rotor, blade::Blade, env::Environment, tvec;
     pitch=0.0, solver::Solver=RK4(), verbose::Bool=false, speakiter=100,
     azimuth0=0.0)


    if verbose
        println("WATT.jl finding initial solution...")
    end

    # if isnothing(inittype)
    #     inittype = find_inittype(blade.airfoils[1].c, blade.twist[1])
    # end #TODO: I should probably check if the passed in type is a valid type. 

    #TODO: This will break without dsmodel exposed. 
    # if isa(dsmodel.detype, DS.Functional)
    #     error("WATT isn't set up to simulate with functional forms of the dsmodel yet, choose the iterative model.")
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

    @unpack y_ds, p_ds, xds_idxs = mesh

    airfoils = blade.airfoils

    ### Initial Condition
    azimuth[1] = azimuth0

    phi0 = view(phi, 1, :)
    alpha0 = view(alpha, 1, :)
    W0 = view(W, 1, :)

    ### Initialize BEM solution
    for j = 1:na 
        Vx, Vy = get_aero_velocities(rotor, blade, env, t0, j, azimuth[1])

        ccout = solve_BEMT(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc) #TODO: Does mesh.xcc allocate a new vector? 
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
    extract_ds_loads!(airfoils, xds0, xds_idxs, phi0, y_ds, p_ds, Cx0, Cy0, Cm0)

    Fx0 = view(Fx, 1, :)
    Fy0 = view(Fy, 1, :)
    Mx0 = view(Mx, 1, :)

    dimensionalize!(Fx0, Fy0, Mx0, Cx0, Cy0, Cm0, blade, env, W0) 


    if verbose
        println("WATT.jl starting simulation...")
    end

    ### Iterate through time #Note: It would be beneficial to save AD across a single time step and reuse (if reverse mode), but that only works if the dynamic stall model doesn't have any branching code, which it does. 
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
        println("WATT.jl simulation complete.")
    end
end


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