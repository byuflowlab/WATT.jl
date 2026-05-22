#=
A fixed point solver that iterates back and forth between CCBlade and GXBeam until it converges on a solution. A dynamic stall model isn't included, so the solver assumes that the blade is at some sort of non-oscillitory steady state. 

7/9/24 Adam Cardoza
=#



"""
    initialize_static(blade, assembly; structural_damping=true, linear=false, kwargs...) -> aerostates, gxstates, mesh

Pre-allocate buffers for the steady-state fixed-point aero-structural solver
[`fixedpoint!`](@ref). No time dimension and no dynamic-stall state.

**Arguments**
- `blade::Blade`
- `assembly::GXBeam.Assembly`

**Returns**
- `aerostates::StaticAeroStates`: Per-node BEM/load outputs at the equilibrium.
- `gxstates::Vector`: GXBeam state buffer of length `2*length(system.x)`.
- `mesh::StaticMesh`: Coupling buffers.
"""
function initialize_static(blade::Blade, assembly::GXBeam.Assembly; verbose::Bool=false, p=nothing, pfunc = (p,t) -> (;), xpfunc=nothing, structural_damping::Bool=true, linear::Bool=false)

    if verbose
        println("WATT.jl initializing solution...")
    end

    # if warnings
    #     checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    # end

    
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types.

    ### Initialization information
    na = length(blade.rR)

    inittype = find_inittype(blade.c[1], blade.twist[1])



    ### ----- Prepare data storage for aerodynamic models ----- ###
    aerostates = StaticAeroStates{inittype}(undef, na)

    # A vector that CCBlade uses for solving.
    xcc = Vector{inittype}(undef, 11)

    

    ### ----- Allocate the GXBeam Data ----- ###

    system = GXBeam.DynamicSystem(assembly)

    ngx = length(system.x)

    gxstates = Array{inittype}(undef, 2*ngx)

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{Float64}}()

    # Allocate the prescribed conditions 
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)) # root section is fixed

    
    

    # The currently un-used GXBeam constants
    point_masses=Dict{Int,PointMass{Float64}}()
    linear_velocity=(@SVector zeros(3))
    angular_velocity=(@SVector zeros(3))
    # structural_damping=true #todo: I might add this in later. 
    two_dimensional=false
    xpfunc = nothing
    



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 


    delta = Vector{SVector{3, inittype}}(undef, na) #todo: What is this used for? 
    def_theta = Vector{SVector{3, inittype}}(undef, na) #todo. What is this used for
    #The structural velocities interpolated to the aerodynamic nodes.
    aerov = Vector{SVector{3, inittype}}(undef, na) #Todo: How are these suppose to be initialized? Like their first value? 

    for i = 1:na
        delta[i] = SVector(0.0, 0.0, 0.0)
        def_theta[i] = SVector(0.0, 0.0, 0.0)
        aerov[i] = SVector(0.0, 0.0, 0.0)
    end

    
    mesh = StaticMesh(interpolationpoints, delta, def_theta, aerov, xcc,
                assembly, system, prescribed_conditions, distributed_loads,
                point_masses, linear_velocity, angular_velocity,
                xpfunc, pfunc, two_dimensional, structural_damping, linear)

    return aerostates, gxstates, mesh
end



"""
    fixedpoint!(aerostates, gxstates, azimuth0, rotor, env, blade, mesh, pitch; iterations=1, kwargs...) -> aerostates, gxstates, mesh

Solve the steady-state aero-structural fixed point by partitioned iteration
between CCBlade and `GXBeam.steady_state_analysis!`. No dynamic-stall model.

**Arguments**
- `aerostates::StaticAeroStates`, `gxstates::Vector`, `mesh::StaticMesh`: From
  [`initialize_static`](@ref).
- `azimuth0::Real`: Fixed blade azimuth (rad).
- `rotor::Rotor`, `env::Environment`, `blade::Blade`
- `pitch::Real`: Blade pitch (rad).

**Keyword Arguments**
- `iterations::Int = 1`: Number of partitioned aero ↔ structural sweeps.
- `atol::Real = 1e-3`: Reserved for a future convergence-tracking exit.
- `g::Real = 9.81`: Gravity magnitude.
- `solver::Solver = RK4()`, `prepp`, `p`, `pfunc`, `gxflag`, `verbose`

**Notes**
The first iteration matches a direct CCBlade call at the undeflected
position; convergence is empirically near machine precision by ~10 iterations
on NREL 5MW at rated conditions.
"""
function fixedpoint!(aerostates, gxstates, azimuth0, rotor::Rotor, env::Environment, blade, mesh, pitch; iterations=1, atol=1e-3, g=9.81, gxflag=nothing, verbose::Bool=false, solver::Solver=RK4(), prepp=nothing, p=nothing, pfunc=nothing) #Todo: Maybe add in something to cut off the solution if it acheives a given tolerance...-> requires some extra memory? Unless maybe I can use the mesh or something? ... Naw, I think it'll require extra memory. 
    #TODO: Maybe include the coupled aerostructural solution as an option? 

    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. -> These might best go in the blade constructor. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). -> I'm not entirely sure what I meant by that. 

    @unpack assembly, system, structural_damping, linear, pfunc, distributed_loads, prescribed_conditions = mesh

    @unpack phi, alpha, W, Cx, Cy, Cm, Cx, Cy, Fx, Fy, Mx = aerostates

    na = length(blade.r)

    verbose ? println("Solving fixed point static solution...") : nothing

    initial_condition_checks(gxflag)

    angular_velocity = SVector(0.0, 0.0, -env.RS(0.0)) 
    gravity = SVector(-g*cos(azimuth0), -g*sin(azimuth0), 0.0)

    airfoils = blade.airfoils

    for i = 1:iterations
        ###BEMT solution 
        for j = 1:na
            ### Update base inflow velocities
            Vx, Vy = get_aerostructural_velocities(rotor, blade, env, 0.0, j, azimuth0, mesh.delta[j], mesh.def_theta[j], mesh.aerov[j])
            
            ccout = solve_BEMT(rotor, blade, env, j, Vx, Vy, pitch, mesh.xcc) 
            # ccout = solve_BEMT(rotor, blade, env, phi_old[j], j, Vx, Vy, pitch, mesh.xcc) #todo: Need to create some sort of fail safe for not converging. 

            phi[j] = ccout.phi
            alpha[j] = ccout.alpha
            W[j] = ccout.W

            sphi, cphi = sincos(phi[j])
            Cx[j] = ccout.cl*cphi + ccout.cd*sphi
            Cy[j] =  -(ccout.cl*sphi - ccout.cd*cphi)
            Cm[j] = 0.0 #airfoils[j].cm(alpha[j]) #todo: Is this applied correctly? 
        end

        #Note: I don't use take_aero_step because these functions are different. 
        dimensionalize!(Fx, Fy, Mx, Cx, Cy, Cm, blade, env, W) 
        update_forces!(distributed_loads, Fx, Fy, Mx, blade, assembly)


        ### GXBeam initial solution 
        if isnothing(prepp)
            p = nothing
        else
            prepp(p, Fx, Fy, Mx)
        end

        system, gxstates, converged = GXBeam.steady_state_analysis!(system, assembly; prescribed_conditions, distributed_loads, angular_velocity, gravity, pfunc, p, linear)  

        ### Update mesh transfer variables
        update_mesh!(blade, mesh, assembly, gxstates, env, 0., na)
    end
    

    return aerostates, gxstates, mesh #todo: I don't know if gxstate will update as I pass in. -> It appears to work, but I should double check. 
end