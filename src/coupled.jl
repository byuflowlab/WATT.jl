#=


=#







"""
    find_paired_points(points, r)

Finds the two structural points before and after the given aerodynamic node. 

### Inputs
- points - the points vector from the assembly object. Each point is a SVector{3}. 
- r - the radial location of the aerodynamic point of interest

### Outputs
- pair - a tuple containing the indices of the point before and after the aerodynamic point (before, after)
"""
function find_paired_points(points, r)
    ngx = length(points)
    rgx = [points[i][1] for i = 1:ngx]
    for i = 1:length(points)-1
        if rgx[i]>r
            return (i-1, i)
        end
    end
    return (ngx-1, ngx)
end

function find_interpolation_value(points, pair, r)
    ngx = length(points)
    rgx = [points[i][1] for i = 1:ngx]
    a = r - rgx[pair[1]]
    b = rgx[pair[2]] - r
    L = a + b
    return a/L
end

"""
    InterpolationPoint(pair, percent)

A struct to quickly interpolate the new radial location of an aerodynamic node. 
"""
struct InterpolationPoint
    pair::Tuple{Int64, Int64}
    percent::Float64
end

"""
    create_interpolationpoints(assembly, rvec)

### Inputs
- assembly - A GXBeam assembly
- rvec - The radial location of the aerodynamic nodes

### Outputs
- interpvals

### Notes
- This assumes that the structural points are in order and go from root to tip. Additionally, the aerodynamic points follow a similar order. 
"""
function create_interpolationpoints(assembly, rvec)
    na = length(rvec)
    pairs = [find_paired_points(assembly.points, rvec[i]) for i = 1:na]
    percents = [find_interpolation_value(assembly.points, pairs[i], rvec[i]) for i =1:na]

    return [InterpolationPoint(pairs[i], percents[i]) for i = 1:na]
end


"""
    interpolate_position(ip, assembly, state)

Get the position of aerodynamic nodes from the interpolation struct and the structural structs. 
"""
function interpolate_position(ip, assembly, state)  
    rgx = [assembly.points[i][1] + state.points[i].u[1] for i = 1:length(assembly.points)]
    return (1-ip.percent)*rgx[ip.pair[1]] + ip.percent*rgx[ip.pair[2]]
end

"""
    interpolate_velocity(ip, assembly, state, env)

Get the relative velocity of an aerodynamic point. 
"""
function interpolate_velocity(ip, assembly, state) #Note: This is probably a very poor approximation. But it'll have to do for now. Maybe there is a better way to get the velocity of the points. -> Taylor said that a linear interpolation should work fine. 
    ne = length(assembly.elements)
    np = ne + 1

    V = state.elements[ip.pair[1]].V  #Get the velocity of the element 
    Omega = state.elements[ip.pair[1]].Omega #Question: If I'm doing what I'm currently doing.... why don't I just calculate the distance from the element node to the aerodynamic node? 
    p1 = assembly.points[ip.pair[1]] #TODO: Should I include the deflected point location? 
    p2 = assembly.points[ip.pair[2]]
    e = assembly.elements[ip.pair[1]].x

    r1 = p1 - e #The position vector from the element point to the starting point
    r2 = p2 - e #The position vector from the element point to the stopping point

    # ra = assembly.elements[ip.pair[1]].L*ip.percent + assembly.points[ip.pair[1]][1]
    
    v1 = V + cross(Omega, r1)
    v2 = V + cross(Omega, r2)

    return (1-ip.percent)*v1 + ip.percent*v2
end

function interpolate_angle(ip, assembly, state) #Todo: I need to see if I should be interpolating or doing what is being done here. 
    # thetagx = [state.points[i].theta[1] for i = 1:length(assembly.points)]
    # return (1-ip.percent)*thetagx[ip.pair[1]] + ip.percent*thetagx[ip.pair[2]]

    return state.elements[ip.pair[1]].theta[1] #Assume that the twist at the closest elemental node is the twist at the aerodynamic node. -> Didn't make much of a difference.
end


function update_structural_forces!(distributed_loads, assembly, gxstate, N, T, env, rvec, rgx, t, b, rhub, rtip)
    Fzfit = Akima(vcat(rhub, rvec, rtip), vcat(0, -N, 0) ) #I use the loads of the previous time. We want everything to be calculated based off of the previous time step to update the next time step. Then we'll move to that step.  
    Fyfit = Akima(rvec, T)

    f_follower = @SVector zeros(3)
    m = @SVector zeros(3)
    m_follower = @SVector zeros(3)
    
    for ielem = 1:length(rgx)
        V = gxstate.elements[ielem].V
        L = assembly.elements[ielem].L
        r_deflected = rgx[ielem] + gxstate.elements[ielem].u[1]
        Vy = V[2] + env.RS(t)*r_deflected #Need to take the rotation velocity out of the structural velocity.
        
        fx = -b*V[1]/L #structural damping in the x direction.
        fy = Fyfit(rgx[ielem]) - b*Vy/L #Aerodynamic load minus the damping #Dividing the damping force by length so that it is a distributed load. 
        fz = Fzfit(rgx[ielem]) - b*V[3]/L #Structural damping in the z direction. 
        # if ielem==length(rgx)
        #     println("y damping: ",  b*Vy/L)
        # end

        f = SVector(fx, fy, fz)
        distributed_loads[ielem] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end
end


function get_aerostructural_velocities(env::Environment, aeroV, t, r, azimuth, precone, tilt, yaw, hubht)

    ### Extract the aero velocities due to the environment, precone, tilt, taw, and shear. 
    vxenv, vyenv = get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

    ### Get the relative structural velocities
    vxrel = aeroV[3] 
    vyrel = aeroV[2] + env.RS(t)*r  

    ### Return the total x and y velocities (x and y in the aerodynamic frame)
    # If the structure is moving in the negative y direction, then the airfoil should see a faster velocity. 
    return vxenv+vxrel, vyenv-vyrel
end


function checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    if minimum(rvec)<=rhub
        @warn("A member of rvec is less than or equal to rhub. This will cause problems with the hub corrections.")
    end
    if maximum(rvec)>=rtip
        @warn("A member of rvec is greater than or equal to rtip. This will cause problems with the tip corrections.")
    end
    if pitch>(pi/4)
        @warn("The pitch value is greater than 45 degrees. Ensure that the pitch is in radians and not degrees.")
    end
    if precone>(pi/4)
        @warn("The precone value is greater than 45 degrees. Ensure that the precone is in radians and not degrees.")
    end
    if tilt>(pi/4)
        @warn("The tilt value is greater than 45 degrees. Ensure that the tilt is in radians and not degrees.")
    end
    if yaw>(pi/4)
        @warn("The yaw value is greater than 45 degrees. Ensure that the yaw is in radians and not degrees.")
    end
    if maximum(twistvec)>(pi/2)
        @warn("The maximum twist value is greater than 90 degrees. Ensure that the twist distribution is given in radians.")
    end
end

#TODO: It might be a good idea to make a version that is completely in place. (Pass in data storage and riso ode.)

#TODO: Function headers
function simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, precone, tilt, yaw, blade::Blade, env::Environment, assembly::GXBeam.Assembly, tvec; turbine::Bool=true, dsmodelinit::DSmodelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false, speakiter=100, warnings::Bool=true, azimuth0=0.0, b=0.01)

    if verbose
        println("Rotors.jl initializing solution...")
    end

    if warnings
        checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
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

    ### Prepare data storage
    Vxvec = Array{eltype(chordvec)}(undef, na)
    Vyvec = Array{eltype(chordvec)}(undef, na)
    azimuth = Array{eltype(chordvec)}(undef, nt)
    Cl = Array{eltype(rvec)}(undef,(nt, na))
    Cd = Array{eltype(rvec)}(undef,(nt, na))
    Cn = Array{eltype(rvec)}(undef,(nt, na))
    Ct = Array{eltype(rvec)}(undef,(nt, na))
    N = Array{eltype(rvec)}(undef,(nt, na))
    T = Array{eltype(rvec)}(undef,(nt, na))
    Fx = Array{eltype(rvec)}(undef,(nt, na))
    Fy = Array{eltype(rvec)}(undef,(nt, na))

    cchistory = Array{Array{CCBlade.Outputs{eltype(rvec)}, 1}, 1}(undef, nt) #[cchistory[1] for i = 1:nt]
    gxhistory = Array{GXBeam.AssemblyState{eltype(rvec), Vector{GXBeam.PointState{eltype(rvec)}}, Vector{GXBeam.ElementState{eltype(rvec)}}}}(undef, nt)




    ### Initialize BEM solution
    azimuth[1] = azimuth0

    rotor = Rotor(rhub, rtip, B; precone, turbine)
    sections = [CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i = 1:na]


    for i = 1:na
        Vxvec[i], Vyvec[i] = get_aero_velocities(env, t0, rvec[i], azimuth[1], precone, tilt, yaw, hubht) #Note: I'm not sure that I need to keep the Vxvec and Vyvec, I could probably just update the operating points inside this vector and initialize the operating points at the beginning. 
    end
    
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na]

    cchistory[1] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)






    ### Initialize DS solution
    risoode = createrisoode(blade)

    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2) for i in 1:na] #TODO: I probably need to update if there is precone, tilt, yaw, etc. -> Maybe I'll make a function to do this. 

    p_ds = vcat(chordvec, twistvec, cchistory[1].phi, cchistory[1].W, Wdotvec, pitch) #p = [chordvec, twistvec, phivec, Wvec, Wdotvec, pitch]

    xds = initializeDSmodel(dsmodelinit, nt, na, p_ds, risoode)




    N[1,:], T[1,:], Cn[1,:], Ct[1,:], Cl[1,:], Cd[1,:] = extractloads(xds[1,:], cchistory[1], t0, rvec, chordvec, twistvec, pitch, blade, env)




    ### Extract CCBlade Loads and create a distributed load #TODO: I feel like there's probably a faster way to do this. Maybe have N include the hub and the tip and CCBlade just updates the nodes inbetween. -> Or maybe there is a better fit function where I get to declare the behavior outside the bounds. 
    Fzfit = Akima(vcat(rhub, rvec, rtip), vcat(0, -N[1,:].*cos(precone), 0))  
    Fyfit = Akima(vcat(rhub, rvec, rtip), vcat(0, T[1,:].*cos(precone), 0))

    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(rvec)}}()
    f_follower = @SVector zeros(3)
    m = @SVector zeros(3)
    m_follower = @SVector zeros(3)

    for ielem = 1:nelem #Iterate through the elements and apply the distributed load at every element.  
        f = SVector(0.0, Fyfit(rgx[ielem]), Fzfit(rgx[ielem])) 
        distributed_loads[ielem] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end 



    ### Prepare GXBeam inputs  
    Omega = SVector(0.0, 0.0, -env.RS(t0))
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed



    ### TODO: I'm still debating whether or not I should include precone, tilt, and yaw in the assembly... what benefits and drawbacks would there be? What things would I have to change if I did? What things should be here if I didn't? 

    #TODO: Create a function to initialize the behavior. One for starting from no loading, one from steady state as below. Look at that, I wrote that down already. 

    ### GXBeam initial solution #TODO: I might make it an option to initialize from rest, or from steady state at the intial conditions (as current).
    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega) 

    gxhistory[1] = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)





    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, rvec) 
    aeroV = [interpolate_velocity(interpolationpoints[j], assembly, gxhistory[1]) for j = 1:na]  
    def_thetax = [interpolate_angle(interpolationpoints[j], assembly, gxhistory[1]) for j = 1:na] 




    # @show Vxvec[end], Vyvec[end]
    # @show maximum(N[1,:]), maximum(T[1,:])


    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. #TODO: Maybe do a better integration like a RK4 or something. 

        ### Update BEM inputs
        for j = 1:na
            ### Update sections with twist
            twist_displaced = twistvec[j] #- def_thetax[j] #Todo: Is this correct? I think this is what I had in the fixed point solution. 
            sections[j] = CCBlade.Section(rvec[j], chordvec[j], twist_displaced, blade.airfoils[j]) #Not displacing rvec because CCBlade uses that to calculate the tip and hub corrections, and the value should be based on relative to the distance along the blade to the hub or tip. (Although in dynamic movement, the tip losses would change dramatically.)

            ### Update base inflow velocities
            Vxvec[j], Vyvec[j] = get_aerostructural_velocities(env, aeroV[j], t, rvec[j], azimuth[i], precone, tilt, yaw, hubht)
            operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
        end

        ### Solve BEM
        cchistory[i] = CCBlade.solve.(Ref(rotor), sections, operatingpoints) #TODO: Write a solver that is initialized with the previous inflow angle. 


        ### Update Dynamic Stall model inputs #TODO: This can probably be made into a function. 
        p_ds[2*na+1:3*na] = cchistory[i].phi #Update phi
        p_ds[3*na+1:4*na] = cchistory[i].W #Update the inflow velocity
        for j = 1:na
            idx = 4*na + j
            p_ds[idx] = sqrt(env.Vinfdot(t)^2 + (env.RSdot(t)*rvec[j])^2) #Update Wdotvec
        end



        ### Integrate Dynamic Stall model
        xds[i,:] = solver(risoode, xds[i-1,:], p_ds, t, dt)



        ### Extract loads
        N[i,:], T[i,:], Cn[i,:], Ct[i,:], Cl[i,:], Cd[i,:] = extractloads(xds[i,:], cchistory[i], t, rvec, chordvec, twistvec, pitch, blade, env)
        
        





        ### Update GXBeam loads 
        update_structural_forces!(distributed_loads, assembly, gxhistory[i-1], N[i,:], T[i,:], env, rvec, rgx, t, b, rhub, rtip)  

        Omega = SVector(0.0, 0.0, -env.RS(t))


        ### Solve GXBeam for time step
        system, localhistory, converged = time_domain_analysis!(system, assembly, tvec[i-1:i]; prescribed_conditions=prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega, reset_state=false, initialize=false) #TODO: Add gravity. 



        ### Extract GXBeam outputs
        gxhistory[i] = localhistory[end] 




        ### Update aero inputs from structures.
        for j = 1:na
            aeroV[j] = interpolate_velocity(interpolationpoints[j], assembly, localhistory[end])
            def_thetax[j] = interpolate_angle(interpolationpoints[j], assembly, localhistory[end]) 
        end



        if verbose & (mod(i, speakiter)==0)
            println("")
            println("Simulation time: ", t)
            # @show Vxvec[end], Vyvec[end]
            # @show maximum(N[i,:]), maximum(T[i,:])
        end
    end

    Fx = N.*cos(precone) #TODO: The N and T might need some rotating for precone. I'm not sure that I'm doing this right.
    Fy = T.*cos(precone) #I copied Dr. Ning's code to get the thrust and torque, but it seems odd that they are both multiplied by the cosine of precone. Which makes sense. 

    return (N=N, T=T, Fx=Fx, Fy=Fy), cchistory, xds, gxhistory  
end
