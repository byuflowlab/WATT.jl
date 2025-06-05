export environment

abstract type Environment end

struct SimpleEnvironment{TF, F<:Function, G<:Function, H<:Function, J<:Function, K<:Function, L<:Function, M<:Function, N<:Function} <: Environment
    rho::TF #Fluid Density (kg/m^3)
    mu::TF #Fluid dynamic viscosity
    a::TF #Fluid Speed of Sound (m/s)
    shearexp::TF #Shear exponent
    U::F #Freestream Velocity (m/s) (Ux, Uy, Uz) #TODO: Maybe switch to F1, F2, .... That should be a little cleaner. 
    Omega::G #Freestream Swirling Velocity (rads/s) (Omega_x, Omega_y, Omega_z)
    Udot::H #Derivative of the freestream velocity w.r.t. time
    Omegadot::J # Derivative of the rotation rate of the turbine w.r.t. time
    Vinf::K #Magnitude of the freestream velocity (m/s)
    RS::L #Rotation rate of Turbine (rads/s)
    Vinfdot::M
    RSdot::N
end

"""
    environment(rho, mu, a, U, Omega, shearexp)

Steady inflow conditions, velocities prescribed at the rotor plane 
(unvarying in horizontal space), with a shear profile.
Fluid properties are assumed constant throughout the field (incompressible). 

# Arguments
- `rho::Number`: The flow field density
- `mu::Number`: The flow field kinematic viscosity
- `a::Number`: The flow field speed of sound 
- `U::Number`: The freestream velocity
- `Omega::Number`: The rotation rate of the rotor
- `shearexp::Number`: The shear exponent of the wind shear profile

# Output
- SimpleEnvironment::Environment 
"""
function environment(rho::Number, mu::Number, a::Number, U::Number, Omega::Number, shearexp::Number)  #TODO: Originally using types for multiple dispatch. Need to change back. 
    ufun(t) = SVector(U, 0.0, 0.0)
    omegafun(t) = SVector(0.0, 0.0, 0.0)
    udotfun(t) = SVector(0.0, 0.0, 0.0)
    omegadotfun(t) = SVector(0.0, 0.0, 0.0)
    Vinf(t) = U
    RS(t) = Omega
    Vinfdot(t) = 0.0
    RSdot(t) = 0.0
    return SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end

function environment(filename::String, rho::Number, mu::Number, a::Number, Omega::Number, shearexp::Number; fit=Akima)  #TODO: Originally using types for multiple dispatch. Need to change back. 

    turb = readdlm(filename, skipstart=4)
    n, m = size(turb)
    tvec = range(turb[1, 1], turb[n, 1], length=n) #Because the file doesn't get the time vector correctly. 
    Ufit = fit(tvec, turb[:, 2])
    Vfit = fit(tvec, turb[:, 5])
    Wfit = fit(tvec, turb[:, 6])

    ufun(t) = SVector(Ufit(t), Vfit(t), Wfit(t))
    omegafun(t) = SVector(0.0, 0.0, 0.0)
    udotfun(t) = SVector(0.0, 0.0, 0.0)
    omegadotfun(t) = SVector(0.0, 0.0, 0.0)
    Vinf(t) = Ufit(t)
    RS(t) = Omega
    Vinfdot(t) = 0.0
    RSdot(t) = 0.0
    return SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end



"""
    evaluate_flowfield_velocity(env, hubht, x, y, z, t)

A function to retrieve the free-stream velocity vector of a simple environment. Note that this function is part of a family that operates on Environment types. 

*Arguments*
- `env::SimpleEnvironment`: 
- `hubht::Number`: The hub height of the wind turbine. 
- `x`, `y`, or `z::Number`: The position coordinates in the global reference frame of the point of interest.
- `t::Number`: The time at which to evaluate the 
"""
function evaluate_flowfield_velocity(env::SimpleEnvironment, hubht, x, y, z, t)
    factor = (z/hubht)^env.shearexp

    return env.U(t).*factor
end







"""
    get_aero_velocities(rotor, blade, env, t, idx, azimuth)

Calculate the velocities in the airfoil reference frame based on the flow-field, and blade geometry at the given time step. Note that in this function the blade is assumed to be infinitely stiff, so there is no structural deflection or motion. 

*Arguments*
- `rotor::Rotor`: The rotor being analyzed. 
- `blade::Blade`:
- `env::Environment`:
- `t::Number`: The time of that you want evaluated. 
- `idx::Int`: The index of the blade node that you want the velocity of. 
- `azimuth::Number`: The azimuthal position of the blade (radians).
"""
function get_aero_velocities(rotor::Rotor, blade::Blade, env::Environment, t, idx, azimuth)

    ### Unpack
    yaw = -rotor.yaw
    tilt = rotor.tilt
    hubht = rotor.hubht

    #extract the local node coordinates (in the Hub-rotating reference frame)
    rbc_x = blade.rx[idx] #Leadlag
    rbc_y = blade.ry[idx] #Flapwise
    rbc_z = blade.rz[idx] #Radial
    
    sweep = -blade.thetax[idx]  #The sweep is negative in the given reference frame 
    curve = blade.thetay[idx]
    precone = blade.precone 


    ### Get the velocities from the freestream. 
    #Rotate the blade positions from the blade center reference frame to the global. 
    rgx, rgy, rgz = transform_BC_G(rbc_x, rbc_y, rbc_z, azimuth, precone, tilt, yaw)


    # Retrieve the flow field velocities
    Ug = evaluate_flowfield_velocity(env, hubht, rgx, rgy, rgz + hubht, t) 
    # I added hubht to rgz to translate the vector to the top of the tower. 

    #Rotate the velocity into the local frame 
    ulx_wind, uly_wind, _ = transform_G_L(Ug..., azimuth, curve, precone, sweep, tilt, yaw)



    ### Get the rotational velocities
    rhr_x, rhr_y, rhr_z = transform_BC_HR(rbc_x, rbc_y, rbc_z, precone)
    Omega_hr = (env.RS(t), 0, 0) #Angular velocity in hub-rotating reference frame

    #Convert angular velocity to linear velocity 
    vx_rot, vy_rot, vz_rot = cross(Omega_hr, (rhr_x, rhr_y, rhr_z))


    #Convert from the hub-rotating frame to the local airfoil frame
    ulx_rot, uly_rot, _ = transform_HR_L(vx_rot, vy_rot, vz_rot, curve, sweep, precone)

    

    ### Sum the velocities in the airfoil reference frame. 
    Vx = ulx_wind - ulx_rot
    Vy = uly_wind - uly_rot
    # Vz = ulz_wind - ulz_rot
    #Note: The rotational velocities are subtracted rather than added because that's
    #converting from structural velocity to the aerodynamic velocity. 

    return Vx, Vy
end

function get_aerostructural_velocities(env::Environment, aerov, t, r, azimuth, precone, tilt, yaw, hubht)
    #Todo: What is this doing? Is it needed? 
    ### Extract the aero velocities due to the environment, precone, tilt, taw, and shear. 
    vxenv, vyenv = get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

    ### Get the relative structural velocities
    vxrel = aerov[3] 
    vyrel = aerov[2] + env.RS(t)*r  #The rotational velocities are included in the GXBeam output velocity. So to get the relative velocities, I need to subtract out the rotational component. But the rotational component is negative in the structural frame, so subtract a negative is add the rotational velocity back in. 
    

    ### Return the total x and y velocities (x and y in the aerodynamic frame)
    # If the structure is moving in the negative y direction, then the airfoil should see a faster velocity. 
    return vxenv+vxrel, vyenv-vyrel
end


"""
    get_aerostructural_velocities(rotor, blade, env, t, idx, azimuth, dx, dy, dz, vsx, vsy, vsz) -> Vx, Vy

Calculate the velocities in the airfoil reference frame based on the flow-field, blade geometry, and structural deflections and motions at the given time step. 

*Arguments*
- `rotor::Rotor`: The rotor being analyzed. 
- `blade::Blade`:
- `env::Environment`:
- `t::Number`: The time of that you want evaluated. 
- `told::Number`: The time of the previous time step (to remove the rotational velocity from )
- `idx::Int`: The index of the blade node that you want the velocity of. 
- `azimuth::Number`: The azimuthal position of the blade (radians).
- `delta::StaticArray{Number, 3}`: The structural defelctions of the node. 
- `Vs::StaticArray{Number, 3}`: The structural velocities of the node. 
"""
function get_aerostructural_velocities(rotor::Rotor, blade::Blade, env::Environment, t, idx, azimuth, delta, delta_theta, Vs)

    ### Unpack
    yaw = rotor.yaw
    tilt = rotor.tilt
    hubht = rotor.hubht

    #Transform the structural deflections into the aerodynamic reference frame. 
    dx = -delta[3]
    dy = delta[2]
    dz = delta[1]


    rbc_x = blade.rx[idx] + dx #Leadlag
    rbc_y = blade.ry[idx] + dy #Flapwise
    rbc_z = blade.rz[idx] + dz #Radial
    
    sweep = -(blade.thetax[idx] - delta_theta[3])  #The sweep is negative in the given reference frame 
    curve = blade.thetay[idx] + delta_theta[2]
    precone = blade.precone 


    ### Get the velocities from the freestream. 
    #Rotate the blade positions from the hub reference frame to the global. 
    rgx, rgy, rgz = transform_BC_G(rbc_x, rbc_y, rbc_z, azimuth, precone, tilt, yaw)

    # Retrieve the flow field velocities
    Ug = evaluate_flowfield_velocity(env, hubht, rgx, rgy, rgz + hubht, t) 
    # I added hubht to rgz to translate the vector to the top of the tower.

    #Rotate the velocity into the local frame 
    # ulx_wind, uly_wind, _ = transform_G_L(Ug..., azimuth, curve, precone, sweep, tilt, yaw)
    ulx_wind, uly_wind, _ = transform_G_L(Ug..., azimuth, curve, precone-delta_theta[2], sweep, tilt, yaw) #Note: Didn't have a large effect on the solution. 

    


    ### Get the rotational velocities
    rhr_x, rhr_y, rhr_z = transform_BC_HR(rbc_x, rbc_y, rbc_z, precone) #TODO: I'm not sure that I need to rotate by the precone... because the deflections should move the aero nodes to there given position in the hub-rotating reference frame. 
    Omega_hr = (env.RS(t), 0, 0) #Angular velocity in hub-rotating reference frame

    #Convert angular velocity to linear velocity 
    vx_rot, vy_rot, vz_rot = cross(Omega_hr, (rhr_x, rhr_y, rhr_z))


    #Convert from the hub-rotating frame to the local airfoil frame
    # ulx_rot, uly_rot, _ = transform_HR_L(vx_rot, vy_rot, vz_rot, curve, sweep, precone)
    ulx_rot, uly_rot, _ = transform_HR_L(vx_rot, vy_rot, vz_rot, curve, sweep, precone-delta_theta[2]) #Note: Didn't have a large effect on the solution. 



    ### Transform the structural velocities from the hub-rotating into the local frame
    # usx, usy, _ = transform_HR_L(Vs..., curve, sweep, precone)
    usx, usy, _ = transform_HR_L(Vs..., curve, sweep, precone-delta_theta[2]) #Note: Didn't have a large effect on the solution. 


    ### Sum the velocities in the airfoil reference frame. 
    Vx = ulx_wind - ulx_rot + usx
    Vy = uly_wind - uly_rot + usy
    # Vz = ulz_wind - ulz_rot + usz
    #Note: The rotational velocities are subtracted rather than added because that's
    #converting from structural velocity to the aerodynamic velocity. 

    # return Vx, Vy

    ### Trying to apply change in twist
    st, ct = sincos(delta_theta[1])

    Vxnew = Vx*ct - Vy*st
    Vynew = -Vx*st + Vy*ct

    return Vxnew, Vynew
end