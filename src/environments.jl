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


function get_aero_velocities(env::Environment, t, r, azimuth, precone, tilt, yaw, hubht)
    sine_yaw = sin(yaw)
    cos_yaw = cos(yaw)
    sin_tilt = sin(tilt)
    cos_tilt = cos(tilt)
    sin_azimuth = sin(azimuth)
    cos_azimuth = cos(azimuth)
    sin_precone = sin(precone)
    cos_precone = cos(precone)

    # coordinate in azimuthal coordinate system
    x_az = -r*sin(precone)
    z_az = r*cos(precone)
    y_az = 0.0  # could omit (the more general case allows for presweep so this is nonzero) #Todo: If I include those velocities, I'd need to augment this. 
    # @show x_az, z_az, y_az

    # get section heights in wind-aligned coordinate system
    heightFromHub = (y_az*sin_azimuth + z_az*cos_azimuth)*cos_tilt - x_az*sin_tilt

    ### velocity with shear
    V = env.Vinf(t)*(1 + heightFromHub/hubht)^env.shearexp
    # @show V

    # transform wind to blade c.s.
    Vwind_x = V * ((cos_yaw*sin_tilt*cos_azimuth + sine_yaw*sin_azimuth)*sin_precone + cos_yaw*cos_tilt*cos_precone)
    Vwind_y = V * (cos_yaw*sin_tilt*sin_azimuth - sine_yaw*cos_azimuth)

    @show Vwind_x, Vwind_y

    # wind from rotation to blade c.s.
    Vrot_x = -env.RS(t)*y_az*sin_precone
    Vrot_y = env.RS(t)*z_az

    # @show Vrot_x, Vrot_y

    # total velocity
    Vx = Vwind_x + Vrot_x
    Vy = Vwind_y + Vrot_y

    return Vx, Vy
end

function evaluate_flowfield_velocity(env::SimpleEnvironment, hubht, x, y, z, t)
    factor = (y/hubht)^env.shearexp
    # factor = (1 + y/hubht)^env.shearexp
    return env.U(t).*factor
end


function get_aero_velocities(rotor::Rotor, blade::Blade, env::Environment, t, idx, azimuth, dx, dy, dz)

    ### Unpack
    yaw = -rotor.yaw
    tilt = rotor.tilt
    hubht = rotor.hubht

    rax = blade.rx[idx] + dx
    ray = blade.ry[idx] + dy
    raz = blade.rz[idx] + dz
    # @show rax, ray, raz

    omega = env.RS(t)

    ### Get the velocities from the freestream. 
    # Rotate the positions into the global reference frame
    rgx, rgy, rgz = rotate_vector(rax, ray, raz, azimuth, yaw, tilt) #Translate the blade up to the location where it is attached to the top of the tower. 
    # rgx, rgy, rgz = rotate_vector(rax, ray, raz, azimuth, 0.0, 0.0) 
    # rgx, rgy, rgz = rotate_vector(rgx, rgy, rgz, 0.0, yaw, tilt) 
    # rgx, rgy, rgz = rotate_vector(rgx, rgy, rgz, azimuth, 0.0, 0.0) 
    rgy += hubht 

    # Retrieve the flow field velocities
    Ug = evaluate_flowfield_velocity(env, hubht, rgx, rgy, rgz, t)

    # @show Ug 

    # Rotate the freestream velocities into the aerodynamic frame (where the YZ plane aligns with the rotor plane)
    uax_wind, uay_wind, uaz_wind = rotate_vector(Ug[1], Ug[2], Ug[3], azimuth, yaw, tilt; forward=false) #Todo: I'm a couple of millimeters of on the z velocity which is propogating down (When I do multiple rotations at once.)
    # uax_wind, uay_wind, uaz_wind = rotate_vector(Ug[1], Ug[2], Ug[3], azimuth, 0.0, 0.0; forward=false)
    # uax_wind, uay_wind, uaz_wind = rotate_vector(uax_wind, uay_wind, uay_wind, 0.0, 0.0, tilt; forward=false)
    # uax_wind, uay_wind, uaz_wind = rotate_vector(uax_wind, uay_wind, uay_wind, 0.0, yaw, 0.0; forward=false)

    # @show uax_wind, uay_wind, uaz_wind


    ### Get the velocities from the rotation.
    # uax_rot = 0.0
    uay_rot = -omega*raz #Todo: Changing the sign of this doesn't seem to affect the outcome... :|
    uaz_rot = omega*ray

    # @show uay_rot, uaz_rot

    ### Sum and rotate into the airfoil reference frame. 
    uax = uax_wind
    uay = uay_wind + uay_rot
    uaz = uaz_wind + uaz_rot

    theta_x = atan(raz, ray) #TODO: Potentially store these to avoid calculating them every iteration. 
    theta_y = 0.0
    theta_z = atan(rax, ray)

    # @show uax, uay, uaz
    # @show theta_x, theta_y, theta_z

    Vx, _, Vy = rotate_vector(uax, uay, uaz, theta_x, theta_y, theta_z; forward=true) #Normal velocity, radial velocity, tangential velocity
    #Todo: Changing forward to true or to false doesn't have any affect on the outcome!!! -> That's because the section that I've been testing it on doesn't have precone... but I don't know how it passed the precone tests.... but... 
    
    return Vx, Vy
end

