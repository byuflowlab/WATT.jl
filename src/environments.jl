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

    @show azimuth

    # coordinate in azimuthal coordinate system
    x_az = -r*sin(precone)
    z_az = r*cos(precone)
    y_az = 0.0  # could omit (the more general case allows for presweep so this is nonzero) #Todo: If I include those velocities, I'd need to augment this. 
    @show x_az, y_az, z_az

    # get section heights in wind-aligned coordinate system
    heightFromHub = (y_az*sin_azimuth + z_az*cos_azimuth)*cos_tilt - x_az*sin_tilt

    @show heightFromHub

    ### velocity with shear
    factor = (1 + heightFromHub/hubht)^env.shearexp
    V = env.Vinf(t)*factor
    @show V, factor

    # transform wind to blade c.s.
    Vwind_x = V * ((cos_yaw*sin_tilt*cos_azimuth + sine_yaw*sin_azimuth)*sin_precone + cos_yaw*cos_tilt*cos_precone)
    Vwind_y = V * (cos_yaw*sin_tilt*sin_azimuth - sine_yaw*cos_azimuth)

    @show Vwind_x, Vwind_y

    # wind from rotation to blade c.s.
    Vrot_x = -env.RS(t)*y_az*sin_precone
    Vrot_y = env.RS(t)*z_az

    @show Vrot_x, Vrot_y

    # total velocity
    Vx = Vwind_x + Vrot_x
    Vy = Vwind_y + Vrot_y

    return Vx, Vy
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
    # @show y/hubht
    @show factor

    # dely = sqrt(x^2 + y^2 + z^2)
    # factor = (1 + dely/hubht)^env.shearexp
    # @show factor

    # dely = abs(y-hubht)
    # factor = (1 + dely/hubht)^env.shearexp
    # @show factor

    # dely = y-hubht
    # factor = (1 + dely/hubht)^env.shearexp
    # @show factor

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

    rax = blade.rx[idx] #Leadlag
    ray = blade.ry[idx] #Flapwise
    raz = blade.rz[idx] #Radial
    
    # @show azimuth

    # omega = SVector(env.RS(t), 0, 0)
    omega = env.RS(t)
    sweep = -atan(ray, raz) #TODO: Potentially store these to avoid calculating them every iteration. 
    precone = atan(rax, raz)
    # theta_z = 0.0 #pi/4

    @show precone, sweep, yaw, tilt
    @show azimuth, azimuth*180/pi
    # @show blade_azimuth
    # @show tilt, yaw

    R_azimuth = rotate_x(azimuth)
    R_sweep = rotate_x(sweep)
    R_tilt = rotate_y(tilt)
    R_yaw = rotate_z(yaw)
    R_cone = rotate_y(precone)

    ### Get the velocities from the freestream. 
    #Rotate the blade positions from the hub reference frame to the global. 
    #Note: Precone and sweep should be inherent in blade rx, ry, and rz

    #= Tests passed
    A1, B1, C1 = 17/25
    A4, B7, C1 = 17/25
    A6, B7, C1 = 18/25
    A6, B8, C1 = 17/25
    A6, B1, C1 = 17/25
    A5, B8, C1 = 16/25
    A6, B5, C1 = 16/25
    A6, B9, C1 = 23/25
    A7, B9, C1 = 23/25
    =#


    ###### Rotation A
    # R_rothub_g = R_yaw*R_tilt*R_azimuth #1
    # R_rothub_g = R_yaw*R_azimuth*R_tilt #2
    # R_rothub_g = R_tilt*R_azimuth*R_yaw #3
    # R_rothub_g = R_tilt*R_yaw*R_azimuth #4
    # R_rothub_g = R_azimuth*R_yaw*R_tilt #5
    # R_rothub_g = R_azimuth*R_tilt*R_yaw #6
    R_rothub_g = R_yaw'*R_tilt'*R_azimuth' #7
    


    rrh = [rax, ray, raz]

    rg = R_rothub_g*rrh

    rgx, rgy, rgz = rg

    # @show rax, ray, raz
    # rgx, rgy, rgz = rotate_x(rax, ray, raz, azimuth; T=false)
    # @show "azimuthed: ", rgx, rgy, rgz
    # rgx, rgy, rgz = rotate_y(rgx, rgy, rgz, tilt; T=false)
    # @show "tilted: ", rgx, rgy, rgz
    # rgx, rgy, rgz = rotate_z(rgx, rgy, rgz, yaw; T=false)
    # @show "yawed: ", rgx, rgy, rgz

    # @show hubht
    @show rax, ray, raz
    @show rgx, rgy, rgz
    @show rgz+hubht



    # Retrieve the flow field velocities
    Ug = evaluate_flowfield_velocity(env, hubht, rgx, rgy, rgz + hubht, t) 
    # I added hubht to rgz to translate the vector to the top of the tower. 

    @show Ug 



    ###Rotate the velocity into the local frame 
    ######### Rotation B
    # R_g_l = R_yaw*R_tilt*R_cone*R_sweep*R_azimuth #1
    # R_g_l = R_yaw*R_tilt*R_cone*R_azimuth*R_sweep #2
    # R_g_l = R_yaw*R_tilt*R_sweep*R_cone*R_azimuth #3
    # R_g_l = R_yaw*R_tilt*R_sweep*R_azimuth*R_cone #4
    # R_g_l = R_yaw*R_tilt*R_azimuth*R_cone*R_sweep #5
    # R_g_l = R_yaw*R_tilt*R_azimuth*R_sweep*R_cone #6
    # R_g_l = R_tilt*R_yaw*R_cone*R_sweep*R_azimuth #7
    # R_g_l = R_tilt*R_yaw*R_azimuth*R_cone*R_sweep #8
    R_g_l = R_cone'*R_sweep'*R_azimuth*R_tilt*R_yaw #9

    # R_g_l = (R_tilt*R_yaw*R_azimuth)' 

    ul_wind = R_g_l*Ug

    ulx_wind, uly_wind, ulz_wind = ul_wind

    # uax_w, uay_w, uaz_w = rotate_x(Ug..., azimuth; T=true)
    # println("Azimuthed flow: ", uax_w, ", ",  uay_w, ", ", uaz_w)

    # uax_w, uay_w, uaz_w = rotate_x(uax_w, uay_w, uaz_w, sweep; T=true)
    # println("Swept flow: ", uax_w, ", ",  uay_w, ", ", uaz_w)

    # uax_w, uay_w, uaz_w = rotate_y(uax_w, uay_w, uaz_w, precone; T=true)
    # println("Coned flow: ", uax_w, ", ",  uay_w, ", ", uaz_w)

    # uax_w, uay_w, uaz_w = rotate_y(uax_w, uay_w, uaz_w, tilt; T=true)
    # println("Tilted flow: ", uax_w, ", ",  uay_w, ", ", uaz_w)

    # ulx_wind, uly_wind, ulz_wind = rotate_z(uax_w, uay_w, uaz_w, yaw; T=true) # (u local x wind)



    @show ulx_wind, uly_wind, ulz_wind



    ##### -------- Rotating Omega
    Omega_a = (env.RS(t), 0, 0)

    vx_rot, vy_rot, vz_rot = cross(Omega_a, (rax, ray, raz))
    @show vx_rot, vy_rot, vz_rot

    ########## Rotation C
    R_rh_l = R_cone*R_sweep                 #1
    # R_rh_l = R_cone*R_sweep*R_azimuth     #2
    # R_rh_l = R_cone*R_azimuth*R_sweep     #3
    # R_rh_l = R_sweep*R_cone               #4
    # R_rh_l = R_sweep*R_azimuth*R_cone     #5
    # R_rh_l = R_sweep*R_cone*R_azimuth     #6
    # R_rh_l = R_azimuth*R_cone*R_sweep     #7
    # R_rh_l = R_azimuth*R_sweep*R_cone     #8

    # R_rh_l = I

    v_rot = [vx_rot, vy_rot, vz_rot]
    ul_rot = R_rh_l*v_rot

    ulx_rot, uly_rot, ulz_rot = ul_rot

    # vgx, vgy, vgz = rotate_x(vx_rot, vy_rot, vz_rot, azimuth; T=true)
    # println("Azimuthed flow: ", vgx, ", ",  vgy, ", ", vgz)

    # vgx, vgy, vgz = rotate_x(vx_rot, vy_rot, vz_rot, sweep; T=false)
    # println("Swept flow: ", vgx, ", ",  vgy, ", ", vgz)

    # ulx_rot, uly_rot, ulz_rot = rotate_y(vgx, vgy, vgz, precone; T=false)
    # println("Coned flow: ", ulx_rot, ", ",  uly_rot, ", ", ulz_rot)




    @show ulx_rot, uly_rot, ulz_rot
    

    ### Sum and rotate into the airfoil reference frame. 
    Vx = ulx_wind - ulx_rot
    Vy = uly_wind - uly_rot
    Vz = ulz_wind - ulz_rot
    #Note: The rotational velocities are subtracted rather than added because that's
    #converting from structural velocity to the aerodynamic velocity. 

    # V = [Vx, Vy, Vz]
    # @show V
    # R_h_l = R_cone*R_sweep
    # Vl = R_h_l*V
    # Vx, Vy, Vz = Vl

    @show Vx, Vy, Vz

    
    
    return Vx, Vy
end

function get_aerostructural_velocities(env::Environment, aeroV, t, r, azimuth, precone, tilt, yaw, hubht)

    ### Extract the aero velocities due to the environment, precone, tilt, taw, and shear. 
    vxenv, vyenv = get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

    ### Get the relative structural velocities
    vxrel = aeroV[3] 
    vyrel = aeroV[2] + env.RS(t)*r  #Todo. Why do I add the rotational velocity into the structural? -> Oh, that's right, the rotational velocities are included in the GXBeam output velocity. So to get the relative velocities, I need to subtract out the rotational component. But the rotational component is negative in the structural frame, so subtract a negative is add the rotational velocity back in. 
    

    ### Return the total x and y velocities (x and y in the aerodynamic frame)
    # If the structure is moving in the negative y direction, then the airfoil should see a faster velocity. 
    return vxenv+vxrel, vyenv-vyrel
end


"""
    get_aerostructural_velocities(rotor, blade, env, t, idx, azimuth, dx, dy, dz, vsx, vsy, vsz)
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
function get_aerostructural_velocities(rotor::Rotor, blade::Blade, env::Environment, t, idx, azimuth, delta, Vs)

    ### Unpack
    yaw = rotor.yaw
    tilt = rotor.tilt
    hubht = rotor.hubht

    #Transform the structural deflections into the aerodynamic reference frame. 
    dx = -delta[3]
    dy = delta[2]
    dz = delta[1]

    rax = blade.rx[idx] + dx #Leadlag
    ray = blade.ry[idx] + dy #Flapwise
    raz = blade.rz[idx] + dz #Radial
    # @show rax, ray, raz

    omega = env.RS(t)
    azimuth -= pi/2

    # @show azimuth

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
    uax_wind, uay_wind, uaz_wind = rotate_vector(Ug[1], Ug[2], Ug[3], -azimuth, -yaw, -tilt; forward=false) #Todo: I'm a couple of millimeters of on the z velocity which is propogating down (When I do multiple rotations at once.)
    # uax_wind, uay_wind, uaz_wind = rotate_vector(Ug[1], Ug[2], Ug[3], azimuth, 0.0, 0.0; forward=false)
    # uax_wind, uay_wind, uaz_wind = rotate_vector(uax_wind, uay_wind, uay_wind, 0.0, 0.0, tilt; forward=false)
    # uax_wind, uay_wind, uaz_wind = rotate_vector(uax_wind, uay_wind, uay_wind, 0.0, yaw, 0.0; forward=false)

    # @show uax_wind, uay_wind, uaz_wind


    ### Get the velocities from the rotation.
    uax_rot = 0.0
    uay_rot = omega*raz #Todo: Changing the sign of this doesn't seem to affect the outcome... :|
    uaz_rot = omega*ray

    # @show uay_rot, uaz_rot

    # ### Transform the structural velocities into the aero frame
    # usx = Vs[3]
    # usy = -Vs[2] 
    # usz = -Vs[1]

    ### Sum and rotate into the airfoil reference frame. 
    uax = uax_wind + uax_rot + Vs[1]
    uay = uay_wind + uay_rot + Vs[2]
    uaz = uaz_wind + uaz_rot + Vs[3]

    theta_x = atan(ray, raz) #TODO: Potentially store these to avoid calculating them every iteration. 
    theta_y = atan(rax, raz)
    theta_z = 0.0

    # @show uax, uay, uaz
    # @show theta_x, theta_y, theta_z

    Vx, Vy, _ = rotate_vector(uax, uay, uaz, theta_x, theta_y, theta_z; forward=true) #Normal velocity, tangential velocity, radial velocity
    #Todo: Changing forward to true or to false doesn't have any affect on the outcome!!! -> That's because the section that I've been testing it on doesn't have precone... but I don't know how it passed the precone tests.... but... 
    
    return Vx, Vy
end