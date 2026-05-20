abstract type Environment end


"""
    Constant{T}

Callable struct wrapping a scalar or `SVector` value. Calling it with any
argument returns the stored value.

**Fields**
- `val::T`: The value to return.
"""
struct Constant{T}
    val::T
end
(c::Constant)(_) = c.val


"""
    TurbulentInflow{A}

Callable struct that evaluates three interpolants (one per velocity component)
at a given time and returns the result as an `SVector{3}`. Replaces the
anonymous closures that previously held an Akima fit per component.

**Fields**
- `Ufit::A`, `Vfit::A`, `Wfit::A`: Per-component interpolants (typically `FLOWMath.Akima`).
"""
struct TurbulentInflow{A}
    Ufit::A
    Vfit::A
    Wfit::A
end
(t::TurbulentInflow)(s) = SVector(t.Ufit(s), t.Vfit(s), t.Wfit(s))


"""
    ScalarFit{A}

Callable struct wrapping a 1D interpolant.

**Fields**
- `fit::A`: The interpolant (typically `FLOWMath.Akima`).
"""
struct ScalarFit{A}
    fit::A
end
(f::ScalarFit)(t) = f.fit(t)


"""
    SimpleEnvironment{TF, F, G, H, J, K, L, M, N}

Inflow conditions and fluid properties at the rotor plane. Each of the eight
velocity/rotation fields is a callable struct of `t` — calling
`env.Vinf(t)` returns the freestream magnitude at time `t`, etc. Using
named callable structs (rather than anonymous closures) keeps the
type-parameter count stable across distinct environments and lets JLD2
roundtrip an env across sessions.

**Fields**
- `rho::TF`: Fluid density (kg/m³).
- `mu::TF`: Fluid dynamic viscosity.
- `a::TF`: Speed of sound (m/s).
- `shearexp::TF`: Power-law shear exponent.
- `U::F`: `t -> SVector{3}` freestream velocity vector.
- `Omega::G`: `t -> SVector{3}` freestream rotation rate.
- `Udot::H`, `Omegadot::J`: Their time derivatives.
- `Vinf::K`: `t -> TF` magnitude of freestream velocity.
- `RS::L`: `t -> TF` rotor angular speed (rad/s).
- `Vinfdot::M`, `RSdot::N`: Time derivatives.
"""
struct SimpleEnvironment{TF, F, G, H, J, K, L, M, N} <: Environment
    rho::TF
    mu::TF
    a::TF
    shearexp::TF
    U::F
    Omega::G
    Udot::H
    Omegadot::J
    Vinf::K
    RS::L
    Vinfdot::M
    RSdot::N
end

"""
    environment(rho, mu, a, U, Omega, shearexp) -> SimpleEnvironment

Steady inflow conditions with a power-law shear profile. Fluid properties are
treated as constant in space (incompressible).

**Arguments**
- `rho::Number`: Density.
- `mu::Number`: Dynamic viscosity.
- `a::Number`: Speed of sound.
- `U::Number`: Freestream velocity magnitude.
- `Omega::Number`: Rotor angular speed (rad/s).
- `shearexp::Number`: Power-law shear exponent.
"""
function environment(rho::Number, mu::Number, a::Number, U::Number, Omega::Number, shearexp::Number)
    ufun       = Constant(SVector(U, 0.0, 0.0))
    omegafun   = Constant(SVector(0.0, 0.0, 0.0))
    udotfun    = Constant(SVector(0.0, 0.0, 0.0))
    omegadotfun = Constant(SVector(0.0, 0.0, 0.0))
    Vinf       = Constant(U)
    RS         = Constant(Omega)
    Vinfdot    = Constant(0.0)
    RSdot      = Constant(0.0)
    return SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end

"""
    environment(filename, rho, mu, a, Omega, shearexp; fit=Akima) -> SimpleEnvironment

Construct a `SimpleEnvironment` from a turbulent inflow file. The full time
range of the file is used.
"""
function environment(filename::String, rho::Number, mu::Number, a::Number, Omega::Number, shearexp::Number; fit=Akima)
    turb = readdlm(filename, skipstart=4)
    n, _ = size(turb)
    tvec = range(turb[1, 1], turb[n, 1], length=n)
    Ufit = fit(tvec, turb[:, 2])
    Vfit = fit(tvec, turb[:, 5])
    Wfit = fit(tvec, turb[:, 6])

    ufun       = TurbulentInflow(Ufit, Vfit, Wfit)
    omegafun   = Constant(SVector(0.0, 0.0, 0.0))
    udotfun    = Constant(SVector(0.0, 0.0, 0.0))
    omegadotfun = Constant(SVector(0.0, 0.0, 0.0))
    Vinf       = ScalarFit(Ufit)
    RS         = Constant(Omega)
    Vinfdot    = Constant(0.0)
    RSdot      = Constant(0.0)
    return SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end

"""
    environment(filename, time, rho, mu, a, Omega, shearexp; fit=Akima) -> SimpleEnvironment

Construct a `SimpleEnvironment` from a turbulent inflow file. If the file is
shorter than `time`, it is tiled forward-then-reversed to fill `time` without
introducing a discontinuity at the seam.
"""
function environment(filename::String, time::Number, rho::Number, mu::Number, a::Number, Omega::Number, shearexp::Number; fit=Akima)
    turb = readdlm(filename, skipstart=4)
    n, _ = size(turb)
    if turb[n, 1] >= time
        nn = findfirst(x -> x >= time, turb[:, 1])
        tvec = range(turb[1, 1], turb[nn, 1], length=nn)
        Ufit = fit(tvec, turb[1:nn, 2])
        Vfit = fit(tvec, turb[1:nn, 5])
        Wfit = fit(tvec, turb[1:nn, 6])
    else
        tmax = turb[n, 1]
        nrepeat = ceil(Int, time/tmax)
        nn = nrepeat*n
        tvec = range(turb[1, 1], turb[n, 1]*nrepeat, length=nn)

        Uvec = vcat([isodd(i) ? turb[:, 2] : reverse(turb[:, 2]) for i in 1:nrepeat]...)
        Vvec = vcat([isodd(i) ? turb[:, 5] : reverse(turb[:, 5]) for i in 1:nrepeat]...)
        Wvec = vcat([isodd(i) ? turb[:, 6] : reverse(turb[:, 6]) for i in 1:nrepeat]...)
        Ufit = fit(tvec, Uvec)
        Vfit = fit(tvec, Vvec)
        Wfit = fit(tvec, Wvec)
    end

    ufun       = TurbulentInflow(Ufit, Vfit, Wfit)
    omegafun   = Constant(SVector(0.0, 0.0, 0.0))
    udotfun    = Constant(SVector(0.0, 0.0, 0.0))
    omegadotfun = Constant(SVector(0.0, 0.0, 0.0))
    Vinf       = ScalarFit(Ufit)
    RS         = Constant(Omega)
    Vinfdot    = Constant(0.0)
    RSdot      = Constant(0.0)
    return SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end



"""
    evaluate_flowfield_velocity(env, hubht, x, y, z, t) -> SVector{3}

Free-stream velocity vector at position `(x, y, z)` and time `t` for a
`SimpleEnvironment` with a power-law shear profile.

**Arguments**
- `env::SimpleEnvironment`: Environment object.
- `hubht::Number`: Hub height.
- `x`, `y`, `z::Number`: Global-frame coordinates of the query point.
- `t::Number`: Time.
"""
function evaluate_flowfield_velocity(env::SimpleEnvironment, hubht, x, y, z, t)
    factor = (z/hubht)^env.shearexp

    return env.U(t).*factor
end




"""
    get_aero_velocities(rotor, blade, env, t, idx, azimuth) -> (Vx, Vy)

Velocities in the airfoil reference frame given the flow field and rigid blade
geometry at a single time step. No structural deflection.

**Arguments**
- `rotor::Rotor`
- `blade::Blade`
- `env::Environment`
- `t::Number`: Time of evaluation.
- `idx::Int`: Blade-node index.
- `azimuth::Number`: Azimuth (rad).
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
    #Note: The rotational velocities are subtracted rather than added because that's
    #converting from structural velocity to the aerodynamic velocity.

    return Vx, Vy
end

"""
    get_aerostructural_velocities(rotor, blade, env, t, idx, azimuth, delta, delta_theta, Vs) -> (Vx, Vy)

Velocities in the airfoil reference frame given the flow field, rigid blade
geometry, **and** structural deflection/motion at a single time step.

**Arguments**
- `rotor::Rotor`
- `blade::Blade`
- `env::Environment`
- `t::Number`: Time of evaluation.
- `idx::Int`: Blade-node index.
- `azimuth::Number`: Azimuth (rad).
- `delta::SVector{3}`: Linear structural deflection.
- `delta_theta::SVector{3}`: Angular structural deflection.
- `Vs::SVector{3}`: Structural velocity at this node.
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

    sweep = -(blade.thetax[idx] - delta_theta[3])
    curve = blade.thetay[idx] + delta_theta[2]
    precone = blade.precone


    ### Get the velocities from the freestream.
    rgx, rgy, rgz = transform_BC_G(rbc_x, rbc_y, rbc_z, azimuth, precone, tilt, yaw)
    Ug = evaluate_flowfield_velocity(env, hubht, rgx, rgy, rgz + hubht, t)
    ulx_wind, uly_wind, _ = transform_G_L(Ug..., azimuth, curve, precone-delta_theta[2], sweep, tilt, yaw)


    ### Get the rotational velocities
    rhr_x, rhr_y, rhr_z = transform_BC_HR(rbc_x, rbc_y, rbc_z, precone)
    Omega_hr = (env.RS(t), 0, 0)
    vx_rot, vy_rot, vz_rot = cross(Omega_hr, (rhr_x, rhr_y, rhr_z))
    ulx_rot, uly_rot, _ = transform_HR_L(vx_rot, vy_rot, vz_rot, curve, sweep, precone-delta_theta[2])


    ### Transform the structural velocities from the hub-rotating into the local frame
    usx, usy, _ = transform_HR_L(Vs..., curve, sweep, precone-delta_theta[2])


    ### Sum the velocities in the airfoil reference frame.
    Vx = ulx_wind - ulx_rot + usx
    Vy = uly_wind - uly_rot + usy

    ### Apply change in twist
    st, ct = sincos(delta_theta[1])
    Vxnew = Vx*ct - Vy*st
    Vynew = -Vx*st + Vy*ct

    return Vxnew, Vynew
end
