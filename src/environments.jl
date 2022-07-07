struct Environment{TF, F<:Function, G<:Function, H<:Function, J<:Function, K<:Function, L<:Function, M<:Function, N<:Function} 
    rho::TF #Fluid Density (kg/m^3)
    mu::TF #Fluid dynamic viscosity
    a::TF #Fluid Speed of Sound (m/s)
    U::F #Freestream Velocity (m/s) (Ux, Uy, Uz) #TODO: Maybe switch to F1, F2, .... That should be a little cleaner. 
    Omega::G #Freestream Swirling Velocity (rads/s) (Omega_x, Omega_y, Omega_z)
    Udot::H #Derivative of the freestream velocity w.r.t. time
    Omegadot::J # Derivative of the rotation rate of the turbine w.r.t. time
    Vinf::K #Magnitude of the freestream velocity (m/s)
    RS::L #Rotation rate of Turbine (rads/s)
    Vinfdot::M
    RSdot::N
end


function environment(rho, mu, a, U::Float64, Omega::Float64)  
    ufun(t) = SVector(U, 0.0, 0.0)
    omegafun(t) = SVector(0.0, 0.0, 0.0)
    udotfun(t) = SVector(0.0, 0.0, 0.0)
    omegadotfun(t) = SVector(0.0, 0.0, 0.0)
    Vinf(t) = U
    RS(t) = Omega
    Vinfdot(t) = 0.0
    RSdot(t) = 0.0
    return Environment(rho, mu, a, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
end


function get_structural_velocity(env::Environment, t)
    Ua = env.U(t)
    return SVector(Ua[3], Ua[2], -Ua[1])
end