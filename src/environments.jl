struct Environment{TF, Tf}
    rho::TF #Fluid Density (kg/m^3)
    mu::TF #Fluid dynamic viscosity
    a::TF #Fluid Speed of Sound (m/s)
    Uinf::Tf #Freestream Velocity (m/s)
    Omega::Tf #Rotation rate of Turbine (rads/s)
end