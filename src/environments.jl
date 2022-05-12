struct Environment{TF, F<:Function, G<:Function, H<:Function, J<:Function}
    rho::TF #Fluid Density (kg/m^3)
    mu::TF #Fluid dynamic viscosity
    a::TF #Fluid Speed of Sound (m/s)
    U::F #Freestream Velocity (m/s) #Todo: How to type this when I expect Fs to be passed in here? I threw in a temporary fix. I'm not sure what the correct fix is. 
    # vinf #Updraft velocity (m/s) 
    Omega::G #Rotation rate of Turbine (rads/s) #TODO: Will this need to come out? Is Omega design variable? -> If so, it might be worthwhile to actually turn this into a model instead of an exterior structure... Although, if I recall, we usually run the turbine for a set of TSRs which is essentially the same as a set of Omega. 
    Udot::H #Derivative of the freestream velocity w.r.t. time
    Omegadot::J # Derivative of the rotation rate of the turbine w.r.t. time
end


#then Todo: update Riso to use function velocities. 

function environment(rho, mu, a, U::Float64, Omega::Float64, Udot::Float64, Omegadot::Float64) #Todo: I'm not sure that declaring these Floats is the best way to change from single values (like a float) to a function. -> Dr. Ning says typing like this isn't neccesary, only when using multiple dispatch. I'm not using multiple dispatch right now because there is a difference in function spelling. 
    ufun(t) = U
    omegafun(t) = Omega
    udotfun(t) = 0.0
    omegadotfun(t) = 0.0
    return Environment(rho, mu, a, ufun, omegafun, udotfun, omegadotfun)
end