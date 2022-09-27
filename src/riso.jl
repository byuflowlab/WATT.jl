#=
Code to interact with the DynamicStallModels Package, specifically for the Risø model. 

=#



function initializeDSmodel(dsmodel::DS.Riso, dsmodelinit::Hansen, solver::RK4, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

    if turbine
        thetavec = [-((twistvec[i] + pitch) - phivec[i]) for i = 1:na] 
    else
        thetavec = [((twistvec[i] + pitch) - phivec[i]) for i = 1:na]
    end

    pds = vcat(Vxvec, Vxdotvec, thetavec, zeros(na), chordvec)

    # ode = createrisoode(blade)

    xds = zeros(nt, 4*na) #Note: I thought about initializing this as an undefined array, but... the first vector needs to be mostly zeros, and I feel like this will be the easiest. 
    xds[1,4:4:4*na] .= 1

    return dsmodel, xds, pds
end

function initializeDSmodel(dsmodel::DS.Riso, dsmodelinit::Steady, solver::RK4, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

    if turbine
        thetavec = [-((twistvec[i] + pitch) - phivec[i]) for i = 1:na] 
    else
        thetavec = [((twistvec[i] + pitch) - phivec[i]) for i = 1:na]
    end

    pds = vcat(Vxvec, Vxdotvec, thetavec, zeros(na), chordvec)
    
    # ode = createrisoode(blade)

    xds = zeros(nt, 4*na) 
    xds[1,4:4:4*na] .= 1

    prob = DE.SteadyStateProblem{false}(dsmodel, xds[1,:], pds) #Todo: I don't know if DifferntialEquations will be okay with this. 
    sol = DE.solve(prob)

    xds[1,:] = sol.u

    return dsmodel, xds, pds
end

function initializeDSmodel(dsmodel::DS.Riso, dsmodelinit::Hansen, solver::DiffEQ, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

    if turbine
        thetavec = [-((twistvec[i] + pitch) - phivec[i]) for i = 1:na] 
    else
        thetavec = [((twistvec[i] + pitch) - phivec[i]) for i = 1:na]
    end

    pds = vcat(Vxvec, Vxdotvec, thetavec, zeros(na), chordvec)


    # ode = createrisoode(blade)

    xds = zeros(nt, 4*na)
    xds[1,4:4:4*na] .= 1

    prob = DE.ODEProblem{false}(dsmodel, xds[1,:], (tvec[1], tvec[end]), pds)
    integrator = DE.init(prob, solver.algorithm)

    return integrator, xds, pds
end

function initializeDSmodel(dsmodel::DS.Riso, dsmodelinit::Steady, solver::DiffEQ, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

    if turbine
        thetavec = [-((twistvec[i] + pitch) - phivec[i]) for i = 1:na] 
    else
        thetavec = [((twistvec[i] + pitch) - phivec[i]) for i = 1:na]
    end

    pds = vcat(Vxvec, Vxdotvec, thetavec, zeros(na), chordvec)

    # ode = createrisoode(blade)

    xds = zeros(nt, 4*na) 
    xds[1,4:4:4*na] .= 1

    prob = DE.SteadyStateProblem{false}(dsmodel, xds[1,:], pds)
    sol = DE.solve(prob)

    xds[1,:] = sol.u

    prob = DE.ODEProblem{false}(dsmodel, xds[1,:], (tvec[1], tvec[end]), pds)
    integrator = DE.init(prob, solver.algorithm)

    return integrator, xds, pds
end

function update_aero_parameters!(dsmodel::DS.Riso, turbine::Bool, pds, na, rvec, Wvec, phivec, twistvec, pitch, env, t)
    # pds = vcat(Vxvec, Vxdotvec, v, vdot, thetavec, thetadot, chordvec)
        for j = 1:na
            pds[j] = Wvec[j]
            idx = na + j
            pds[idx] = sqrt(env.Vinfdot(t)^2 + (env.RSdot(t)*rvec[j])^2) #Update Wdotvec
            idx = 2na + j
            pds[idx] = ((twistvec[j] + pitch) - phivec[j])
            if turbine
                pds[idx] *= -1
            end
        end
end


function extractloads(dsmodel::DS.Riso, x, ccout, t, rvec, chordvec, twistvec, pitch, blade::Blade, env::Environment) #TODO: This should probably be an inplace function. 
    n = length(blade.airfoils)

    Cl = Array{eltype(chordvec)}(undef, n)
    Cd = Array{eltype(chordvec)}(undef, n)
    Cn = Array{eltype(chordvec)}(undef, n)
    Ct = Array{eltype(chordvec)}(undef, n)
    N = Array{eltype(chordvec)}(undef, n)
    T = Array{eltype(chordvec)}(undef, n)


    for i = 1:n
        idx = 4*(i-1)
        xs = x[1+idx:idx+4]

        u = ccout.W[i]
        # udot = sqrt(env.Vinfdot(t)^2 + (env.RSdot(t)*rvec[i])^2)
        alpha = -((twistvec[i] + pitch) - ccout.phi[i]) #TODO: Make this work for a turbine or a propeller. 
        alphadot = 0.0

        # ys = [u, udot, 0.0, 0.0, theta, 0.0]
    
        Cl[i], Cd[i] = DS.riso_coefficients(xs, u, alpha, alphadot, chordvec[i], blade.airfoils[i])

        cphi = cos(ccout.phi[i])
        sphi = sin(ccout.phi[i])

        Cn[i] = Cl[i]*cphi - Cd[i]*sphi
        Ct[i] = Cl[i]*sphi + Cd[i]*cphi

        N[i] = Cn[i]*0.5*env.rho*ccout.W[i]^2*chordvec[i]
        T[i] = Ct[i]*0.5*env.rho*ccout.W[i]^2*chordvec[i]
    end
    return N, T, Cn, Ct, Cl, Cd
end


function parsesolution(dsmodel::DS.Riso, xds, W, phi, tvec, chordvec, twistvec, pitch, blade::Blade)
    m = length(tvec)
    n = length(blade.airfoils)

    Cl = zeros(m, n)
    Cd = zeros(m, n)

    for j = 1:m
        ut = xds[j,:]

        for i = 1:n
            idx = 4*(i-1)
            xs = ut[1+idx:idx+4]

            u = W[i]
            # udot = sqrt(env.Vinfdot(tvec[j])^2 + (env.RSdot(tvec[j])*rvec[i])^2)
            alpha = -((twistvec[i] + pitch) - phi[i]) 
            alphadot = 0.0

    
            Cl[j, i], Cd[j, i] = DS.riso_coefficients(xs, u, alpha, alphadot, chordvec[i], blade.airfoils[i])
        end
    end
    return Cl, Cd
end














