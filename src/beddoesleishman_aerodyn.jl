#=
Code to interact with the DynamicStallModels Package, specifically for the Beddoes-Leishman model, the AeroDyn version. 

=#



function initializeBLA(dsmodel::DS.BeddoesLeishman, dsmodelinit::BeddoesLeishman, solver::Solver, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

    if turbine
        thetavec = [-((twistvec[i] + pitch) - phivec[i]) for i = 1:na] #TODO: I wonder if I should do this in the simulate section. 
    else
        thetavec = [((twistvec[i] + pitch) - phivec[i]) for i = 1:na]
    end

    ns = DS.numberofstates(dsmodel)
    # @show ns #This looks correct. 
    xds = Array{eltype(chordvec), 2}(undef, nt, ns)
    pds = Array{eltype(chordvec), 1}(undef, 18*na)

    println("Initializing BLA BL init")

    for i = 1:na
        xidx = 22*(i-1)+1:i*22
        pidx = 18*(i-1)+1:i*18
        # @show xidx, pidx 
        xds[1,xidx], _, pds[pidx] = DS.initialize_ADO([Vxvec[i]], [thetavec[i]], tvec, dsmodel.airfoils[i], chordvec[i], a) #It appears that I'm inputing the correct p vector. 
    end

    # @show xds #This looks correct. 
    return dsmodel, xds, pds
end

function initializeBLA(dsmodel::DS.BeddoesLeishman, dsmodelinit::Steady, solver::Solver, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) #TODO: I don't think that I need to pass na and nt in. na is found in dsmodel. nt is found from the length of tvec. 
    @warn("The steady initialization is not yet prepared, doing the BeddoesLeishman initialization.")

    if turbine
        thetavec = [-((twistvec[i] + pitch) - phivec[i]) for i = 1:na] 
    else
        thetavec = [((twistvec[i] + pitch) - phivec[i]) for i = 1:na]
    end

    xds = Array{eltype(chordvec), 2}(undef, nt, DS.numberofstates(dsmodel))
    pds = Array{eltype(chordvec), 1}(undef, 18*na)
    for i = 1:na
        xidx = 22*(i-1)+1:i*22
        pidx = 18*(i-1)+1:i*18
        xds[1,xidx], _, pds[pidx] = DS.initialize_ADO([thetavec[i]], tvec, dsmodel.airfoils[i], chordvec[i], a)  #Todo: Broken. 
    end

    # @show xds
    return dsmodel, xds, pds
end



function update_BLA_parameters!(dsmodel::DS.BeddoesLeishman, turbine::Bool, pds, na, rvec, Wvec, phivec, twistvec, pitch, env, t)
        for j = 1:na
            ### Update inflow velocity
            idx = 18*(j-1) + 17 #TODO: There is a function is DSM that does most of what this function is doing, except calculate the angle of attack from the inflow angle... should I use that function here? 
            pds[idx] = Wvec[j]

            ### Update angle of attack. 
            idx = 18*(j-1) + 18
            pds[idx] = ((twistvec[j] + pitch) - phivec[j])
            if turbine
                pds[idx] *= -1
            end
        end
end


function extractloads_BLA(dsmodel::DS.BeddoesLeishman, x, ccout, chordvec, twistvec, pitch, blade::Blade, env::Environment) #TODO: This should probably be an inplace function. #TODO: The blade seems like it has repetitive data, the airfoils and radial node locations.... So I either need to rely on it more, or just get rid of it. 
    n = dsmodel.n

    # println("Using BLA extract loads")

    Cl = Array{eltype(chordvec)}(undef, n)
    Cd = Array{eltype(chordvec)}(undef, n)
    Cn = Array{eltype(chordvec)}(undef, n)
    Ct = Array{eltype(chordvec)}(undef, n)
    Cx = Array{eltype(chordvec)}(undef, n)
    Cy = Array{eltype(chordvec)}(undef, n)

    # @show x
    for i = 1:n
        idx = 22*(i-1)
        xs = x[1+idx:idx+22] 
        af = dsmodel.airfoils[i]

        u = ccout.W[i] 
        # alpha = -((twistvec[i] + pitch) - ccout.phi[i]) #TODO: Make this work for a turbine or a propeller. 
        # if i==2
        #     @show alpha #This appears to be slightly off. -> It could be off cause I might be starting from a different azimuthal angle. 
        # end

        Cn[i], Ct[i], Cl[i], Cd[i], _ = DS.BLAD_coefficients(dsmodel, xs, u, chordvec[i], af, env.a)

        # if i==n
        #     @show u, alpha, Cn[i], Ct[i]
        # end

        cphi = cos(ccout.phi[i])
        sphi = sin(ccout.phi[i])
        Cx[i] = Cl[i]*cphi + Cd[i]*sphi
        Cy[i] = Cl[i]*sphi - Cd[i]*cphi

        # N[i] = Cn[i]*0.5*env.rho*u^2*chordvec[i] #Todo. Is this going to need to be dimensionalized by the actual velocity (including induced velocities)? -> yes. yes it is. And ccout.W isn't the actual inflow velocity, it is just the sum of the.... wait.... it is... Line 304 of CCBlade.jl calculates the actual inflow velocity. 
        # T[i] = Ct[i]*0.5*env.rho*u^2*chordvec[i]
    end
    return Cx, Cy, Cn, Ct, Cl, Cd
end


# function parsesolution(xds, W, phi, tvec, chordvec, twistvec, pitch, blade::Blade)
#     m = length(tvec)
#     n = length(blade.airfoils)

#     Cl = zeros(m, n)
#     Cd = zeros(m, n)

#     for j = 1:m
#         ut = xds[j,:]

#         for i = 1:n
#             idx = 4*(i-1)
#             xs = ut[1+idx:idx+4]

#             u = W[i]
#             # udot = sqrt(env.Vinfdot(tvec[j])^2 + (env.RSdot(tvec[j])*rvec[i])^2)
#             alpha = -((twistvec[i] + pitch) - phi[i]) 
#             alphadot = 0.0

    
#             Cl[j, i], Cd[j, i] = DS.riso_coefficients(xs, u, alpha, alphadot, chordvec[i], blade.airfoils[i])
#         end
#     end
#     return Cl, Cd
# end
