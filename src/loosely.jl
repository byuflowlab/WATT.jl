#=

Recreating what they do in AeroDyn, which is solve the BEM, feed the inflow angle to the dynamic stall model, then calculate the loading. 

Adam Cardoza 8/6/22
=# 

function createrisoode(blade)  
    risoODE = function(x, p, t)
        ### Unpack inputs. 
        n = length(blade.airfoils)

        #p = [chordvec, twistvec, phivec, Wvec, Wdotvec, pitch]
        chordvec = view(p, 1:n)
        twistvec = view(p, n+1:2*n)
        phivec = view(p, 2n+1:3n)
        Wvec = view(p, 3n+1:4n)
        Wdotvec = view(p, 4n+1:5n)
        pitch = p[end]


        ### Iterate through the nodes and calculate the state rates. 
        dx = zeros(4*n)

        for i = 1:n
            idx = 4*(i-1)
            xs = x[1+idx:idx+4]

            u = Wvec[i] 
            v = 0 

            udot = Wdotvec[i]  
            vdot = 0.0

            theta = -((twistvec[i] + pitch) - phivec[i]) 

            # if i==n
            #     @show theta, u, udot, t
            # end

            thetadot = 0.0 #Assuming that the airfoil isn't in the act of turning???? TODO: What are my other options?

    
            dx[1+idx:4+idx] = riso_states(xs, u, udot, v, vdot, theta, thetadot, chordvec[i], blade.airfoils[i])

            if isnanvec(dx[1+idx:4+idx])
                println("Riso State Rates: ", [dx[1+idx], dx[2+idx], dx[3+idx], dx[4+idx]])
                error("Nan in Riso")
            end

        end
        return dx
    end
    return risoODE
end 

abstract type DSmodelInit end

struct Hansen <: DSmodelInit
end

function initializeDSmodel(dsmodelinit::Hansen, nt, na)
    xds = zeros(nt, 4*na) #Note: I thought about initializing this as an undefined array, but... the first vector needs to be mostly zeros, and I feel like this will be the easiest. 
    xds[1,4:4:4*na] .= 1
    return xds
end

function extractloads(x, ccout, t, rvec, chordvec, twistvec, pitch, blade::Blade, env::Environment) #TODO: This should probably be an inplace function. 
    n = length(blade.airfoils)

    Cl = Array{eltype(rvec)}(undef, n)
    Cd = Array{eltype(rvec)}(undef, n)
    Cn = Array{eltype(rvec)}(undef, n)
    Ct = Array{eltype(rvec)}(undef, n)
    N = Array{eltype(rvec)}(undef, n)
    T = Array{eltype(rvec)}(undef, n)


    for i = 1:n
        idx = 4*(i-1)
        xs = x[1+idx:idx+4]

        u = ccout.W[i]
        udot = sqrt(env.Vinfdot(t)^2 + (env.RSdot(t)*rvec[i])^2)
        theta = -((twistvec[i] + pitch) - ccout.phi[i]) #TODO: Make this work for a turbine or a propeller. 

        ys = [u, udot, 0.0, 0.0, theta, 0.0]
    
        Cl[i], Cd[i] = riso_coefs(xs, ys, chordvec[i], blade.airfoils[i])

        cphi = cos(ccout.phi[i])
        sphi = sin(ccout.phi[i])

        Cn[i] = Cl[i]*cphi - Cd[i]*sphi
        Ct[i] = Cl[i]*sphi + Cd[i]*cphi

        N[i] = Cn[i]*0.5*env.rho*ccout.W[i]^2*chordvec[i]
        T[i] = Ct[i]*0.5*env.rho*ccout.W[i]^2*chordvec[i]
    end
    return N, T, Cn, Ct, Cl, Cd
end


function simulate(rvec, chordvec, twistvec, blade::Blade, env::Environment, tvec, rhub, rtip, B, precone; turbine::Bool=true, dsmodelinit::DSmodelInit=Hansen(), solver::Solver=RK4(), verbose::Bool=false)

    if verbose
        println("Rotors.jl initializing solution...")
    end

    ### Initialization information
    na = length(rvec)
    nt = length(tvec)

    t0 = tvec[1]


    ### Initialize BEM solution
    rotor = Rotor(rhub, rtip, B; precone, turbine)
    sections = [CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i = 1:na]

    Vxvec = [env.Vinf(t0) for i in 1:na]
    Vyvec = [env.RS(t0)*rvec[i]*cos(precone) for i in 1:na]
    # @show Vxvec
    # @show Vyvec
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na]

    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)
    # @show ccout.phi

    ### Initialize DS solution
    risoode = createrisoode(blade)

    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i])^2) for i in 1:na]

    p_ds = vcat(chordvec, twistvec, ccout.phi, ccout.W, Wdotvec, pitch) #p = [chordvec, twistvec, phivec, Wvec, Wdotvec, pitch]


    xds = initializeDSmodel(dsmodelinit, nt, na)


    ### Prepare data storage
    cchistory = Array{Array{Outputs{eltype(rvec)}, 1}, 1}(undef, nt) #[ccout for i = 1:nt]
    Cl = Array{eltype(rvec)}(undef,(nt, na))
    Cd = Array{eltype(rvec)}(undef,(nt, na))
    Cn = Array{eltype(rvec)}(undef,(nt, na))
    Ct = Array{eltype(rvec)}(undef,(nt, na))
    N = Array{eltype(rvec)}(undef,(nt, na))
    T = Array{eltype(rvec)}(undef,(nt, na))

    cchistory[1] = ccout
    N[1,:], T[1,:], Cn[1,:], Ct[1,:], Cl[1,:], Cd[1,:] = extractloads(xds[1,:], cchistory[1], t0, rvec, chordvec, twistvec, pitch, blade, env)
    

    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        ### Update BEM inputs
        for j = 1:na
            Vxvec[j] = env.Vinf(t)
            Vyvec[j] = env.RS(t)*rvec[j]*cos(precone)
            operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
        end

        ### Solve BEM
        cchistory[i] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)


        ### Update Dynamic Stall model inputs
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
        if verbose
            println("Simulation time: ", t)
        end
    end
    return N, T, cchistory, xds
end



function parsesolution(xds, cchistory, tvec, rvec, chordvec, twistvec, pitch, blade::Blade, env::Environment)
    m = length(tvec)
    n = length(blade.airfoils)

    Cl = zeros(m, n)
    Cd = zeros(m, n)

    for j = 1:m
        ut = xds[j,:]

        for i = 1:n
            idx = 4*(i-1)
            xs = ut[1+idx:idx+4]

            u = cchistory[j].W[i]
            udot = sqrt(env.Vinfdot(tvec[j])^2 + (env.RSdot(tvec[j])*rvec[i])^2)
            theta = -((twistvec[i] + pitch) - cchistory[j].phi[i]) 

            ys = [u, udot, 0.0, 0.0, theta, 0.0]
    
            Cl[j, i], Cd[j, i] = riso_coefs(xs, ys, chordvec[i], blade.airfoils[i])
        end
    end
    return Cl, Cd
end