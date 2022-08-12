#=

A dynamic BEM - dynamic stall (DS) model solver, using a time-marching approach. Option for fixed-point iteration at every time step. DS model is a indicial approach.  

My goal is to include everything I need inside this file, except for blade.jl and environments.jl. I'll move stuff out later. -> I need more space and organization. I'm starting some new files. 

Adam Cardoza 8/4/22
=# 




function initializeBEM(rvec, chordvec, twistvec, blade, env, t0, rhub, rtip, B, precone, turbine)
    rotor = Rotor(rhub, rtip, B; precone, turbine)
    sections = [CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i = 1:na]

    Vxvec = [env.Vinf(t0) for i in 1:na]
    Vyvec = [env.RS(t0)*rvec[i]*cos(precone) for i in 1:na]
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na]

    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

    return ccout.phi, Vxvec, Vyvec #Not passing out W because I want the inflow directly into the dynamic stall model. 
end

function initializeDS()
end


function simulate(rvec, chordvec, twistvec, blade, env, tvec, rhub, rtip, B, precone; inneriterations::Int=3, turbine::Bool=true)

    ### Initialization information
    na = length(rvec)
    nt = length(tvec)

    t0 = tvec[1]

    phi = Array{eltype(rvec)}(undef, (nt, na))

    ### Initialize BEM solution
    phi[1,:], Vxvec, Vyvec = initializeBEM(rvec, chordvec, twistvec, blade, env, t0, rhub, rtip, B, precone, turbine)


    ### Initialize DS solution


    ### Prepare data storage


    ### Iterate through time 
    for i = 2:nt
        t = tvec[i-1]
        dt = tvec[i] - tvec[i-1]

        ### Fixed point iteration
        for j = 1:inneriterations

        end
    end

end
