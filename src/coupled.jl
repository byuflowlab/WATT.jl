
struct RisoState{TF}
    x::Array{TF, 1}
end


function createrisoode(blade, chordvec, twistvec, pitch, phivec, Vxvec, Vyvec, Vxdotvec, Vydotvec)
    risoODE = function(x, p, t)
        n = length(blade.airfoils)
        dx = zeros(4*n)
        for i = 1:n
            idx = 4*(i-1)
            xs = x[1+idx:idx+4]

            # Vx = 1 #Do I want to pass this in? Or do I calculate it here. I think I'll pass it in. 
            # Vy = 1

            u = sqrt(Vxvec[i]^2 + Vyvec[i]^2) #Condense the inflow velocity to a single value.
            v = 0 

            udot = sqrt(Vxdotvec[i]^2 + Vydotvec[i]^2)  
            vdot = 0.0

            theta = -((twistvec[i] + pitch) - phivec[i]) #Calculate the angle of attack that the airfoil see
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

#=
    Vxvec = [env.Vinf(t0) for i in 1:na]
    Vdotvec = zeros(na)
    risoode = createrisoode(blade, chordvec, twistvec, pitch, ccout.phi, Vxvec, Vyvec, Vdotvec, Vdotvec)
    x0 = zeros(4*na) #TODO: I think that one of these states might need to be initialized as phi.
    x0[4:4:4*na] .= 1.0 

    prob = SteadyStateProblem(risoode, x0)
    sol = DifferentialEquations.solve(prob)

=#


function simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, t0, assembly)

    ### Create CCBlade inputs
    rotor = CCBlade.Rotor(rhub, rtip, 1.0, precone=0.0, turbine=true)
    na = length(blade.airfoils)
    sections = [Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i in 1:na]

    ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
    Vx = env.Vinf(t0)
    Vyvec = [env.RS(t0)*rvec[i] for i in 1:na]
    pitch = 0.0 #TODO: I'm not sure if this pitch does something or not. And I don't know if I want to pitch yet. 
    operatingpoints = [CCBlade.OperatingPoint(Vx, Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na] 
    #TODO: Dr. Ning uses mu=1 and a=1 in his example. 
    #TODO: Might get rid of the gxmodel, bemmodel, and dsmodel as pass ins. 

     
    ### Solve BEM residual
    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

    



    ### Get steady state Riso states
    #I need a way to solve the Riso states. I'd love to use DifferentialEquations. But... I'm not quite sure how I'll pass in all the information that I want. Additionally. I don't want to create a function everytime I call this function. That'll slow things down. So I think it'd be best to use a functor. Which means that I'll need to rearrange what goes in the Riso object... or use another object. Nothing is inside the Riso object currently, so it wouldn't be bad to put something in there... That or I could use a massive p vector like Taylor did. And just suck it up. I can rearrange what the heck I put in there... because currently it isn't super useful. The final option is I could use my own solver. Then I can pass in what ever the heck I'd like. 

    #It might be better to use the p vector... because the twist, phi, and other inputs might change every iteration.... which I guess is problematic if I only get to pass in x, p, and t. .... hmmm. 

    #Okay, so I did a quick little test to see if creating a method on a struct is faster than creating a function. It is faster to create a struct than it is to create the function (by 100x in my test.. but I don't know how that scales.. but we're talking 1.5 ns vs 100 ns.). But it was faster to solve the ode problem using the compiled function rather than the method on a struct (by like 8%, but I don't know how that scales either. It seems like it stays the same order of magnitude). I guess in this case... every iteration I'd be creating a new instance. At least the way I'm approaching it now. But with the fact that the function doesn't take that much more time to compile, then I can run that for now (especially since I might have that already). 

    Vxvec = [env.Vinf(t0) for i in 1:na]
    Vdotvec = zeros(na)
    risoode = createrisoode(blade, chordvec, twistvec, pitch, ccout.phi, Vxvec, Vyvec, Vdotvec, Vdotvec)
    x0 = zeros(4*na) #TODO: I think that one of these states might need to be initialized as phi.
    x0[4:4:4*na] .= 1.0 

    prob = SteadyStateProblem(risoode, x0)
    sol = DifferentialEquations.solve(prob)




    ### Create GXBeam State & constant inputs -> Might be a static solve.... I need to see what options are available. 
    ## Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed 

    ## Create point masses
    point_masses = Dict(1 => PointMass(0.0, SVector(0.0, 0.0, 0.0), @SMatrix zeros(3,3)))

    ### Extract CCBlade Loads and create a distributed load
    Fzfit = Akima(rvec, -ccout.Np )
    Fyfit = Akima(rvec, ccout.Tp)

    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]
    distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(rgx[ielem]), fz= (s) -> Fzfit(rgx[ielem])) for ielem in 1:nelem) #TODO: Could make this a fit that I apply... I think.. 
    

    Omega = SVector(0.0, 0.0, omega)

    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega)

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)


    # ### Iterate through time steps
    # for i = 1:nt
        ### Solve Riso states based on inflow and deflection
        #Todo: Make sure that my ODE solver works




        ### Create changing CCBlade inputs (based on GXBeam state)
        #Todo: Figure out how to create a function that returns the lift and drag based on the current DS states




        ### Solve CCBlade residual 




        ### Update GXBeam loads



        ### Solve GXBeam for time step




    # end

    return ccout, state
end