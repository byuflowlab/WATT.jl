#=
This is an attempt at not tightly coupling the models. The tight coupling is proving difficult to solve, especially with the damping as part of the structural model. So we're just going to scale back. Not tightly coupled. Just passing loads and deflections every time step. Maybe in the future I'll allow for some interior iterations (in a fixed point solution). 





##### Work Notes
- I just had the thought. I should go simpler than this. I should do a solution of just BEM to GXBeam across time. Simple one-way coupling. Then do a two-way coupling. Then add in the dynamic stall model. Then add in the possibility of having multiple interior iterations.

=#


struct RisoState{TF}
    x::Array{TF, 1}
    u::TF
    c::TF
    airfoil::Airfoil
end

function (state::RisoState)(alpha, Re , Mach)
    y = [state.u, 0.0, 0.0, 0.0, alpha, 0.0]
    return riso_coefs(state.x, y, state.c, state.airfoil) #TODO: I don't know if the fact that this function returns a Static Vector will throw a wrench in the works, but we'll just have to see... won't we? 
end

function afeval(state::RisoState, alpha, Re, Mach)
    return state(alpha, Re, Mach)
end

"""
createrisoode(blade)

Create the riso ode for this simulation. 

### Inputs
- blade::Blade - the blade description

### Outputs
- risoODE::Function - a function that calculates the state rates for the riso model. 

### Notes
- The states (x) are the four riso states for each node concatonated together. [point 1 state 1-4, point 2 state 1-4, point 3.... ]. 
- The parameters (p) for this function are the following vectors concatonated together: chordvec, twistvec, phivec, Vxvec, Vyvec, Vxdotvec, Vydotvec, pitch. 
"""
function createrisoode(blade)  
    risoODE = function(x, p, t)
        ### Unpack inputs. 
        n = length(blade.airfoils)

        chordvec = view(p, 1:n)
        twistvec = view(p, n+1:2*n)
        phivec = view(p, 2n+1:3n)
        Wvec = view(p, 3n+1:4n)
        Wdotvec = view(p, 4n+1:5n)
        # Vxvec = view(p, 3n+1:4n)
        # Vyvec = view(p, 4n+1:5n)
        # Vxdotvec = view(p, 5n+1:6n)
        # Vydotvec = view(p, 6n+1:7n)
        pitch = p[end]


        ### Iterate through the nodes and calculate the state rates. 
        dx = zeros(4*n)
        for i = 1:n
            idx = 4*(i-1)
            xs = x[1+idx:idx+4]

            u = Wvec[i] #sqrt(Vxvec[i]^2 + Vyvec[i]^2) #Condense the inflow velocity to a single value.
            v = 0 

            udot = Wdotvec[i] #sqrt(Vxdotvec[i]^2 + Vydotvec[i]^2)  
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

"""
update_p_a!(p_a, deformed_twist)

Update the aerodynamic parameters that change from iteration to iteration (for Riso). 

### Inputs

### Outputs
- p_a::Array{TF, 1} - updated vector of aerodynamic parameters for Riso. 

### Notes
- Recall that the order of p_a is: [chordvec, twistvec, phivec, Vxvec, Vyvec, Vxdotvec, Vydotvec, pitch]. Changed to vcat(chordvec, twistvec, phivec, Wvec, Wdotvec, pitch)
"""
function update_p_a!(p_a, deformed_twist, phivec, Wvec, Wvecdot; n=Int(length(deformed_twist)-1)/5)
    p_a[n+1:2*n] = deformed_twist
    p_a[2n+1:3n] = phivec
    p_a[3n+1:4n] = Wvec
    p_a[4n+1:5n] = Wvecdot
end




"""
find_paired_points(points, r)

Finds the two structural points before and after the given aerodynamic node. 

### Inputs
- points - the points vector from the assembly object. Each point is a SVector{3}. 
- r - the radial location of the aerodynamic point of interest

### Outputs
- pair - a tuple containing the indices of the point before and after the aerodynamic point (before, after)
"""
function find_paired_points(points, r)
    ngx = length(points)
    rgx = [points[i][1] for i = 1:ngx]
    for i = 1:length(points)-1
        if rgx[i]>r
            return (i-1, i)
        end
    end
    return (ngx-1, ngx)
end

function find_interpolation_value(points, pair, r)
    ngx = length(points)
    rgx = [points[i][1] for i = 1:ngx]
    a = r - rgx[pair[1]]
    b = rgx[pair[2]] - r
    L = a + b
    return a/L
end

"""
InterpolationPoint(pair, percent)

A struct to quickly interpolate the new radial location of an aerodynamic node. 
"""
struct InterpolationPoint
    pair::Tuple{Int64, Int64}
    percent::Float64
end

"""
create_interpolationpoints(assembly, rvec)

### Inputs
- assembly - A GXBeam assembly
- rvec - The radial location of the aerodynamic nodes

### Outputs
- interpvals

### Notes
- This assumes that the structural points are in order and go from root to tip. Additionally, the aerodynamic points follow a similar order. 
"""
function create_interpolationpoints(assembly, rvec)
    na = length(rvec)
    pairs = [find_paired_points(assembly.points, rvec[i]) for i = 1:na]
    percents = [find_interpolation_value(assembly.points, pairs[i], rvec[i]) for i =1:na]

    return [InterpolationPoint(pairs[i], percents[i]) for i = 1:na]
end


"""
interpolate_position(ip, assembly, state)

Get the position of aerodynamic nodes from the interpolation struct and the structural structs. 
"""
function interpolate_position(ip, assembly, state)  
    rgx = [assembly.points[i][1] + state.points[i].u[1] for i = 1:length(assembly.points)]
    return (1-ip.percent)*rgx[ip.pair[1]] + ip.percent*rgx[ip.pair[2]]
end

"""
interpolate_velocity(ip, assembly, state, env)

Get the relative velocity of an aerodynamic point. 
"""
function interpolate_velocity(ip, assembly, state) #Note: This is probably a very poor approximation. But It'll have to do for now. Maybe there is a better way to get the velocity of the points. 
    ne = length(assembly.elements)
    np = ne + 1

    V = state.elements[ip.pair[1]].V  #Get the velocity of the element 
    Omega = state.elements[ip.pair[1]].Omega
    r1 = assembly.points[ip.pair[1]] - assembly.elements[ip.pair[1]].x #The position vector from the element point to the starting point
    r2 = assembly.elements[ip.pair[1]].x - assembly.points[ip.pair[2]] #The position vector from the element point to the stopping point
    
    v1 = V + cross(Omega, r1)
    v2 = V + cross(Omega, r2)

    # if ip.pair[2] == 18
    #     @show V #V is what is making the velocity go super wonky. 
    # end

    return (1-ip.percent)*v1 + ip.percent*v2
end




"""
simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, assembly, tspan, dt; solver=RK4())

Simulates the described wind turbine across the given time span using a BEM (I am using CCBlade.) that is extended by the Riso dynamic stall method. A structural response is simulated by GXBeam. The BEM is assumed to be constant across a given time step. The BEM is influenced by the current state of the dynamic stall method. The dynamic stall and structural models are then integrated forward based on the BEM solution. The dynamic stall model is integrated using the provided solver (defaults to a Runge-Kutta 4th order) and the structural model is solved using it's internal dynamic solver. 

### Inputs
- rvec::Array{TF, 1} - A vector holding the aerodynamic 

### Outpus


### Notes
- The structural model does not include damping. 
- The initial condition is assumed to be the loads from the undeflected solution (to add the option to start from the fixed point iteration solution later). Those loads are then used to find the steady state deflections and inflow velocities. The inflow velocities are used to calculate the steady state states for the Riso model.
"""
function simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, assembly, tspan, dt; solver=RK4(), verbose=false)

    t0 = tspan[1]
    tvec = tspan[1]:dt:tspan[2]
    nt = length(tvec)

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

     
    ### Solve BEM residual
    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

    



    ### Get steady state Riso states
    #I need a way to solve the Riso states. I'd love to use DifferentialEquations. But... I'm not quite sure how I'll pass in all the information that I want. Additionally. I don't want to create a function everytime I call this function. That'll slow things down. So I think it'd be best to use a functor. Which means that I'll need to rearrange what goes in the Riso object... or use another object. Nothing is inside the Riso object currently, so it wouldn't be bad to put something in there... That or I could use a massive p vector like Taylor did. And just suck it up. I can rearrange what the heck I put in there... because currently it isn't super useful. The final option is I could use my own solver. Then I can pass in what ever the heck I'd like. 

    #It might be better to use the p vector... because the twist, phi, and other inputs might change every iteration.... which I guess is problematic if I only get to pass in x, p, and t. .... hmmm. 

    #Okay, so I did a quick little test to see if creating a method on a struct is faster than creating a function. It is faster to create a struct than it is to create the function (by 100x in my test.. but I don't know how that scales.. but we're talking 1.5 ns vs 100 ns.). But it was faster to solve the ode problem using the compiled function rather than the method on a struct (by like 8%, but I don't know how that scales either. It seems like it stays the same order of magnitude). I guess in this case... every iteration I'd be creating a new instance. At least the way I'm approaching it now. But with the fact that the function doesn't take that much more time to compile, then I can run that for now (especially since I might have that already). 

    Vxvec = [env.Vinf(t0) for i in 1:na]
    # Wvec = [sqrt(Vxvec[i]^2 + Vyvec[i]^2) for i = 1:na]
    # Wdotvec = zeros(na)

    # risoode = createrisoode(blade)

    # p_a = vcat(chordvec, twistvec, phivec, Wvec, Wdotvec, pitch) #Create p_a #Todo: I don't think that I've done the twist & inflow relationship correctly for the dynamic stall model. 

    # x0 = zeros(4*na) #TODO: I think that one of these states might need to be initialized as phi.
    # x0[4:4:4*na] .= 1.0 

    # prob = SteadyStateProblem(risoode, x0, p = p_a)
    # sol = DifferentialEquations.solve(prob)

    # dsstates = [RisoState(sol[4*(i-1)+1:4*(i-1)+4], Wvec[i], chordvec[i], blade.airfoils[i]) for i = 1:na]


    #### I'm going to use my own solver for the meantime. 
    # odeprob = ODEProblem(risoode, sol.u, tspan, p_a) #Note: I think this is going to be bad and either soak up memory or time. I'll probably have to switch to an in house integrator to minimize costs. Who knows though. It'd be nice to have 
    # integrator = init(odeprob, Tsit5())

    # x_old = sol.u






    ### Create GXBeam State & constant inputs -> Might be a static solve.... I need to see what options are available. 
    ## Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed 

    ## Create point masses
    # point_masses = Dict(1 => PointMass(0.0, SVector(0.0, 0.0, 0.0), @SMatrix zeros(3,3)))

    ### Extract CCBlade Loads and create a distributed load
    Fzfit = Akima(rvec, -ccout.Np )
    Fyfit = Akima(rvec, ccout.Tp)

    # plt = plot(rvec, ccout.Np, lab="Normal", legend=:topleft) #This plot looks fairly normal. I don't know why gxbeam is having such trouble with it. 
    # plot!(rvec, ccout.Tp, lab="Tangent")
    # display(plt)

    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]
    distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(rgx[ielem]), fz= (s) -> Fzfit(rgx[ielem])) for ielem in 1:nelem) #TODO: Could make this a fit that I apply... I think.. #Todo. What type of load is this? A follower load? -> It's not a follower load. 

    # @show distributed_loads[nelem]
    
    # @show distributed_loads

    Omega = SVector(0.0, 0.0, env.RS(t0))

    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega) #I don't know why, but this is taking a really long time (1.6 seconds w/ 106.6k allocations). Additionally, the deflections were super wonky. -> Maybe I try it without the damping. Will that break my create_gxbeam_assembly??? 
    # It took 1.06 seconds with 72.66k allocations with the main branch of GXBeam.... which isn't great... but is better. 
    # The first run with damping took 123.46 seconds with 116.2 M allocations (98.85% compilation time). And again... the deflections are wonky. 
    #Well... heck. Options? - bail on damping, go forward. Either add damping later... or Fork GXBeam and add it later. - Figure out external damping now. - Figure out internal damping now.  -> Well... I think I'm going to go with bail on damping... and then I'll add it in later. I want results sooner... and I'm tired of GXBeam not working. 

    gxstate = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)



    cchistory = [ccout for i = 1:nt]
    gxhistory = [gxstate for i = 1:nt]

    # def_thetax = [gxstate.elements[ielem].theta[1] for ielem = 1:n]

    interpolationpoints = create_interpolationpoints(assembly, rvec) #Todo. I might have just made this and it might not be useful. ... Wait I'm using it to interpolate the aerodynamic load onto the structural... dang... this doesn't work as well... well.. at least I have an r vector that corresponds with the aerodynamic loads. Now I just need to create a function to create the distributed loads from a fit. (Could easily create the fit inside the function as well)
    aeroV = [interpolate_velocity(interpolationpoints[i], assembly, gxstate) for i = 1:na]

    ### Iterate through time steps
    for i = 1:nt-1
        ## Update timestep
        t = tvec[i]


        ## Update CCBlade inputs (based on GXBeam state)
        #Todo. Figure out how to create a function that returns the lift and drag based on the current DS states. I'm going to create a struct that holds the local riso states, then I'll create a method on that struct and extend afeval. I might need to stick a lot more than just the state in the struct in order to evaluate.... We'll just go with it for now. 
        # twist_displaced = twistvec #.- def_thetax

        rotor = CCBlade.Rotor(rhub, rtip, 1.0, precone=0.0, turbine=true) #TODO: We use these values for tip corrections... when the beam deflects, then they nodes will be farther from these values... which will decrease tip losses. If I just shrink the tip value then the more center region will feel the effect of tip losses creater than they normally would. How do I approach this? -> Maybe it doesn't matter. I don't change the section radial value... then it might not make a difference in the tip correction. The velocity is already getting corrected, so I should be good. 

        # for i = 1:na   #Note: Don't need to update sections while 
        #     sections[i] = Section(rvec[i], chordvec[i], twist_displaced[i], dsstates[i])
        # end

        # just finished the dsstates struct, now I need to update the deflected velocity. Which means that I need to figure out where the aerodynamic nodes are in relationship to the structural nodes. (Man this would be so much easier if they were just colocated...) Maybe I'll create a vector before this that gives the percentage distance between structural nodes that the aerodynamic nodes are.... but that's assuming there are as many structural nodes as there are aerodynamic loads... maybe I'll brainstorm with Tyler. 

        # Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
        
        Vxvec .= [env.Vinf(t) for i in 1:na] #[env.Vinf(t) + aeroV[i][3] for i in 1:na] #Note: I'm not sure that this is actually going to save allocations. 
        Vyvec .= [env.RS(t)*rvec[i] for i in 1:na]  #[-aeroV[i][2] for i in 1:na] 
        operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na] 
        #TODO: Dr. Ning uses mu=1 and a=1 in his example. 


        # Solve BEM residual
        cchistory[i+1] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

        # @show ccout.Np[end] #Stays the same, as expected. 

        #### I need to revamp the risoode so that all of the things that change come through the p vector so that I don't have to recreate the problem every time step. -> Revamped the function. It needs to be tested. Additionally, I need to figure out a way to update p_a. -> I created a function that updates the p_a vector every time step. 

        ## Solve Riso states based on inflow and deflection
        # Wvec .= [sqrt(Vxvec[i]^2 + Vyvec[i]^2) for i = 1:na]
        # Wdotvec .= 0.0 #Note: I think this will need to change, but I'm not 100% how to get the acceleration out of the structural model. 

        # p_a = update_p_a!(p_a, twist_displaced, ccout.phi, Wvec, Wvecdot) #TODO: Update p_a
     
        # remake(odeprob, p=p_a) #Note. I think this might be costly. I'm not sure how to do this properly. Todo. Additionally, I need to make sure that the correct states gets passed to the integrator at every step... I think this is probably going to be simpler to use my own solvers... but I'd love to have all the algorithms from DifferentialEquations available. 
        # integrator = init(odeprob, Tsit5())
        # step!(integrator, dt, true)

        # x_new = solver(risoode, x_old, p_a, t, dt) 
        


        ### I need to add in the GXBeam stuffs now. Although, I just had the thought. I should go simpler than this. I should do a solution of just BEM to GXBeam across time. Simple one-way coupling. Then do a two-way coupling. Then add in the dynamic stall model. Then add in the possibility of having multiple interior iterations. 


        ## Update GXBeam loads
        Fzfit = Akima(rvec, -cchistory[i+1].Np ) #TODO: I need to plot the aerodynamic loads and the interpolated loads that are applied to GXBeam. 
        Fyfit = Akima(rvec, cchistory[i+1].Tp)
        # r_aero = [interpolate_position(interpolationpoints[i], assembly, state) for i in 1:na]
        # Fzfit = Akima(r_aero, -ccout.Np )
        # Fyfit = Akima(r_aero, ccout.Tp)

        rgx_deformed = [assembly.elements[i].x[1] for i in 1:nelem] #[assembly.elements[i].x[1] + state.elements[i].u[1] for i in 1:nelem] #Question: Will I ever need to actually deflect the aerodynamic nodes? Because I could always use the undeflected r_gx values to grab the loads from the BEM fit. So even though in reality those nodes have been deflected.... it wouldn't matter, because I know how the original structural nodes correlate to the original aerodynamic nodes. So... If I use that relationship, no matter how the structural nodes have moved, that should preserve the same relationship to the aerodynamic nodes. (At least fit-wise)
        distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(rgx_deformed[ielem]), fz= (s) -> Fzfit(rgx_deformed[ielem])) for ielem in 1:nelem)

        Omega = SVector(0.0, 0.0, env.RS(t))


        ## Solve GXBeam for time step
        #Question. Can I start a time simulation off of a steady state simulation system? -> Yes. Yes I can. 
        # println("sys 1: ", system.x[16:20])
        system, localhistory, converged = time_domain_analysis!(system, assembly, tvec[i:i+1]; prescribed_conditions=prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega, reset_state=false, initialize=false) #The system doesn't appear to be changing. 
        # println("sys 2: ", system.x[16:20])
        
        # @show length(localhistory) #Length is 2 as expected. 

        ## Extract GXBeam outputs
        # cchistory[i+1] = ccout #cchistory is constant... as expected (nothing should be changing it.)
        gxhistory[i+1] = localhistory[end] #Question. (Why is this bad boy not deflecting?) -> Does it happen to be since I've begun it from a steady state solution? -> That doesn't make sense... because it isn't deflected at all. I'd expect it to be at some deflected length... and stay there (since I'm starting from a steady state solution on GXBeam's half.). -> I wasn't initializing the beam properly. 

        for j = 1:na
            aeroV[j] = interpolate_velocity(interpolationpoints[j], assembly, localhistory[end])
        end
        # @show aeroV[end] #The first tip velocities are super small (9.3e-5 m/s)... Then it goes super large (2.0 e8 m/s)... then super large in the other direction (-4.66e8 m/s). Well on our way to lightspeed. 

        if verbose
            println("t=", round(t, digits=4))
        end
    end

    if verbose
        println("t=", round(tvec[end], digits=4))
    end

    return cchistory, gxhistory, tvec, gxstate
end

function simulate2(bemmodel, gxmodel, env, blade, p; maxiterations=1000, tolerance=1e-3, g=0.0, verbose=false)


    ### Seperate the parameters and read in any pertainent data
    n = gxmodel.ne
    pa = view(p, 1:(7*n))
    ps = view(p, (7*n+1):length(p))

    ### Constant problem dependent parameters
    pitch = pa[4]
    rhub = pa[5]
    rtip = pa[6]
    hubHt = pa[7]

    ### Environmental variables
    Omega = SVector(0.0, 0.0, -env.RS(0.0))
    grav = SVector(0.0, -g, 0.0)
    U = env.Vinf(0.0)


    ### Base vectors for turbine. 
    radii_idx = 1:7:(7*n)
    chord_idx = 2:7:(7*n)
    twist_idx = 3:7:(7*n)
    rvec = view(pa, radii_idx)  
    chordvec = view(pa, chord_idx)
    twistvec = view(pa, twist_idx)

    ### Unchanging constants
    t = 0.0
    precone = 0.0 #Precone. -> I'm going to change the radius manually. -> I could have a base precone, then allow for changes on top of that? 
    B = 1 #Number of blades
    yaw = 0.0
    tilt = 0.0
    azimuth = 0.0





    ### Create storage structures for solution
    resids = zeros(maxiterations, 3) #columns of resids is deflection y, force y, force z





    #### Create CCBlade inputs
    ### Create Rotor
    if bemmodel.tipcorrection
        rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bemmodel.turbine)
    else
        rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=bemmodel.turbine, tip=nothing)
    end

    # @show rvec
    # @show chordvec
    # @show twistvec

    ### Create Airfoils
    # airfoils = [CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", env.rho*U*rvec[i]/env.mu, U/env.a, Akima(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2]), Akima(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,3])) for i in 1:n] #Todo: Failing for some reason. 

    airfoils = [CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3]) for i in 1:n]

    ### Create Section
    sections = CCBlade.Section.(rvec, chordvec, twistvec, airfoils)

    ### Create Operating Point
    operatingpoints = CCBlade.windturbine_op.(U, env.RS(t), pitch, rvec, precone, yaw, tilt, azimuth, hubHt, bemmodel.shearexp, env.rho)



    #### Create GXBeam inputs
    ### Create Assembly
    assembly = create_gxbeam_assembly(gxmodel, ps)
    elements = view(assembly.elements, :) 


    ### Create Prescribed Conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))



    ##### Run initial solution
    ### Run CCBlade
    outs = CCBlade.solve.(Ref(rotor), sections, operatingpoints) 
    outshistory = [outs]

    ### Extract CCBlade Loads
    Fz = -outs.Np  
    Fy = outs.Tp

    ### Create distributed_load
    distributed_load = Dict{Int64, GXBeam.DistributedLoads{eltype(p)}}()
    f = @SVector zeros(3)
    m = @SVector zeros(3)
    m_follower = @SVector zeros(3)
    for i = 1:n #Iterate through the elements and apply the distributed load at every element. 
        f_follower = SVector(0.0, Fy[i]/elements[i].L, Fz[i]/elements[i].L) #Dividing by the length of the element so the force is distributed across the element. 
        distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end

    ### Run GXBeam
    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)
    history = [state]
    syshistory = [system]

    ### Obtain the deflections
    def_x = [history[1].elements[ielem].u[1] for ielem = 1:n]
    def_y = [history[1].elements[ielem].u[2] for ielem = 1:n]
    def_z = [history[1].elements[ielem].u[3] for ielem = 1:n]
    def_thetax = [history[1].elements[ielem].theta[1] for ielem = 1:n]
    Vx = [env.Vinf(t) for ielem = 1:n]
    Vy = [-history[1].elements[ielem].V[2] for ielem = 1:n]


    ### Store old solution 
    Fz_new = deepcopy(Fz)
    Fz_old = deepcopy(Fz)

    Fy_new = deepcopy(Fy)
    Fy_old = deepcopy(Fy)

    deflection_new = deepcopy(def_y)
    deflection_old = deepcopy(def_y)

    # @show Fz_new
    # @show deflection_new

    ### Create solution flags
    converged_flag = false
    def_flag = false
    y_flag = false
    z_flag = false
    iter = 1

    ##### Iterate to solve
    while !converged_flag

        ### Iterate through the aerodynamic stations and solve
        ## Create Section
        r_displaced = rvec #.+ def_x #Todo: Why is there deflection in the positive x? #I'm not including the displacement in the r direction, because that'll mess with the tip and hub corrections. -> 
        twist_displaced = twistvec .- def_thetax
        sects = CCBlade.Section.(r_displaced, chordvec, twist_displaced, airfoils)  

        ## Create Operating Point
        # ops = CCBlade.windturbine_op.(U, env.RS(t), pitch, r_displaced, precone, yaw, tilt, azimuth, hubHt, bemmodel.shearexp, env.rho)
        for i = 1:n #Iterate throught the nodes and update the operating points
            operatingpoints[i] = CCBlade.OperatingPoint(Vx[i], Vy[i], env.rho, pitch, env.mu, env.a)
        end
        # @show operatingpoints[end].Vx operatingpoints[end].Vy

        ## Run CCBlade
        outs = CCBlade.solve.(Ref(rotor), sects, operatingpoints)

        ## Extract CCBlade Loads
        Fz_new .= -outs.Np
        Fy_new .= outs.Tp

        outshistory[1] = outs

        # @show Fz_new

        ### Solve GXBeam
        ## Update GXBeam Loads
        for i = 1:n #Iterate through the elements and apply the distributed load at every element. 
            f_follower = SVector(0.0, Fy_new[i]/elements[i].L, Fz_new[i]/elements[i].L) #Dividing by the length of the element so the force is distributed across the element. 
            distributed_load[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
        end

        ## Run GXBeam
        syshistory[1], converged = steady_state_analysis!(system, assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_load, linear = false, angular_velocity = Omega, gravity=grav)

        history[1] = AssemblyState(syshistory[1], assembly; prescribed_conditions = prescribed_conditions)

        ## Obtain the deflections
        def_x .= [history[1].elements[ielem].u[1] for ielem = 1:n]
        def_y .= [history[1].elements[ielem].u[2] for ielem = 1:n]
        def_z .= [history[1].elements[ielem].u[3] for ielem = 1:n]
        def_thetax .= [history[1].elements[ielem].theta[1] for ielem = 1:n]

        ## Update the velocities
        for i = 1:n
            V_elem = history[1].elements[i].V # Element linear velocity
            Vx[i] = env.Vinf(t) + V_elem[3] #Freestream velocity  
            Vy[i] = - V_elem[2]
        end


        ### Update solution
        deflection_new .= def_y


        ### Calculate residuals
        resids[iter, 1] = maximum(abs.(deflection_old.-deflection_new))
        resids[iter, 2] = maximum(abs.(Fy_old.-Fy_new))
        resids[iter, 3] = maximum(abs.(Fz_old.-Fz_new))

        ### Check if the solution has converged 
        if resids[iter,1]<=tolerance
            if verbose*!def_flag
                println("Y deflection converged after $iter iterations.")
            end
            def_flag = true
        else
            def_flag = false
        end

        if resids[iter,2]<=tolerance
            if verbose*!y_flag
                println("Y loadings converged after $iter iterations.")
            end
            y_flag = true
        else
            y_flag = false
        end

        if resids[iter,3]<=tolerance
            if verbose*!z_flag
                println("Z loadings converged after $iter iterations.")
            end
            z_flag = true
        else
            z_flag = false
        end

        if def_flag*y_flag*z_flag
            converged_flag = true
            break
        end

        if iter>=maxiterations
            break
        end
        iter += 1


        ### Update old vectors
        deflection_old .= deflection_new
        Fy_old .= Fy_new
        Fz_old .= Fz_new
    end
    # iter -= 1

    return outshistory[1], history[1], syshistory[1], assembly, prescribed_conditions, converged, iter, resids[1:iter,:] #Todo. Is the state and system that are getting passed out the state and system that were most recently calculated... or the one before the loop. I'm guessing it's the one before the loop. ... but then how/why does my dynamic solution look like it is converging to the fixed-point solution? -> I know that it's passing out the correct solution now. I don't think it was doing it incorrectly in the first place. 
end