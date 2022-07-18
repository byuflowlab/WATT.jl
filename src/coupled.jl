#=
This is an attempt at not tightly coupling the models. The tight coupling is proving difficult to solve, especially with the damping as part of the structural model. So we're just going to scale back. Not tightly coupled. Just passing loads and deflections every time step. Maybe in the future I'll allow for some interior iterations (in a fixed point solution). 





##### Work Notes
- I just had the thought. I should go simpler than this. I should do a solution of just BEM to GXBeam across time. Simple one-way coupling. Then do a two-way coupling. Then add in the dynamic stall model. Then add in the possibility of having multiple interior iterations.
    - I've done the simple one-way coupling. I haven't validated it, but it looks good. 
    - I got the two-way coupling done. I'm not sure if it's 100% correct, but it looks good. It doesn't seem to oscillate about the steady state solution... but I'm not crazy worried about that. It might be due to the fact that the blade is massive. I'll see what happens once I add dynamic stall and damping. 

=#

abstract type Coupling end

struct OneWay <: Coupling
end

struct TwoWay <: Coupling
end

struct ThreeWay <: Coupling
end


struct RisoState{TF}
    x::Array{TF, 1}
    u::TF
    c::TF
    airfoil::Airfoil
end

function (state::RisoState)(alpha, Re , Mach) #TODO: I might be able to combine this with the function below.... which might be good. 
    y = [state.u, 0.0, 0.0, 0.0, alpha, 0.0]
    # return riso_coefs(state.x, y, state.c, state.airfoil) #TODO: I don't know if the fact that this function returns a Static Vector will throw a wrench in the works, but we'll just have to see... won't we? 
    cf = riso_coefs(state.x, y, state.c, state.airfoil)
    # if cf[1]<0
    #     @show cf
    # end
    return cf
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
        # @show length(dx)
        for i = 1:n
            idx = 4*(i-1)
            xs = x[1+idx:idx+4]

            u = Wvec[i] #sqrt(Vxvec[i]^2 + Vyvec[i]^2) #Condense the inflow velocity to a single value.
            v = 0 

            udot = Wdotvec[i] #sqrt(Vxdotvec[i]^2 + Vydotvec[i]^2)  
            vdot = 0.0

            theta = -((twistvec[i] + pitch) - phivec[i]) #((twistvec[i] + pitch) - phivec[i]) #Calculate the angle of attack that the airfoil see #Note: In CCBlade's out struct, it has the alpha of the airfoil section... and for a wind turbine... all the values are positive. All of the inflow values are also positive.... So I don't know... I'll try not negating this whole thing and see what happens. -> Except when you check the source code, it evaluates the negative of the angle of attack. -> And it stores the negative in the angle of attack. Which means.... that the aoa was negative in the first place... which tells me that I've got something backwards somewhere. What I'm using (with the negative out front) matches what I would expect using the diagram for wind turbine operation from CCBlade. 
            thetadot = 0.0 #Assuming that the airfoil isn't in the act of turning???? TODO: What are my other options?

            # if theta <0 #Theta is pretty much always negative. Which tells me that whatever this angle of attack is.... it is likely wrong. Should my angle of attack just be phi? 
            #     @show t, i, theta
            # end
    
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
function update_p_a!(p_a, deformed_twist, phivec, Wvec, Wvecdot; n=length(deformed_twist))
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
function interpolate_velocity(ip, assembly, state) #Note: This is probably a very poor approximation. But it'll have to do for now. Maybe there is a better way to get the velocity of the points.  
    ne = length(assembly.elements)
    np = ne + 1

    V = state.elements[ip.pair[1]].V  #Get the velocity of the element 
    Omega = state.elements[ip.pair[1]].Omega #Question: If I'm doing what I'm currently doing.... why don't I just calculate the distance from the element node to the aerodynamic node? 
    p1 = assembly.points[ip.pair[1]]
    p2 = assembly.points[ip.pair[2]]
    e = assembly.elements[ip.pair[1]].x

    r1 = p1 - e #The position vector from the element point to the starting point
    r2 = p2 - e #The position vector from the element point to the stopping point

    # ra = assembly.elements[ip.pair[1]].L*ip.percent + assembly.points[ip.pair[1]][1]
    
    v1 = V + cross(Omega, r1)
    v2 = V + cross(Omega, r2)

    # if ip.pair[2]==np
    #     println("")
    #     @show r1, r2, Omega
    #     @show cross(Omega, r1), cross(Omega, r2)
    # end

    # if ip.pair[2] == 18
    #     @show V #V is what is making the velocity go super wonky. 
    # end

    return (1-ip.percent)*v1 + ip.percent*v2
end

function interpolate_angle(ip, assembly, state)
    thetagx = [state.points[i].theta[1] for i = 1:length(assembly.points)]
    return (1-ip.percent)*thetagx[ip.pair[1]] + ip.percent*thetagx[ip.pair[2]]
end


"""
update_aero_inputs!(coupling::OneWay, sections, operatingpoints, Vxvec, Vyvec, env, dsstates, rvec, aeroV, twist_displaced, t) #TODO: I need to rethink the name of this function. the _blah_ syntax underlines blah in a docstring... which isn't useful in a title. 

    Updates the inputs for the aerodynamic models for a one-way coupling (BEM to GXBeam). Specifically it updates the sections, operatingpoints, Vxvec, Vyvec, and Wvec. 

### Inputs
- coupling::OneWay - Coupling struct. Used for multiple dispatch. 
- sections::Array{CCBlade.Section, 1} - Vector of CCBlade Sections. Length is equal to the number of aerodynamic nodes. 
- operatingpoints::Array{CCBlade.OperatingPoint, 1} - Vector of CCBlade OperatingPoints. Length is equal to the number of aerodynamic nodes.
- Vxvec:: - Vector of velocities in freestream direction
- Vyvec:: - Vector of velocities due to rotation. 
- Wvec::Array - Vector of total inflow velocity.
- env::Environment
- rvec::Array{TF, 1}
- aeroV:: - Vector of velocities due to structural deflection. Unused in this version of the function. 
- twist_displaced:: - Vector of twist at every aerodynamic node. Includes structural displacements. 
- t::TF - current simulation time.
- pitch::TF - current blade pitch

### Outputs
- inplace function... so no outputs. but sections, operatingpoints, Vxvec, and Vyvec are all updated. 

### Notes
- This version of the function is for a one-way coupling, meaning that none of the structural outputs are actually input back into the aerodynamic model every time step. This means that for this function, there will be several inputs that are unused. They are just there for ease of multiple dispatch. 
- I decided not to update the rvec for ease of coupling between the aerodynamic mesh and the structural mesh. I just use the original location of the aerodynamic mesh to create a fit, then use the origianl location of the structural mesh to determine where the load should be applied. -> TODO: A lot of content from this note really belongs in a header comment. 
"""
function update_aero_inputs!(coupling::OneWay, sections, operatingpoints, Vxvec, Vyvec, Wvec, Wdotvec, env, blade, dsstates, rvec, aeroV, twist_displaced, t, pitch)
    na = length(sections)

    for i = 1:na
        Vxvec[i] = env.Vinf(t) 
        Vyvec[i] = env.RS(t)*rvec[i]  
        Wvec[i] = sqrt(Vxvec[i]^2 + Vyvec[i]^2)
        operatingpoints[i] = CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a)
    end
end

function update_aero_inputs!(coupling::TwoWay, sections, operatingpoints, Vxvec, Vyvec, Wvec, Wdotvec, env, blade, dsstates, rvec, aeroV, twist_displaced, t, pitch) 
    na = length(sections) 

    # Wdotvec .= 0.0 #Note: I think this will need to change, but I'm not 100% how to get the acceleration out of the structural model.

    for i = 1:na
        sections[i] = Section(rvec[i], chordvec[i], twist_displaced[i], blade.airfoils[i])  

        Vxvec[i] = env.Vinf(t) + aeroV[i][3] 
        Vyvec[i] = aeroV[i][2] #For some reason this was negative before and isn't anymore. 
        Wvec[i] = sqrt(Vxvec[i]^2 + Vyvec[i]^2)

        operatingpoints[i] = CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a)
    end
end

function update_aero_inputs!(coupling::ThreeWay, sections, operatingpoints, Vxvec, Vyvec, Wvec, Wdotvec, env, blade, dsstates, rvec, aeroV, twist_displaced, t, pitch) #TODO: Will need to add Wvec and other Riso inputs in here. 
    na = length(sections) 

    # Wdotvec .= 0.0 #Note: I think this will need to change, but I'm not 100% how to get the acceleration out of the structural model.

    for i = 1:na
        sections[i] = Section(rvec[i], chordvec[i], twist_displaced[i], dsstates[i]) #Section(rvec[i], chordvec[i], twist_displaced[i], blade.airfoils[i])  

        Vxvec[i] = env.Vinf(t) + aeroV[i][3] 
        Vyvec[i] = aeroV[i][2] #For some reason this was negative before and isn't anymore. 
        Wvec[i] = sqrt(Vxvec[i]^2 + Vyvec[i]^2)

        operatingpoints[i] = CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a)
    end
end


"""
update_dsstates!(x, dsstates, Wvec, chordvec, blade)

    Updates the dynamic stall state object vector based on the new state vector x. 

### Inputs
- x::Array{TF, 1} - A vector of Riso states. Four states for each aerodynamic node (4*na). Arranged with the four states concatonated together [x1, x2, x3, x4] for each node, then vcat that sucker. 
- dsstates::Array{RisoState, 1} - A vector holding the RisoState objects, one for each aerodynamic node. 
- Wvec::Array{TF, 1} - A vector holding the total inflow velocity for each aerodynamic node. 
- chordvec::Array{TF, 1} - A vector holding the chord length for each aerodynamic node. 
- blade::Blade - the blade info. 
"""
function update_dsstates!(x, dsstates, Wvec, chordvec, blade)
    na = length(dsstates)
    for i = 1:na
        dsstates[i] = RisoState(x[4*(i-1)+1:4*(i-1)+4], Wvec[i], chordvec[i], blade.airfoils[i])
    end
end


"""
simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, assembly, tspan, dt; solver=RK4(), verbose=false, coupling=OneWay())

Simulates the described wind turbine across the given time span using a BEM (I am using CCBlade.) that is extended by the Riso dynamic stall method. A structural response is simulated by GXBeam. The BEM is assumed to be constant across a given time step. The BEM is influenced by the current state of the dynamic stall method. The dynamic stall and structural models are then integrated forward based on the BEM solution. The dynamic stall model is integrated using the provided solver (defaults to a Runge-Kutta 4th order) and the structural model is solved using it's internal dynamic solver. 

### Inputs
- rvec::Array{TF, 1} - A vector holding the radial distance from the center of rotation of all of the aerodynamic nodes (meters)
- chordvec::Array{TF, 1} - A vector holding the chord lengths of the blade at all of the aerodyanmic nodes (meters)
- twistvec::Array{TF, 1} - A vector holding the aerodynamic twist at every aerodynamic node (radians)
- rhub::TF - the radial location of the hub (meters) #TODO: Should be able to remove this since it's in the blade. 
- rtip::TF - the radial location of the tip (meters) # ""
- blade::Blade - a blade struct
- env::Environment
- assembly
- tspan::Tuple{TF, TF}
- dt::TF
- solver::Solver - the type of integration method used for integrating the Riso states. Defaults to RK4()
- verbose::Bool - Whether or not to print out information at every time step. 
- coupling::Coupling - a struct to declare to the model which type of coupling used. OneWay() or TwoWay(). 

### Outpus


### Notes
- The structural model does not include damping. 
- The initial condition is assumed to be the loads from the undeflected solution (to add the option to start from the fixed point iteration solution later). Those loads are then used to find the steady state deflections and inflow velocities. The inflow velocities are used to calculate the steady state states for the Riso model.
"""
function simulate(rvec, chordvec, twistvec, rhub, rtip, blade, env, assembly, tspan, dt; solver=RK4(), verbose=false, coupling=OneWay())

    t0 = tspan[1]
    tvec = tspan[1]:dt:tspan[2]
    nt = length(tvec)

    ### Create CCBlade inputs
    rotor = CCBlade.Rotor(rhub, rtip, 1.0, precone=0.0, turbine=true)
    na = length(blade.airfoils)
    sections = [Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i in 1:na]

    ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
    Vxvec = [env.Vinf(t0) for i in 1:na] #Freestream velocity
    Vyvec = [env.RS(t0)*rvec[i] for i in 1:na]
    pitch = 0.0 #TODO: I'm not sure if this pitch does something or not. And I don't know if I want to pitch yet. 
    operatingpoints = [CCBlade.OperatingPoint(Vxvec[i], Vyvec[i], env.rho, pitch, env.mu, env.a) for i in 1:na] 
    #TODO: Dr. Ning uses mu=1 and a=1 in his example.  

     
    ### Solve BEM residual
    ccout = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

    

    ### Get steady state Riso states
    Wvec = [sqrt(Vxvec[i]^2 + Vyvec[i]^2) for i = 1:na] #Total inflow velocity (assumed) -> Question: will the inflow velocity from CCBlade be different? 
    Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i])^2) for i = 1:na] #Derivative of the total inflow velocity

    risoode = createrisoode(blade) 

    p_a = vcat(chordvec, twistvec, ccout.phi, Wvec, Wdotvec, pitch) #Create p_a 

    x0 = zeros(4*na) #TODO: I think that one of these states might need to be initialized as phi.
    x0[4:4:4*na] .= 1.0 

    prob = SteadyStateProblem(risoode, x0, p = p_a)
    sol = DifferentialEquations.solve(prob)
    # @show sol.u
    xds = sol.u #Array(sol)' #TODO. I think this will work. -> The array didn't work. It was returning an array instead of a vector. 

    dsstates = [RisoState(sol[4*(i-1)+1:4*(i-1)+4], Wvec[i], chordvec[i], blade.airfoils[i]) for i = 1:na]
    # sections = [Section(rvec[i], chordvec[i], twistvec[i], dsstates[i]) for i in 1:na] #TODO: Should probably make a multiple-dispatch function for this bad boy. Or an if statement. :| 
    sections = [Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i]) for i in 1:na]

    ### Create GXBeam State & constant inputs 
    ## Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed 


    ### Extract CCBlade Loads and create a distributed load
    Fzfit = Akima(rvec, -ccout.Np )
    Fyfit = Akima(rvec, ccout.Tp)



    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]
    distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(rgx[ielem]), fz= (s) -> Fzfit(rgx[ielem])) for ielem in 1:nelem)  



    Omega = SVector(0.0, 0.0, env.RS(t0))

    system, converged = steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega) 

    gxstate = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

    
    cchistory = [ccout for i = 1:nt]
    gxhistory = [gxstate for i = 1:nt]
    # dshistory = [dsstates for i = 1:n] #Question: Do I need to save these across time? 


    interpolationpoints = create_interpolationpoints(assembly, rvec) 
    aeroV = [interpolate_velocity(interpolationpoints[j], assembly, gxstate) for j = 1:na]
    def_thetax = [interpolate_angle(interpolationpoints[j], assembly, gxstate) for j = 1:na] 

    ### Iterate through time steps
    for i = 1:nt-1
        ## Update timestep
        t = tvec[i]


        ## Update CCBlade inputs (based on GXBeam state)
        twist_displaced = twistvec .- def_thetax 
        update_aero_inputs!(coupling, sections, operatingpoints, Vxvec, Vyvec, Wvec, Wdotvec, env, blade, dsstates, rvec, aeroV, twist_displaced, t, pitch)
        # @show t, Vyvec[end]

        # Solve BEM residual
        cchistory[i+1] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

        # if i ==5 
            # @show cchistory[i+1].W #Todo: ccout.W and Wvec are very different from each other..... 
            # @show Wvec
            # @show twist_displaced #Positive. I think they are in an acceptable range. 
            # @show cchistory[i+1].phi #Positive, ranges from .1 to .1e-5
            # @show cchistory[i+1].alpha #has some negatives.. actually is mostly negative. 
            # @show cchistory[i+1].cl #half negative, half positive
        # end

        ## Solve Riso states based on inflow and deflection
        # update_p_a!(p_a, twist_displaced, ccout.phi, Wvec, Wdotvec) #TODO: Update p_a #WorkLocation: The Riso model is off somehow. -> I wonder if there is something wrong with the BEM... because if I recall... the BEM-GXBeam coupling wasn't in the same force range as the static coupling. -> I'm going to try the two-way coupling without the dynamic stall model again. Make sure that that oscillates about the solution that I expect. 
        # @show size(xds)
        # Note: I could just pass the angle of attack and inflow velocity from CCBlade to the DS model. 
        # xds[:] = solver(risoode, xds, p_a, t, dt) #TODO: I need to convert this to work in place. #Todo. This function isn't working properly. It's creating more states. -> I was passing in an array instead of a vector. 
        # update_dsstates!(xds, dsstates, Wvec, chordvec, blade) #Todo: A lot of these states are returning a negative state. 

        
        ## Update GXBeam loads
        Fzfit = Akima(rvec, -cchistory[i].Np ) #I use the loads of the previous time. We want everything to be calculated based off of the previous time step to update the next time step. Then we'll move to that step.  
        Fyfit = Akima(rvec, cchistory[i].Tp)
        
        distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(rgx[ielem]), fz= (s) -> Fzfit(rgx[ielem])) for ielem in 1:nelem)

        Omega = SVector(0.0, 0.0, env.RS(t))


        ## Solve GXBeam for time step
        system, localhistory, converged = time_domain_analysis!(system, assembly, tvec[i:i+1]; prescribed_conditions=prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega, reset_state=false, initialize=false) #TODO: Add gravity. 



        ## Extract GXBeam outputs
        gxhistory[i+1] = localhistory[end]  


        for j = 1:na
            aeroV[j] = interpolate_velocity(interpolationpoints[j], assembly, localhistory[end])
            def_thetax[j] = interpolate_angle(interpolationpoints[j], assembly, localhistory[end]) 
        end

        if verbose
            println("t=", round(t, digits=4))
        end
    end

    if verbose
        println("t=", round(tvec[end], digits=4))
    end

    return cchistory, gxhistory, tvec
end

