

#=
The states for this is just the aero states (bem-riso) concatonated with the gxbeam states. Really, all of the inputs that we're going to pass to Differential Equations are just the aero and the structural inputs concatonated together. 

X_as = [x_aero, x_gxbeam]

Y_as = [y_aero, y_gxbeam]

P_as = [p_aero, p_gxbeam]

You now that I'm thinking about it... there might be some repeated parameters that we can easily get rid of. We'll see. 


=#
#TODO: Break this up into a bunch of smaller functions that can be reused. (It will need to be broken up in the composing sub models.)


function create_aerostructuralfun(riso::Riso, bem::BEM, gxmodel::gxbeam, blade::Blade, env::Environment; b=0.01, g=9.817)

    ### Aero preamble - Initialize structs and vectors that don't need to be recreated every iteration. 
    na = length(blade.airfoils) #Note: I think I'm going to force the number of aerodynamic and structural states to be the same for now. 

    rarray = Array{CCBlade.Rotor}(undef, 1)
    afarray = Array{CCBlade.AFType}(undef,1)
    secarray = Array{CCBlade.Section}(undef,1)
    oparray = Array{CCBlade.OperatingPoint}(undef, 1)

    

    ### Structural preamble - Initialize structs and vectors that don't need to be recreated every iteration. 
    ## Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

    ## Create point masses
    point_masses = Dict{Int64, GXBeam.PointMass{eltype(p)}}()
    #Dict(1 => PointMass(0.0, SVector(0.0, 0.0, 0.0), @SMatrix zeros(3,3)))

    ## Create GXBeam pass ins.
    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1
    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem = GXBeam.system_indices(start, stop, false; prescribed_points=1:1)





    ######## Create the function to pass to DifferentialEquations 
    asfun! = function(outs, dx, x, p, t)

        ### Split up the states, state rates, and parameters
        ## Aero states 
        idxa = 1:5*na
        outs_a = view(outs, idxa)
        dxa = view(dx, idxa)
        xa = view(x, idxa)
        pa = view(p, 1:7*na)

        #Question: I just realized that the way that I'm doing this doesn't have the y... so will that affect how derivatives are calculated? 

        ## Structural states
        idxs = (5*na +1):length(x)
        outs_s = view(outs, idxs)
        dxs = view(dx, idxs)
        xs = view(x, idxs)
        ps = view(p, (7*na + 1):length(p))





        ### Prepare structs, vectors, etc for this iteration. 
        ## Aerodynamic data
        alphavec = eltype(x).(-pi:.1:pi) #Todo. How to initialize this to be of the same type as Cl? Additionally, I don't know if I need to initialize all of the elements. 
        clvec = ones(eltype(x), length(alphavec)) #TODO: This might be able to be moved outside of the function. I don't think it is required.... but dual numbers might pass through it. So I might need to allow it's type to change. 
        cdvec = ones(eltype(x), length(alphavec))
        

        # Create Rotor
        if bem.tipcorrection  
            rarray[1] = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine)
        else
            rarray[1] = CCBlade.Rotor(p[5], p[6], 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
        end

        ## Structural data
        # Create assembly
        assembly = create_gxbeam_assembly(gxmodel, ps, start, stop) 

        # Create distributed load 
        elements = view(assembly.elements, :) 

        distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(x)}}()

        m = @SVector zeros(3)
        m_follower = @SVector zeros(3)






        ############ Iterate across the sections. Update the aero states and create the loads for GXBeam. 
        for i = 1:n 
            ## Obtain the parameters for this aerodynamic section
            pas = view(pa, 1+7*(i-1):7+7*(i-1)) 

            ## Extract the riso states
            riso_idx = 4*(i-1)
            dxs_riso = view(dxa, 1+riso_idx:riso_idx+4) 
            xs_riso = view(xa, 1+riso_idx:riso_idx+4) 
            ps_riso = pas[2]

            ## Extract the bem states
            bem_idx = 4*n + i
            xbem = view(xa, bem_idx)

            ## Extract the GXBeam element states
            elem_idx = 6 + 18*(i-1)
            dxss = view(dxs, (elem_idx+1):(elem_idx + 18))
            xss = view(xs, (elem_idx+1):(elem_idx + 18))





            ## Get the structural displacements
            # Linear displacement
            delta_elem = SVector(xss[1], xss[2], xss[3])

            # angular displacement
            theta_elem = SVector(xss[4], xss[5], xss[6]) 
            # dtheta_elem = SVector(dxss[3], dxss[4], dxss[5])

            # Linear velocity
            V_elem = SVector(xss[13], xss[14], xss[15])
            # dV_elem = SVector(dxss[12], dxss[13], dxss[14])

            # angular velocity
            Omega_elem = SVector(xss[16], xss[17], xss[18])
            # dOmega_elem = SVector(dus[15], dus[16], dus[17])

            # convert rotation parameter to Wiener-Milenkovic parameters
            scaling = GXBeam.rotation_parameter_scaling(theta_elem)
            theta_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)
            # dtheta_elem *= scaling # angular displacement rate (Wiener-Milenkovic parameters)




            ## Calculate the time varying inputs for the aerodynamic models.
            # Extract inputs
            radius, chord, twist, pitch, _, _, _ = pas
            phi = xbem[1]
            airfoil = blade.airfoils[i]

            # Find velocities in BEM frame
            Vx = env.U(t) + V_elem[3] #Freestream velocity #TODO: I could include the system velocity in the structural frame, but.... I don't think it would accomplish what I want. 
            displaced_radius = radius + delta_elem[3]
            displaced_omega = - Omega_elem[3]
            Vy = displaced_radius*displaced_omega - V_elem[2]
            #TODO: Double check these velocities with someone. (DG or Ning)

            # @show displaced_omega #This is a positive value... which is what I wanted. 


            #=
            - Including the positive V_elem[3] in Vx should account for any preconing deflection. It is positive because if the beam moves in the positve Z direction (structural frame) that would mean that the airfoil is seeing the wind speed plus the speed of the moving beam. 
            - I've calculated the displaced radius. GXBeam will tell me how far my element has moving in the structural x direction, which elongates my beam. If I have a longer beam, then each section will see higher windspeeds. 
            - For the displaced rotational velocity (displaced_omega), I've used the negative of the structural z angular velocity. GXBeam's angular velocity includes displaced velocities and the angular velocity of the entire system. I really want to use the absolute value, but I think it should be the negative, because for an upwind turbine, our structural angular velocity should be about the negative z. 
            - For Vy, I use the displaced radius and displaced omega to get the speed that the blade sees due to spinning. Then if the blade has some positive deformation velocity, that should slow down the velocity that the airfoil sees (BEM) frame.
            =#


            Vxdot = env.Udot(t) #TODO: These need updating to include structural deflections, but I'm not sure how important it'd be. 
            Vydot = env.Omegadot(t)*radius 

            # Find velocities in dynamic stall frame. (Not accounting for induced velocities in Dynamic stall model.)
            u = sqrt(Vy^2 + Vx^2)  
            v = 0 

            udot = sqrt(Vxdot^2 + Vydot^2)
            vdot = 0.0

            alpha = -((twist + pitch) - phi)
            alphadot = 0.0 

            ys_riso = SVector(u, udot, v, vdot, alpha, alphadot)

            # Find force coefficients
            Cl, Cd = riso_coefs(xs_riso, ys_riso, chord, airfoil)

            ys_bem = SVector(Cl, Cd, Vx, Vy) #TODO: I don't think that I need to store these Y values in vectors. I should be able to just used them where I need them. 

            outs_a[1+riso_idx:riso_idx+4] = riso_residual(dxs_riso, xs_riso, ys_riso, ps_riso, t, blade.airfoils[i])

            U = sqrt(ys_bem[3]^2 + ys_bem[4]^2)
            clvec[:] .= ys_bem[1]
            cdvec[:] .= ys_bem[2]
             
            afarray[1] = CCBlade.AlphaAF(alphavec, clvec, cdvec, "", env.rho*U*pas[2]/env.mu, U/env.a) #Note: We definitely want to leave the airfoils in the replacing array, they take a crap ton of space. -> Maybe I can overload the structs to use a function that doesn't take allocations. 
            
            ### Create section object
            section = CCBlade.Section(pas[1], pas[2], pas[3], afarray[1])

            ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
            oparray[1] = CCBlade.windturbine_op(ys_bem[3], ys_bem[4]/pas[1], pitch, pas[1], 0.0, 0.0, 0.0, 0.0, pas[7], bem.shearexp, env.rho, env.mu, env.a)

            get_bem_residual!(outs_a, bem_idx, xbem[1], rarray[1], section, oparray[1])






            ### Create loadings for structures
            qinf = 0.5*env.rho*(env.U(t)^2)*chord #TODO: Which velocity should I use, the displaced (Vx) or the freestream? 
            f = SVector(0.0, Cd*cos(phi) - Cl*sin(phi), -(Cl*cos(phi) + Cd*sin(phi))) #Todo. I'm not sure that this is correct. 
            f = qinf*elements[i].L*f #qinf*elements[i].L*f/2  #TODO. Why divided by 2? -> Because the BEM gives the distributed loa..... wait.... the BEM gives the distributed load.... not the total load. It shouldn't be divided by 2.  
            m = qinf*elements[i].L*m/2


            # Calculate the damping force
            f_follower = SVector(-b*V_elem[1]/elements[i].L, -b*V_elem[2]/elements[i].L, -b*V_elem[3]/elements[i].L) #Fd = -b.*ue

            #TODO: Why am I dividing by L? I should figure out what this should be. -> Oh, my thought process is that -b*u should be the total force of an element... but I think that this force is integrated across the element, so I should divide -b*u by the length of the element.

            # I apply damping as a follower force because A) it would be, and B) I might as well, because I need to define a follower force and it might as well be the damping. (It wouldn't make a difference if it was a follower or not, because it will be updated from time step to time step, which is when a non-following force would change.)
            # @show typeof(V_elem)
            # @show typeof(V_elem[1])
            # @show b
            # @show typeof(elements[i].L)
            # @show typeof(f), typeof(m)
            # @show typeof(f_follower), typeof(m_follower)
            # @show typeof(GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower))
            
            # Store the loads. 
            distributed_loads[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
        end



        ### System velocities and accelerations
        ## create the origin
        x0 = @SVector zeros(3) #System origin
        v0 = @SVector zeros(3) #System linear velocity 
        omega0 = SVector(0.0, 0.0, -env.Omega(t)) #System angular velocity
        a0 = @SVector zeros(3) #System linear acceleration
        alpha0 = SVector(0.0, 0.0, -env.Omegadot(t)) #System angular acceleration

        ## Create gravity vector 
        rotated_angle = env.Omega(t)*t 
        gx = g*sin(rotated_angle)
        gy = -g*cos(rotated_angle)
        gz = 0.0
        gvec = SVector(gx, gy, gz) 

        force_scaling = GXBeam.default_force_scaling(assembly) 

        ## Update GXBeam states
        GXBeam.dynamic_system_residual!(outs_s, dxs, xs, assembly, prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, omega0, a0, alpha0)

        return outs
    end #End asfun!

    return SciMLBase.DAEFunction{true, true}(asfun!)
end

function initialize_riso_states(aoa, dsmodel, blade, env, p; tspan=(0,100.0))

    ### Read in data
    n = length(blade.airfoils)

    ### Initialize states
    x0 = zeros(4*n)
    dx0 = zeros(4*n)

    fun = create_risofun(aoa, blade, env, 0.0, 0.0) #The frequency and amplitude of oscillation is zero. 

    diffvars = differentialvars(dsmodel, n)

    prob = DifferentialEquations.DAEProblem(fun, dx0, x0, tspan, p, differential_vars = diffvars)

    sol = DifferentialEquations.solve(prob, DABDF2())

    x = Array(sol)'
    return x[end,:]
end

function initialize_aerostructural_states(bemmodel, gxmodel, env, blade, p; maxiter=100, verbose=false, tolerance=1e-6)

    n = gxmodel.n

    outs, state, system, assembly, prescribed_conditions, converged, iters, resids = fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations = maxiter, verbose = verbose, tolerance= tolerance)

    # x0_riso = zeros(4*n) #I bet there is a way that I could calculate these. 


    x0_bem = outs.phi

    pa = view(p, 1:7*n)
    chordidx = 2:7:7*n
    chord = view(pa, chordidx)
    twistidx = 3:7:7*n
    twist = view(pa, twistidx)
    pitchidx = 4:7:7*n
    pitch = view(pa, pitchidx)
    aoa = -((twist .+ pitch) .- x0_bem)
    # @show aoa
    x0_riso = initialize_riso_states(aoa, Riso(), blade, env, chord) #Todo: Getting NaN in my states



    x0_gxbeam = convert_assemblystate(state)

    x0 = vcat(x0_riso, x0_bem, x0_gxbeam)

    dx0_riso = zero(x0_riso)
    dx0_bem = zero(x0_bem)
    dx0_gxbeam = zero(x0_gxbeam)

    dx0 = vcat(dx0_riso, dx0_bem, dx0_gxbeam)

    return x0, dx0
end