#=
Assuming that the aerodynamics respond quasi-statically.... which isn't true, but should give me some insight to solving the coupled aerostructural system. 

#### States
The states of this system are just the residual equations from CCBlade and then the GXBeam states... So

X = [x_bem; ... ; x_gxbeam]

where x_bem is just the BEM residual for a given radial node. 

=#


function create_bemgxbeamfun(radiusvec, chordvec, twistvec, rhub, rtip, bem::BEM, blade::Blade, env::Environment; b=0.01, g=9.817)
    ### BEM preamble 
    n = length(blade.airfoils)
    
    rotor = CCBlade.Rotor(rhub, rtip, 1.0, precone=0.0, turbine=bem.turbine)

    airfoils = Array{CCBlade.AFType}(undef,n)
    sections = Array{CCBlade.Section}(undef,n)




    ### Structural preamble - Initialize structs and vectors that don't need to be recreated every iteration. 
    ## Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

    ## Create point masses
    point_masses = Dict{Int64, GXBeam.PointMass{eltype(p)}}()

    ## Create GXBeam pass ins.
    start = 1:gxmodel.ne
    stop = 2:gxmodel.np
    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem = GXBeam.system_indices(start, stop, false; prescribed_points=1:1)


    bemgxfun! = function(outs, dx, x, p, t)
        ### Split up the residuals, state rates, states, and parameters
        ## Bem residual, state rates, states, and parameters
        outs_a = view(outs, 1:n)
        # dxa = view(dx, 1:n)
        xa = view(x, 1:n) #phi
        pa = view(p, 1:7*n) # radius, chord, twist, pitch, rhub, rtip, hubHt = p

        ## GXbeam residuals, state rates, states and parameters
        idxs = (n+1):length(x)
        outs_s = view(outs, idxs)
        dxs = view(dx, idxs)
        xs = view(x, idxs)
        ps = view(p, (7*n + 1):length(p))







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

            ## Extract the bem states
            xbem = view(xa, i)

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
            radius, chord, twist, pitch, rhub, rtip, hubHt = pas 
            phi = xbem[1]
            # airfoil = blade.airfoils[i]
            if isnan(phi)
                @show t, radius
            end

            # @show phi

            # Find velocities in BEM frame
            Vx = env.U(t) + V_elem[3] #Freestream velocity #TODO: I could include the system velocity in the structural frame, but.... I don't think it would accomplish what I want. 
            Vy = - V_elem[2] #Includes the angular velocity*radius, and all the displacements. #Todo: Should this be a negative V_elem. The value should be positive. 

            # ulocal = sqrt(Vx^2 + Vy^2)

            # airfoils[i] = CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", env.rho*ulocal*chord/env.mu, ulocal/env.a) 
            airfoils[i] = CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", 0.0, 0.0)
        
            ### Create section object
            radius_displaced = radius # + delta_elem[1] #Note: This might affect the tip correction. 
            twist_displaced = twist # - theta_elem[1]
            sections[i] = CCBlade.Section(radius_displaced, chord, twist_displaced, airfoils[i]) #TODO: Check how radius affects the Section function in CCBlade. 

            # @show -Omega_elem[3] #This is positive like I want. 

            ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
            # op = CCBlade.windturbine_op(Vx, -Omega_elem[3], pitch, radius, 0.0, 0.0, 0.0, 0.0, hubHt, bem.shearexp, env.rho, env.mu, env.a)

            op = CCBlade.OperatingPoint(Vx, Vy, env.rho, pitch, env.mu, env.a)

            # @show typeof(op)

            # op.Vy = Vy #I want to use the Vy I prepared up above.... but I can't alter an immutable struct. So I might just have to fill the operating point struct on my own. 

            # @show t #Takes a single time step then fails. 

            # get_bem_residual!(outs_a, i, xbem[1], rotor, sections[i], op) ### I can probably directly use the ccblade residual. 
            outs_a[i], ccout = CCBlade.residual(phi, rotor, sections[i], op) #Todo: Solution is getting forced past phi=1





            ### Create loadings for structures #Todo: I need to interpolate the CCBlade loads to the GXBeam nodes. I think that's going to require me to save the loads to a vector. Or something. Because I need to access the loads, or at least the surrounding load nodes. Which, depending on the type of interpolation, I might need the entire load framework. 
         
            f = SVector(0.0, ccout.Tp, -ccout.Np) #TODO: Should I multiply by L? Does it need the section length?
            # # f = qinf*elements[i].L*f #qinf*elements[i].L*f/2  #TODO. Why divided by 2? -> Because the BEM gives the distributed loa..... wait.... the BEM gives the distributed load.... not the total load. It shouldn't be divided by 2.  
            # # m = qinf*elements[i].L*m/2
            # @show f

            #Note: Both the force calculation above and below appear to give the same output. (I randomly checked some elements of the residual vector and they were all the same.)

            # Cl = ccout.cl
            # Cd = ccout.cd
            # qinf = 0.5*env.rho*(env.U(t)^2)*chord #TODO: Which velocity should I use, the displaced (Vx) or the freestream? 
            # f = SVector(0.0, Cd*cos(phi) - Cl*sin(phi), -(Cl*cos(phi) + Cd*sin(phi))) #Todo. I'm not sure that this is correct. 
            # f = qinf*elements[i].L*f #qinf*elements[i].L*f/2  #TODO. Why divided by 2? -> Because the BEM gives the distributed loa..... wait.... the BEM gives the distributed load.... not the total load. It shouldn't be divided by 2.  
            # m = qinf*elements[i].L*m/2


            # @show typeof(b)
            # Calculate the damping force
            f_follower = SVector(-b*V_elem[1]/elements[i].L, -b*V_elem[2]/elements[i].L, -b*V_elem[3]/elements[i].L) #Fd = -b.*ue

            # @show typeof(f) typeof(m)
            # @show typeof(f_follower) typeof(m_follower)
            
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
    end #End bemgxfun!

    return SciMLBase.DAEFunction{true, true}(bemgxfun!)

end

function create_explicitbemgxfun(rhub, rtip, bem::BEM, blade::Blade, env::Environment; b=0.01, g=9.817, damping=false)
        ### BEM preamble 
        n = length(blade.airfoils)
    
        rotor = CCBlade.Rotor(rhub, rtip, 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
    
        airfoils = Array{CCBlade.AFType}(undef,n)
        sections = Array{CCBlade.Section}(undef,n)
    
    
    
    
        ### Structural preamble - Initialize structs and vectors that don't need to be recreated every iteration. 
        ## Create prescribed conditions
        prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    
        ## Create point masses
        point_masses = Dict{Int64, GXBeam.PointMass{eltype(p)}}()
    
        ## Create GXBeam pass ins.
        start = 1:gxmodel.ne
        stop = 2:gxmodel.np
        N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem = GXBeam.system_indices(start, stop, false; prescribed_points=1:1)
    
    
        bemgxfun! = function(outs, dx, x, p, t)
            ### Split up the residuals, state rates, states, and parameters
            ## Bem residual, state rates, states, and parameters
            # outs_a = view(outs, 1:n)
            # dxa = view(dx, 1:n)
            # xa = view(x, 1:n)
            pa = view(p, 1:7*n) # radius, chord, twist, pitch, rhub, rtip, hubHt = p
    
            ## GXbeam residuals, state rates, states and parameters
            # idxs = (n+1):length(x)
            outs_s = view(outs, :)
            dxs = view(dx, :)
            xs = view(x, :)
            ps = view(p, (7*n + 1):length(p))
    
    
    
    
    
    
    
            ## Structural data
            # Create assembly
            assembly = create_gxbeam_assembly(gxmodel, ps, start, stop) 
    
            # Create distributed load 
            elements = view(assembly.elements, :) 
    
            distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(x)}}()
    
            m = @SVector zeros(3)
            m_follower = @SVector zeros(3)
    
    
            # @show eltype(x[1])
            # @show isa(x[1], ForwardDiff.Dual)
            # @show typeof(x[1])
            phivec = Array{eltype(x)}(undef, n)
    
    
            ############ Iterate across the sections. Update the aero states and create the loads for GXBeam. 
            for i = 1:n 
                ## Obtain the parameters for this aerodynamic section
                pas = view(pa, 1+7*(i-1):7+7*(i-1)) 
    
                ## Extract the bem states
                # xbem = view(xa, i)
    
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
                V_elem = SVector(xss[13], xss[14], xss[15]) #Todo: Should this be scaled at all? 
                # @show V_elem
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
                radius, chord, twist, pitch, rhub, rtip, hubHt = pas 
                # phi = xbem[1]
                # # airfoil = blade.airfoils[i]
                # if isnan(phi)
                #     @show t, radius
                # end
    
                # Find velocities in BEM frame
                Vx = env.U(t) + V_elem[3] #Freestream velocity #TODO: I could include the system velocity in the structural frame, but.... I don't think it would accomplish what I want. 
                Vy = - V_elem[2] #Includes the angular velocity*radius, and all the displacements.
    
                # ulocal = sqrt(Vx^2 + Vy^2)
    
                # airfoils[i] = CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", env.rho*ulocal*chord/env.mu, ulocal/env.a) 
                airfoils[i] = CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", 0.0, 0.0)
            
                ### Create section object
                twist_displaced = twist - theta_elem[1]
                sections[i] = CCBlade.Section(radius, chord, twist_displaced, airfoils[i])
    
                # @show -Omega_elem[3] #This is positive like I want. 
    
                ### Create OperatingPoint #TODO: Make general so it can take propeller configurations. 
                # op = CCBlade.windturbine_op(Vx, -Omega_elem[3], pitch, radius, 0.0, 0.0, 0.0, 0.0, hubHt, bem.shearexp, env.rho, env.mu, env.a)
    
                op = CCBlade.OperatingPoint(Vx, Vy, env.rho, pitch, env.mu, env.a)
    
                # @show typeof(op)
    
                # op.Vy = Vy #I want to use the Vy I prepared up above.... but I can't alter an immutable struct. So I might just have to fill the operating point struct on my own. 
    
                # @show t #Takes a single time step then fails. 
    
                # get_bem_residual!(outs_a, i, xbem[1], rotor, sections[i], op) ### I can probably directly use the ccblade residual. 
                # outs_a[i], ccout = CCBlade.residual(phi, rotor, sections[i], op) #Todo: Solution is getting forced past phi=1
                ccout = CCBlade.solve(rotor, sections[i], op)

                phivec[i] = ccout.phi
    
    
    
    
                ### Create loadings for structures
             
                f = SVector(0.0, ccout.Tp, -ccout.Np) #TODO: Should I multiply by L? Does it need the section length?
                # # f = qinf*elements[i].L*f #qinf*elements[i].L*f/2  #TODO. Why divided by 2? -> Because the BEM gives the distributed loa..... wait.... the BEM gives the distributed load.... not the total load. It shouldn't be divided by 2.  
                # # m = qinf*elements[i].L*m/2
                # @show f
    
                #Note: Both the force calculation above and below appear to give the same output. (I randomly checked some elements of the residual vector and they were all the same.)
    
                # Cl = ccout.cl
                # Cd = ccout.cd
                # qinf = 0.5*env.rho*(env.U(t)^2)*chord #TODO: Which velocity should I use, the displaced (Vx) or the freestream? 
                # f = SVector(0.0, Cd*cos(phi) - Cl*sin(phi), -(Cl*cos(phi) + Cd*sin(phi))) #Todo. I'm not sure that this is correct. 
                # f = qinf*elements[i].L*f #qinf*elements[i].L*f/2  #TODO. Why divided by 2? -> Because the BEM gives the distributed loa..... wait.... the BEM gives the distributed load.... not the total load. It shouldn't be divided by 2.  
                # m = qinf*elements[i].L*m/2
    
    
                # @show typeof(b)
                # Calculate the damping force
                if damping #Todo: How do I get the velocity of relative motion? (Not including the velocity from the blade moving through the environment) -> I think that the only velocity that should be affected by environmental variables would be the Y velocity. 
                    radius_displaced = radius + delta_elem[1] #Get the displaced radius so that I can extract the relative velocity from the system velocity. 
                    Vrel = SVector(V_elem[1], V_elem[2] + env.Omega(t)*radius_displaced, V_elem[3]) #TODO: There is a lot of radial velocity of the outer elements. Is this from the blade flipping around? Or is this actual radial displacement? 
                    # @show Vrel
                    f_follower = SVector(-b*Vrel[1]/elements[i].L, -b*Vrel[2]/elements[i].L, -b*Vrel[3]/elements[i].L) #Fd = -b.*ue
                else
                    f_follower = SVector(0.0, 0.0, 0.0)
                end
    
                # @show typeof(f) typeof(m)
                # @show typeof(f_follower) typeof(m_follower)
                
                # Store the loads. 
                distributed_loads[i] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
            end
    
            ### Save phi
            saveextra(t, "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/Rotors/testing/savingfile.csv", phivec...) #Todo: Change to the saving callback. 

            
    
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
        end #End bemgxfun!
    
        return SciMLBase.DAEFunction{true, true}(bemgxfun!)
end

function initialize_bemgx_states(bemmodel, gxmodel, env, blade, p; maxiter=100, verbose=false, tolerance=1e-6)

    n = gxmodel.ne

    outs, state, system, assembly, prescribed_conditions, converged, iters, resids = fixedpoint(bemmodel, gxmodel, env, blade, p; maxiterations = maxiter, verbose = verbose, tolerance= tolerance)

    # x0_riso = zeros(4*n) #I bet there is a way that I could calculate these. 
    # @show state.elements[end].V

    # @show outs.Tp

    x0_bem = outs.phi

    if !converged
        @show converged
    end

    x0_gxbeam = convert_assemblystate(state, assembly)

    x0 = vcat(x0_bem, x0_gxbeam)

    dx0_bem = zero(x0_bem)
    dx0_gxbeam = zero(x0_gxbeam)

    dx0 = vcat(dx0_bem, dx0_gxbeam)

    return x0, dx0
end

function extractBEMloads(filename, states, bem, blade, env, p)

    mat = readextra(filename)
    
    tvec = mat[:,1]
    data = mat[:,2:end]

    nt, nn = size(data)

    pa = view(p, 1:7*n) # radius, chord, twist, pitch, rhub, rtip, hubHt = p #Extract aerodynamic parameters
    ps = view(p, (7*n + 1):length(p)) #Extract structural parameters

    rhub = pa[5]
    rtip = pa[6]

    rotor = CCBlade.Rotor(rhub, rtip, 1.0, precone=0.0, turbine=bem.turbine, tip=nothing)
    airfoils = Array{CCBlade.AFType}(undef,n)
    sections = Array{CCBlade.Section}(undef,n)

    loadsT = Array{eltype(states)}(undef, size(data)...)
    loadsN = Array{eltype(states)}(undef, size(data)...)


    for j = 1:nt
        xs = states[j,:]
        t = tvec[j]
        for i = 1:nn
            ## Extract values out of the parameters vector for this radial position. 
            pas = view(pa, 1+7*(i-1):7+7*(i-1)) 
            radius, chord, twist, pitch, _, _, _ = pas

            ## Extract the GXBeam element states
            elem_idx = 6 + 18*(i-1)
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

            Vx = env.U(t) + V_elem[3] #Freestream velocity 
            Vy = - V_elem[2]
    
    
            airfoils[i] = CCBlade.AlphaAF(blade.airfoils[i].polar[:,1], blade.airfoils[i].polar[:,2], blade.airfoils[i].polar[:,3], "", 0.0, 0.0)
                
            ### Create section object
            radius_displaced = radius #+ delta_elem[1]
            twist_displaced = twist - theta_elem[1]
            sections[i] = CCBlade.Section(radius_displaced, chord, twist_displaced, airfoils[i])

    
            ### Create OperatingPoint #TODO: Make general so it can take propeller configurations.
            op = CCBlade.OperatingPoint(Vx, Vy, env.rho, pitch, env.mu, env.a)
    
            ccout = CCBlade.solve(rotor, sections[i], op)

            loadsT[j,i] = ccout.Tp
            loadsN[j,i] = ccout.Np
        end
    end

    return loadsT, loadsN
end