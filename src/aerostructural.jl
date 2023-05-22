#=


=#



function checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    if minimum(rvec)<=rhub
        @warn("A member of rvec is less than or equal to rhub. This will cause problems with the hub corrections.")
    end
    if maximum(rvec)>=rtip
        @warn("A member of rvec is greater than or equal to rtip. This will cause problems with the tip corrections.")
    end
    if pitch>(pi/4)
        @warn("The pitch value is greater than 45 degrees. Ensure that the pitch is in radians and not degrees.")
    end
    if precone>(pi/4)
        @warn("The precone value is greater than 45 degrees. Ensure that the precone is in radians and not degrees.")
    end
    if tilt>(pi/4)
        @warn("The tilt value is greater than 45 degrees. Ensure that the tilt is in radians and not degrees.")
    end
    if yaw>(pi/4)
        @warn("The yaw value is greater than 45 degrees. Ensure that the yaw is in radians and not degrees.")
    end
    if maximum(twistvec)>(pi/2)
        @warn("The maximum twist value is greater than 90 degrees. Ensure that the twist distribution is given in radians.")
    end
end

#TODO: It might be a good idea to make a version that is completely in place. (Pass in data storage and riso ode.)

#TODO: Function headers
function simulate(rotor::Rotors.Rotor, blade::Blade, env::Environment, assembly::GXBeam.Assembly, tvec; pitch=0.0, turbine::Bool=true, solver::Solver=RK4(), verbose::Bool=false, speakiter::Int=100, warnings::Bool=true, azimuth0=0.0, structural_damping::Bool=true, linear::Bool=false, g=9.81, plotbool::Bool=false, plotiter::Int=speakiter)
    #Todo: Move the DSM initialization type (the current ModelInit struct) to DSM. 
    # println("Why isn't it printing inside the function.")
    if verbose
        println("Rotors.jl initializing solution...")
    end

    # if warnings
    #     checkforwarnings(rvec, twistvec, rhub, rtip, pitch, precone, tilt, yaw)
    # end

    #TODO: I need to put something in to be able to disable the dynamic stall model (near the hub). 
    #TODO: I might put some checks in to help with problems caused by having rvec[i]~=rhub or rvec[i]~=rtip. 
    #TODO: I might put some checks in to see if the airfoils are in degrees or radians. 
    #TODO: I might put some checks in to see if twist or pitch is in degrees or radians. 
    #TODO: I might put something in to make the solution run smoother when the dynamic stall states would get screwed up (i.e. at the root cylinders, when rvec[i]=rhub or rtip). 
    #TODO: It might be a good idea to check rvec, chordvec, and twistvec to get the design variables to get the right types. 

    ### Initialization information
    na = length(blade.rR)
    nt = length(tvec)

    t0 = tvec[1]


    twistvec = blade.twist
    airfoils = blade.airfoils
    chordvec = airfoils.c

    ### Prepare data storage
    Vxvec = Array{eltype(chordvec)}(undef, na) #TODO: These might be too type specific. 
    Vyvec = Array{eltype(chordvec)}(undef, na)
    azimuth = Array{eltype(chordvec)}(undef, nt)
    # Vxvec = @SVector zeros(eltype(chordvec), na)  #Note: I guess static arrays are immutable... even though the docs say they aren't. 
    # Vyvec = @SVector zeros(eltype(chordvec), na) 
    # azimuth = @SVector zeros(eltype(chordvec), nt) 
    Cx = Array{eltype(chordvec)}(undef,(nt, na))
    Cy = Array{eltype(chordvec)}(undef,(nt, na))
    Cm = Array{eltype(chordvec)}(undef,(nt, na))

    Fx = Array{eltype(chordvec)}(undef,(nt, na))
    Fy = Array{eltype(chordvec)}(undef,(nt, na))
    Mx = Array{eltype(chordvec)}(undef,(nt, na)) #Moment about the aero Z axis is the moment about the structural X axis. 

    # sections = Array{CCBlade.Section{eltype(rvec)}, 1}(undef, na)
    # operatingpoints = Array{CCBlade.OperatingPoint{eltype(rvec)}, 1}(undef, na)
    cchistory = Array{CCBlade.Outputs{eltype(chordvec)}, 2}(undef, nt, na) 
    xcc = zeros(11)

    gxhistory = Array{GXBeam.AssemblyState{eltype(chordvec), Vector{GXBeam.PointState{eltype(chordvec)}}, Vector{GXBeam.ElementState{eltype(chordvec)}}}}(undef, nt)




    ### Initialize BEM solution
    azimuth[1] = azimuth0


    for j = 1:na
        # sections[i] = CCBlade.Section(rvec[i], chordvec[i], twistvec[i], blade.airfoils[i])
        Vxvec[j], Vyvec[j] = get_aero_velocities(rotor, blade, env, t0, j, azimuth[1])

        cchistory[1,j] = solve_BEM!(rotor, blade, env, j, Vxvec[j], Vyvec[j], pitch, xcc)
    end

    # cchistory[1,:] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)

    # @show env.RS(t0), env.U(t0)




    ### Initialize DS solution
    # ode = initializeDSmodel(dsmodel, solver) 

    # Wdotvec = SVector{na}(zeros(na)) #TODO: Typing
    Wdotvec = zeros(eltype(chordvec), na)
    alphadotvec = zeros(eltype(chordvec), na)
    # for i = 1:na
    #     Wdotvec[i] = sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2)
    # end

    # Wdotvec = [sqrt(env.Vinfdot(t0)^2 + (env.RSdot(t0)*rvec[i]*cos(precone))^2) for i in 1:na] #TODO: I probably need to update if there is precone, tilt, yaw, etc. -> Maybe I'll make a function to do this. 

    # ode, xds, p_ds = initializeDSmodel(dsmodel, dsmodelinit, solver, turbine, nt, na, tvec, cchistory[1, :], Wdotvec, chordvec, twistvec, pitch, env.a)
    xds, xds_idxs, p_ds = initialize_DS_model(airfoils, turbine, nt, tvec, cchistory[1, :], Wdotvec, alphadotvec, twistvec, pitch)

    # @show xds[1,end-31:end]


    # Cx[1,:], Cy[1,:], Cn[1,:], Ct[1,:], Cl[1,:], Cd[1,:], Cm[1,:] = extractloads(dsmodel, xds[1,:], cchistory[1], chordvec, twistvec, pitch, blade, env) #Todo: I might need to individual rotations for the stations. 
    #Todo. Why am I storing Cl, Cn, and Cx? 
    # extractloads!(dsmodel, xds[1,:], cchistory[1, :], chordvec, twistvec, pitch, blade, env, view(Cx, 1, :), view(Cy, 1, :), view(Cn, 1, :), view(Ct, 1, :), view(Cl, 1, :), view(Cd, 1, :), view(Cm, 1, :))

    extract_ds_loads!(airfoils, view(xds, 1, :), xds_idxs, cchistory[1, :], p_ds, view(Cx, 1, :), view(Cy, 1, :), view(Cm, 1, :)) #Todo: With Precone, and any deflections, then Cx, and Cy will no longer be in the blade root frame. 
    
    ### Dimensionalize #TODO: Could probably make a function out of this. 
    #Todo: These loads may need to be rotated from a local frame to a hub frame. 
    for j = 1:na
        u_1 = cchistory[1, j].W
        qinf_1 = 0.5*env.rho*u_1^2
        # N[1,j] = Cn[1,j]*qinf_1*chordvec[j] 
        # T[1,j] = Ct[1,j]*qinf_1*chordvec[j]
        Fx[1,j] = Cx[1,j]*qinf_1*chordvec[j] 
        Fy[1,j] = Cy[1,j]*qinf_1*chordvec[j] 
        Mx[1,j] = Cm[1,j]*qinf_1*chordvec[j]^2 #The coefficient of moment is positive about the negative Z aero axis, so we need the negative of this to move it to the structural X axis. 
    end



    ### Extract CCBlade Loads and create a distributed load 

    # nelem = length(assembly.elements)
    # rgx = [assembly.elements[i].x[1] for i in 1:nelem]

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(chordvec)}}()
    update_forces!(distributed_loads, view(Fx, 1,:), view(Fy, 1,:), view(Mx, 1,:), blade, assembly)



    ### Prepare GXBeam inputs  
    Omega0 = SVector(0.0, 0.0, -env.RS(t0))
    gravity0 = SVector(g*sin(azimuth0), -g*cos(azimuth0), 0.0)
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed



    ### TODO: I'm still debating whether or not I should include precone, tilt, and yaw in the assembly... what benefits and drawbacks would there be? What things would I have to change if I did? What things should be here if I didn't? => I would have to translate the velocity further. I think it would be easier to include the precone. 

    #TODO: Create a function to initialize the behavior. One for starting from no loading, one from steady state as below. Look at that, I wrote that down already. 

    ### GXBeam initial solution #TODO: I might make it an option to initialize from rest, or from steady state at the intial conditions (as current).


    # system, history0, converged = GXBeam.time_domain_analysis(assembly, tvec[1:1]; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, angular_velocity = Omega, gravity=gravity0, steady=false, initialize=true) 

    # gxhistory[1] = history0[end]



    ### Trying initializing from steady state spin
    system, _, _ = GXBeam.time_domain_analysis(assembly, tvec[1:1]; prescribed_conditions = prescribed_conditions, angular_velocity = Omega0, gravity=gravity0, steady=true, initialize=true, structural_damping) #Todo: Is there going to be some wonky behavior from passing in tvec[1:1]

    system, history0, converged = GXBeam.time_domain_analysis!(system, assembly, tvec[1:1]; prescribed_conditions=prescribed_conditions, distributed_loads, linear, angular_velocity = Omega0, reset_state=false, initialize=false, structural_damping, gravity=gravity0)
    gxhistory[1] = history0[end]



    ### Points to interpolate velocity, deflections, from structural to aerodynamic nodes. 
    interpolationpoints = create_interpolationpoints(assembly, blade) 

    delta = [interpolate_deflection(interpolationpoints[j], assembly, gxhistory[1]) for j = 1:na]
    # delta = [SVector(0.0, 0.0, 0.0) for j = 1:na]

    aeroV = [convert_velocities(blade, env, assembly, gxhistory[1], interpolationpoints, t0, j) for j = 1:na]  
    # aeroV = [SVector(0.0, 0.0, 0.0) for j = 1:na]
    
    def_theta = [interpolate_angle(interpolationpoints[j], assembly, gxhistory[1]) for j = 1:na] 
    # def_theta = [SVector(0.0, 0.0, 0.0) for j = 1:na]

    Vxmat = zeros(nt, na)
    Vymat = zeros(nt, na)

    Vxmat[1,:] .= Vxvec
    Vymat[1,:] .= Vyvec



    ### Iterate through time 
    for i = 2:nt
        t = tvec[i]
        dt = tvec[i] - tvec[i-1]

        #update azimuthal position
        azimuth[i] = env.RS(t)*dt + azimuth[i-1] #Euler step for azimuthal position. #TODO: Maybe do a better integration like a RK4 or something? I don't know if it matters much while I'm assuming the angular velocity is constant. 

        ### Update BEM inputs
        for j = 1:na
            ### Update sections with twist
            # if def_theta[j][1]>1
            #     @show i, j, def_theta[j]
            # end
            twistvec[j] = blade.twist[j]  
            # twistvec[j] = blade.twist[j] - def_theta[j][1] #Todo: Is this correct? I think this is what I had in the fixed point solution. -> A negative smooths the output. 
            # sections[j] = CCBlade.Section(rvec[j], chordvec[j], twist_displaced, blade.airfoils[j]) #Not displacing rvec because CCBlade uses that to calculate the tip and hub corrections, and the value should be based on relative to the distance along the blade to the hub or tip. (Although in dynamic movement, the tip losses would change dramatically.)

            ### Update base inflow velocities
            # Vxvec[j], Vyvec[j] = get_aerostructural_velocities(env, aeroV[j], t, rvec[j], azimuth[i], precone, tilt, yaw, hubht) #Todo. This should probably see the deflected radius. -> Or I could not subtract out the angular portion of the structural velocity.... then it would automatically see the deflected radius. 
            Vxvec[j], Vyvec[j] = Rotors.get_aerostructural_velocities(rotor, blade, env, t, j, azimuth[i], delta[j], def_theta[j], aeroV[j])

            # if j==300
            #     @show Vxvec[j], Vyvec[j]
            # end

            # if Vxvec[j]>330 || Vyvec[j]>330
            #     @show i, j, Vxvec[j]
            # end

            # if j==140
            #     @show Vxvec[j], Vyvec[j]
            # end

            # operatingpoints[j] = CCBlade.OperatingPoint(Vxvec[j], Vyvec[j], env.rho, pitch, env.mu, env.a)
            cchistory[i, j] = solve_BEM!(rotor, blade, env, j, Vxvec[j], Vyvec[j], pitch, xcc; twist = twistvec[j])

            # if i==2&&j==300
            #     @show i, j, cchistory[i,j].W, cchistory[i,j].phi
            #     @show Vxvec[j], Vyvec[j]
            #     @show cchistory[i,j].a
            #     @show cchistory[i,j].ap
            #     @show def_theta[j]
            #     println("")
            # end

            # if cchistory[i,j].W > 330&&j==300
            #     @show i, j, cchistory[i,j].W, cchistory[i,j].phi
            #     @show Vxvec[j], Vyvec[j]
            #     @show cchistory[i,j].a
            #     @show cchistory[i,j].ap
            #     @show def_theta[j]
            #     println("")

            #     #Todo: Problem 1) I don't think that the freestream velocity should ever be near 20, it should be at or below 10, right? Never above... 
            #     #Todo: Problem 2) why the heck is the inflow velocity jumping up to 900+. That violates the cauchy shwarz inequality. -> Crazy twist distribution. 
            # end
            #TODO: Write a solver that is initialized with the previous inflow angle.
        end

        Vxmat[i,:] .= Vxvec
        Vymat[i,:] .= Vyvec

        ### Solve BEM
        # cchistory[i, :] = CCBlade.solve.(Ref(rotor), sections, operatingpoints)  

        # for j = 1:na
        #     if cchistory[i,j].W < 0
        #         @show i, j, cchistory[i,j].W
        #     elseif cchistory[i,j].W > 330
        #         @show i, j, cchistory[i,j].W
        #     end
        # end


        ### Update Dynamic Stall model inputs 
        update_ds_inputs!(airfoils, p_ds, cchistory[i,:].W, cchistory[i,:].phi, twistvec, pitch, dt, turbine)



        
        ### Integrate Dynamic Stall model
        update_ds_states!(solver, airfoils, view(xds, i-1, :), view(xds, i, :), xds_idxs, p_ds, t, dt)
        # println("")
        # @show xds[i,end-31:end]



        ### Extract loads 
        
        extract_ds_loads!(airfoils, view(xds, i, :), xds_idxs, cchistory[i,:], p_ds, view(Cx, i, :), view(Cy, i, :), view(Cm, i, :))
    
        
        ### Dimensionalize
        for j = 1:na
            u_i = cchistory[i,j].W #TODO: Should I be normalizing by the actual nodal velocity, or the undistrubed nodal velocity. 
            qinf = 0.5*env.rho*u_i^2
            Fx[i,j] = Cx[i,j]*qinf*chordvec[j] 
            Fy[i,j] = Cy[i,j]*qinf*chordvec[j] 
            Mx[i,j] = Cm[i,j]*qinf*chordvec[j]^2
        end





        ### Update GXBeam loads  
        update_forces!(distributed_loads, view(Fx, i-1,:), Fy[i-1,:], Mx[i-1,:], blade, assembly) 

        Omega = SVector(0.0, 0.0, -env.RS(t))
        gravity = SVector(g*sin(azimuth[i-1]), -g*cos(azimuth[i-1]), 0.0)


        ### Solve GXBeam for time step #TODO: This function is taking a lot of time. 
        system, localhistory, converged = GXBeam.time_domain_analysis!(system, assembly, tvec[i-1:i]; prescribed_conditions=prescribed_conditions, distributed_loads, linear, angular_velocity = Omega, reset_state=false, initialize=false, structural_damping, gravity) #TODO: I feel like there is a faster way to accomplish this. Like, do I really need to reallocate Omega and gravity every time step? 


        ### Extract GXBeam outputs
        gxhistory[i] = localhistory[end] 




        ### Update aero inputs from structures.
        for j = 1:na
            delta[j] = interpolate_deflection(interpolationpoints[j], assembly, gxhistory[i])

            aeroV[j] = convert_velocities(blade, env, assembly, gxhistory[i], interpolationpoints, t, j)

            # if any(item->item<0, aeroV[j])
            # if aeroV[j][1]<0
            #     @show t, j, aeroV[j]
            # end

            def_theta[j] = interpolate_angle(interpolationpoints[j], assembly, gxhistory[i]) 
        end



        if verbose & (mod(i-1, speakiter)==0)
            println("")
            println("Simulation time: ", t)
        end

        if plotbool & (mod(i-1, plotiter)==0)
            # tipdef_x = [gxhistory[k].points[end].u[1] for k in eachindex(tvec[1:i])]
            # tipdef_y = [gxhistory[k].points[end].u[2] for k in eachindex(tvec[1:i])]
            # tipdef_z = [gxhistory[k].points[end].u[3] for k in eachindex(tvec[1:i])]
            # plt = plot(xaxis="Time (s)", yaxis="Tip Deflection", legend=:outerright)
            # plot!(tvec[1:i], tipdef_x, lab="X deflection")
            # plot!(tvec[1:i], tipdef_y, lab="Y deflection")
            # plot!(tvec[1:i], tipdef_z, lab="Z deflection")
            # display(plt)

            thetamat = zeros(i, 3)
            for k in 1:i
                thetamat[k,:] = gxhistory[k].points[end].theta
            end

            plt = plot(xaxis="Time (s)", yaxis="theta def", legend=:outerright)
            plot!(tvec[1:i], thetamat[:,1], lab="X deflection")
            plot!(tvec[1:i], thetamat[:,2], lab="Y deflection")
            plot!(tvec[1:i], thetamat[:,3], lab="Z deflection")
            display(plt)
        end
    end

    # return (Fx=Fx, Fy=Fy, Mx=Mx), cchistory, xds, gxhistory, def_theta
    return (Fx=Fx, Fy=Fy, Mx=Mx), cchistory, xds, gxhistory, Vxmat, Vymat
end
