

function initialize_ds_model(blade::Blade, nt, TF)

    airfoils = blade.airfoils
    n = length(airfoils) #Number of airfoils
    ns = DS.numberofstates_total(airfoils) #Total number of states

    states = Array{TF, }(undef, nt, ns) 
    y = Array{TF, 1}(undef, 4n)
    p = Vector{TF}(undef, 2n)
    stateidx = Vector{Int}(undef, n)
    tempx = 1
    
    for i in eachindex(airfoils)
        stateidx[i] = tempx 
        tempx += DS.numberofstates(airfoils[i].model)

        paramidx = 2*(i-1)+1:2*i
        p[paramidx] = [blade.c[i], blade.xcp[i]] 
    end
    
    return states, stateidx, y, p
end

function dsmodel_initial_condition!(xds, phi, W, mesh, blade::Blade, turbine::Bool, t0, pitch)

    for i in eachindex(blade.airfoils)
        nsi1, nsi2 = DS.state_indices(blade.airfoils[i].model, mesh.xds_idxs[i])
        
        alpha = (blade.twist[i] + pitch) - phi[i] 
        
        if turbine #Todo: This seems wrong... 
            alpha *= -1
        end

        env_idx = 4*(i-1)+1:4*i 
        ys = view(mesh.y_ds, env_idx)
        ys[1] = W[i] #Uvec[i]
        ys[2] = 0.0 #TODO: This will probably need to be updated later down the line. 
        ys[3] = alpha
        ys[4] = 0.0

        paramidx = 2*(i-1)+1:2*i
        ps = view(mesh.p_ds, paramidx)

        xds[nsi1:nsi2], _ = DS.initialize(blade.airfoils[i], [t0], ys, ps) 
    end
end

function update_ds_states!(solver::Solver, airfoils::AbstractVector{<:DS.Airfoil}, states_old, states_new, xds_idxs, y_ds, p_ds, t, dt)

    # @show y_ds[end-3:end]

    #=Note: The AeroDyn 15 theory book says that if alpha is greater than 90, or less than -90 then they shift it back

    if alpha>90
        alpha = 180-alpha
    elseif alpha<-90
        alpha = -180 - alpha
    end

    =#
    
    airfoils(states_old, states_new, xds_idxs, y_ds, p_ds, dt)

    # if isa(airfoil.model.detype, Indicial) #Indicial #todo: I need a way to either switch between, or enforce that ll of the dsmodels on a blade will be of a similar DEType. 
    #     #Pass in dt
    #     airfoils(states_old, states_new, xds_idxs, y_ds, dt)

    # else #Functional and Iteratives
    #     @error("Rotors isn't set up to handle functionals or iteratives yet. yet.")
    #     #Solve the State rate equations (Pass in t)
    # end

end

function update_ds_inputs!(airfoils::AbstractVector{<:Airfoil}, y_ds, W, phivec, twistvec, pitch, dt, turbine, blade::Blade)
    for j in eachindex(airfoils)
        idx = 4*(j-1)
        Udot = (W[j] - y_ds[idx+1])/dt #Calculate the inflow acceleration. 

        alpha = (twistvec[j] + pitch) - phivec[j]
        if turbine
            alpha *= -1
        end
        alphadot = (alpha - y_ds[idx+3])/dt

        # c = blade.c[j]
        # xcp = blade.xcp[j]
        # @show alpha
        # alpha = alpha + c*alphadot*(1 - 2*xcp) / 2 / W[j] #Theodorsen correction #Causes some hefty errors. 
        # @show alpha
        # println("") 

        y_ds[idx + 1] = W[j] #Update the inflow velocity
        y_ds[idx + 2] = Udot #Update the inflow acceleration
        y_ds[idx + 3] = alpha #Update the aoa
        y_ds[idx + 4] = alphadot #Update the angular velocity
    end
    # println("updated ds inputs")
end

function extract_ds_loads!(airfoils::AbstractVector{<:Airfoil}, states, state_idxs, phi, y_ds, p_ds, Cx, Cy, Cm)


    for j in eachindex(airfoils)
        nsi1, nsi2 = DS.state_indices(airfoils[j].model, state_idxs[j])
        xs = view(states, nsi1:nsi2) #section states
        ys = view(y_ds, 4*(j-1)+1:4*j) #section environment states
        ps = view(p_ds, 2*(j-1)+1:2*j) #section parameters

        
        Cl, Cd, Cm[j] = DS.get_loads(airfoils[j].model, airfoils[j], xs, ys, ps) 

        sphi, cphi = sincos(phi[j])
        Cx[j] = Cl*cphi + Cd*sphi
        Cy[j] = -(Cl*sphi - Cd*cphi)
    end
end

#Todo. I need a function to rotate the loads -> I don't need to because I can just load them into GXBeam as a follower load, which is in the local reference frame (the deformed reference frame). 
# function rotateloads()
# end