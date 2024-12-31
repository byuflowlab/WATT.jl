
#Todo: Move the DSM initialization type (the current ModelInit struct) to DSM. -> I'm not sure what I meant here. I think I had an initialization type previously, probably for how I was initializing the ds model. 
function initialize_ds_model(airfoils::AbstractVector{<:Airfoil}, nt; inittype=typeof(airfoils[1].c))
    n = length(airfoils) #Number of airfoils
    ns = DS.numberofstates_total(airfoils) #Total number of states

    states = Array{inittype, }(undef, nt, ns) 
    y = Array{inittype, 1}(undef, 4n)
    stateidx = Vector{Int}(undef, n)
    tempx = 1
    
    for i in eachindex(airfoils)
        stateidx[i] = tempx 
        tempx += DS.numberofstates(airfoils[i].model)
    end
    
    return states, stateidx, y
end

function dsmodel_initial_condition!(xds, phi, W, mesh, blade::Blade, turbine::Bool, t0, pitch)

    for i in eachindex(blade.airfoils)
        nsi1, nsi2 = DS.state_indices(blade.airfoils[i].model, mesh.xds_idxs[i])
        
        paramidx = 4*(i-1)+1:4*i 

        theta = (blade.twist[i] + pitch) - phi[i] #Todo: Isn't this alpha???
        
        if turbine
            theta *= -1
        end

        ys = view(mesh.p_ds, paramidx)
        ys[1] = W[i] #Uvec[i]
        ys[2] = 0.0 #TODO: This will probably need to be updated later down the line. 
        ys[3] = theta
        ys[4] = 0.0
        
        xds[nsi1:nsi2], _ = DS.initialize(blade.airfoils[i], [t0], ys) 
    end
end

function update_ds_states!(solver::Solver, airfoils::AbstractVector{<:DS.Airfoil}, states_old, states_new, xds_idxs, p_ds, t, dt)

    # @show p_ds[end-3:end]

    #=Note: The AeroDyn 15 theory book says that if alpha is greater than 90, or less than -90 then they shift it back

    if alpha>90
        alpha = 180-alpha
    elseif alpha<-90
        alpha = -180 - alpha
    end

    =#
    
    airfoils(states_old, states_new, xds_idxs, p_ds, dt)

    # if isa(airfoil.model.detype, Indicial) #Indicial #Todo: I need a way to either switch between, or enforce that ll of the dsmodels on a blade will be of a similar DEType. 
    #     #Pass in dt
    #     airfoils(states_old, states_new, xds_idxs, p_ds, dt)

    # else #Functional and Iteratives
    #     @error("Rotors isn't set up to handle functionals or iteratives yet. yet.")
    #     #Solve the State rate equations (Pass in t)
    # end

end

function update_ds_inputs!(airfoils::AbstractVector{<:Airfoil}, p_ds, W, phivec, twistvec, pitch, dt, turbine)
    for j in eachindex(airfoils)
        idx = 4*(j-1)
        Udot = (W[j] - p_ds[idx+1])/dt #Calculate the inflow acceleration. 

        alpha = (twistvec[j] + pitch) - phivec[j]
        if turbine
            alpha *= -1
        end
        alphadot = (alpha - p_ds[idx+3])/dt

        # if W[j]<0
        #     @show j, W[j]
        # end

        # if W[j] < 0
        #     @show j, W[j]
        # end

        p_ds[idx + 1] = W[j] #Update the inflow velocity
        p_ds[idx + 2] = Udot #Update the inflow acceleration
        p_ds[idx + 3] = alpha #Update the aoa
        p_ds[idx + 4] = alphadot #Update the angular velocity
    end
end

function extract_ds_loads!(airfoils::AbstractVector{<:Airfoil}, states, state_idxs, phi, p_ds, Cx, Cy, Cm)

    # loads = [0.0, 0.0, 0.0] #Todo: I need to figure out a way to return just the values. (Inside DSM)
    # @show typeof(ccstate)

    for j in eachindex(airfoils)
        nsi1, nsi2 = DS.state_indices(airfoils[j].model, state_idxs[j])
        xs = view(states, nsi1:nsi2)
        ys = view(p_ds, 4*(j-1)+1:4*j)

        # if ys[1]<1
        #     @show j, ys[1]
        # end

        
        Cl, Cd, Cm[j] = DS.get_loads(airfoils[j].model, airfoils[j], xs, ys) 

        # Cl = loads[1]
        # Cd = loads[2]
        # Cm[j] = loads[3]

        # cphi = cos(phi[j])
        # sphi = sin(phi[j])
        sphi, cphi = sincos(phi[j])
        Cx[j] = Cl*cphi + Cd*sphi
        Cy[j] = -(Cl*sphi - Cd*cphi)
    end
end

#Todo. I need a function to rotate the loads -> I don't need to because I can just load them into GXBeam as a follower load, which is in the local reference frame (the deformed reference frame). 
function rotateloads()
end