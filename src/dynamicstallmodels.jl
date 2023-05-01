

function initialize_DS_model(airfoils::AbstractVector{<:Airfoil}, turbine::Bool, nt, tvec, ccstate, Udotvec, alphadotvec, twistvec, pitch)
    n = length(airfoils) #Number of airfoils
    ns = DS.numberofstates_total(airfoils) #Total number of states

    states = Array{eltype(Udotvec), 2}(undef, nt, ns) #Todo: This is going to need to be changed when we start doing optimization. 
    # y = DS.initialize_environment(Uvec[1], Udotvec[1], alphavec[1], alphadotvec[1], n) 
    y = Array{eltype(Udotvec), 1}(undef, 4n)
    stateidx = Vector{Int}(undef, n)
    tempx = 1
    # loads = Array{eltype(Udotvec), 2}(undef, nt, 3n)
    # println("Got here. ")
    # @show typeof(ccstate)
    # @show fieldnames(typeof(ccstate))
    # @show ccstate.indices
    # @show ccstate.phi
    # @show typeof(ccstate[1])
    # @show size(ccstate)
    for i in eachindex(airfoils)
        stateidx[i] = tempx
        nsi1, nsi2 = DS.state_indices(airfoils[i].model, stateidx[i])
        # loadidx = 3*(i-1)+1:3i
        paramidx = 4*(i-1)+1:4*i 
        # println("Making it here. ")
        if turbine
            theta = -((twistvec[i] + pitch) - ccstate[i].phi) #TODO: I wonder if I should do this in the simulate section. 
        else
            theta = ((twistvec[i] + pitch) - ccstate[i].phi)
        end

        ys = view(y, paramidx)
        ys[1] = ccstate[i].W #Uvec[i]
        ys[2] = Udotvec[i]
        ys[3] = theta
        ys[4] = alphadotvec[i]
        # states[1,nsi1:nsi2], loads[1,loadidx] = DS.initialize(airfoils[i], tvec, ys) 
        states[1,nsi1:nsi2], _ = DS.initialize(airfoils[i], tvec, ys) 
        tempx += DS.numberofstates(airfoils[i].model)
    end
    # println(" Got here. 2.0")
    return states, stateidx, y
end

function update_ds_states!(solver::Solver, airfoils::AbstractVector{<:DS.Airfoil}, states_old, states_new, xds_idxs, p_ds, t, dt)
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

        alpha = ((twistvec[j] + pitch) - phivec[j])
        if turbine
            alpha *= -1
        end
        alphadot = (alpha - p_ds[idx+3])/dt

        # if W[j]<0
        #     @show j, W[j]
        # end

        p_ds[idx + 1] = W[j] #Update the inflow velocity
        p_ds[idx + 2] = Udot #Update the inflow acceleration
        p_ds[idx + 3] = alpha #Update the aoa
        p_ds[idx + 4] = alphadot #Update the angular velocity
    end
end

function extract_ds_loads!(airfoils::AbstractVector{<:Airfoil}, states, state_idxs, ccstate, p_ds, Cx, Cy, Cm)

    # loads = [0.0, 0.0, 0.0] #Todo: I need to figure out a way to return just the values. (Inside DSM)
    # @show typeof(ccstate)

    for j in eachindex(airfoils)
        nsi1, nsi2 = DS.state_indices(airfoils[j].model, state_idxs[j])
        xs = view(states, nsi1:nsi2)
        ys = view(p_ds, 4*(j-1)+1:4*j)

        
        Cl, Cd, Cm[j] = DS.get_loads(airfoils[j].model, airfoils[j], xs, ys) 

        # Cl = loads[1]
        # Cd = loads[2]
        # Cm[j] = loads[3]

        cphi = cos(ccstate[j].phi)
        sphi = sin(ccstate[j].phi)
        Cx[j] = Cl*cphi + Cd*sphi
        Cy[j] = -(Cl*sphi - Cd*cphi)
    end
end

#Todo: I need a function to rotate the loads
function rotateloads()
end