

"""
    initialize_ds_model(blade, nt, TF) -> (states, stateidx, y, p)

Allocate the per-section state, environment, parameter, and index buffers
needed to time-march the dynamic-stall models attached to `blade.airfoils`.

The state width per section depends on the DS model (e.g. Beddoes–Leishman
vs. NoModel), so the total state count is queried via
`DS.numberofstates_total(airfoils)` and the per-section offsets are recorded
in `stateidx`.

**Arguments**
- `blade::Blade`: blade carrying the per-section `DS.Airfoil` instances.
- `nt::Int`: number of time steps to allocate history for.
- `TF::Type`: element type for the allocated buffers (e.g. `Float64`, or a
  `ForwardDiff.Dual` type so the state history can carry duals).

**Returns**
- `states::Matrix{TF}`: `nt × ns` DS state history (`ns` is the total state count).
- `stateidx::Vector{Int}`: starting column index in `states` for each section.
- `y::Vector{TF}`: length-`4n` environment vector (one `[W, Ẇ, α, α̇]` block per section).
- `p::Vector{TF}`: length-`2n` parameter vector (`[chord, xcp]` per section), pre-filled.
"""
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

"""
    dsmodel_initial_condition!(xds, phi, W, mesh, blade, turbine, t0, pitch)

Populate the first row of the DS state history `xds` with each section's
quasi-steady initial condition. For every section, the environment block
in `mesh.y_ds` is set to `[W, 0, α, 0]` (no inflow acceleration, no
α̇), and `DS.initialize` is called on the section's airfoil to compute
the equilibrium DS states at `t0`.

Mutates `xds` and `mesh.y_ds` in place.

**Arguments**
- `xds`: DS state matrix; row 1 is written.
- `phi`: per-section inflow angles (radians).
- `W`: per-section relative-wind magnitudes.
- `mesh`: simulation mesh holding `xds_idxs`, `y_ds`, and `p_ds` buffers.
- `blade::Blade`: blade with airfoils and twist distribution.
- `turbine::Bool`: when `true`, α is negated (turbine sign convention).
- `t0`: initial time.
- `pitch`: blade pitch angle (radians).
"""
function dsmodel_initial_condition!(xds, phi, W, mesh, blade::Blade, turbine::Bool, t0, pitch)

    for i in eachindex(blade.airfoils)
        nsi1, nsi2 = DS.state_indices(blade.airfoils[i].model, mesh.xds_idxs[i])
        
        alpha = (blade.twist[i] + pitch) - phi[i] 
        
        if turbine 
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

"""
    update_ds_states!(solver, airfoils, states_old, states_new, xds_idxs, y_ds, p_ds, t, dt)

Advance the dynamic-stall states from `states_old` to `states_new` over a
single time step `dt`, applying each section's DS model in turn. The
underlying step is dispatched by calling the `airfoils` vector as a
function — DynamicStallModels.jl defines `(::AbstractVector{<:Airfoil})(...)`
to apply the appropriate per-model update on that call.

`solver` and `t` are accepted for signature uniformity with the rest of
WATT's time-stepping API; the current implementation hands off entirely to
the indicial DS update inside DynamicStallModels.

Mutates `states_new` in place.

**Arguments**
- `solver::Solver`: WATT time-integrator (currently informational).
- `airfoils`: per-section `DS.Airfoil` instances.
- `states_old`: previous DS state row.
- `states_new`: DS state row to be written.
- `xds_idxs`: per-section starting indices into the state vector.
- `y_ds`: environment vector (`[W, Ẇ, α, α̇]` per section).
- `p_ds`: parameter vector (`[chord, xcp]` per section).
- `t`: current time.
- `dt`: step size.
"""
function update_ds_states!(solver::Solver, airfoils::AbstractVector{<:DS.Airfoil}, states_old, states_new, xds_idxs, y_ds, p_ds, t, dt)


    #=Note: The AeroDyn 15 theory book says that if alpha is greater than 90, or less than -90 then they shift it back

    if alpha>90
        alpha = 180-alpha
    elseif alpha<-90
        alpha = -180 - alpha
    end

    =#
    
    airfoils(states_old, states_new, xds_idxs, y_ds, p_ds, dt)

    # if isa(airfoil.model.detype, Indicial) #Indicial #todo: I need a way to either switch between, or enforce that all of the dsmodels on a blade will be of a similar DEType. -> the way I've set up the airfoils methods, it might be able to automatically switch... 
    #     #Pass in dt
    #     airfoils(states_old, states_new, xds_idxs, y_ds, dt)

    # else #Functional and Iteratives
    #     @error("WATT isn't set up to handle functionals or iteratives yet. yet.")
    #     #Solve the State rate equations (Pass in t)
    # end

end

"""
    update_ds_inputs!(airfoils, y_ds, W, phivec, twistvec, pitch, dt, turbine, blade)

Refresh the environment vector `y_ds` ahead of a DS state update. For each
section, computes the current angle of attack and uses backward
differences against the previous-step values stored in `y_ds` to estimate
the inflow acceleration `Ẇ` and pitch rate `α̇`. The new
`[W, Ẇ, α, α̇]` block overwrites the previous one in place.

The Theodorsen-style midchord correction
`α + c·α̇·(1 - 2·xcp)/(2·W)` is intentionally disabled (it produced
large residuals in past testing).

**Arguments**
- `airfoils`: per-section `DS.Airfoil` instances (used for dispatch / length).
- `y_ds`: environment vector to mutate.
- `W`: per-section relative-wind magnitudes at the new time.
- `phivec`: per-section inflow angles at the new time.
- `twistvec`: per-section twist angles (radians).
- `pitch`: blade pitch angle (radians).
- `dt`: step size used for the backward-difference rates.
- `turbine`: when `true`, α is negated (turbine sign convention).
- `blade::Blade`: blade (currently unused, kept for signature symmetry).
"""
function update_ds_inputs!(airfoils::AbstractVector{<:Airfoil}, y_ds, W, phivec, twistvec, pitch, dt, turbine, blade::Blade)
    for j in eachindex(airfoils)
        idx = 4*(j-1)
        Udot = (W[j] - y_ds[idx+1])/dt #Calculate the inflow acceleration. 

        alpha = (twistvec[j] + pitch) - phivec[j]
        if turbine
            alpha *= -1
        end
        alphadot = (alpha - y_ds[idx+3])/dt

        # alpha = alpha + c*alphadot*(1 - 2*xcp) / 2 / W[j] #Theodorsen correction #Causes some hefty errors. 

        y_ds[idx + 1] = W[j] #Update the inflow velocity
        y_ds[idx + 2] = Udot #Update the inflow acceleration
        y_ds[idx + 3] = alpha #Update the aoa
        y_ds[idx + 4] = alphadot #Update the angular velocity
    end
end

"""
    extract_ds_loads!(airfoils, states, state_idxs, phi, y_ds, p_ds, Cx, Cy, Cm)

Pull `(Cl, Cd, Cm)` from each section's DS model given the current state
row, then rotate the in-plane force coefficients into the rotor frame
using the inflow angle `phi[j]`:

    Cx = Cl·cos(φ) + Cd·sin(φ)
    Cy = -(Cl·sin(φ) - Cd·cos(φ))

The negation on `Cy` puts tangential loading into the turbine torque sign
convention used elsewhere in WATT.

Mutates `Cx`, `Cy`, and `Cm` in place.

**Arguments**
- `airfoils`: per-section `DS.Airfoil` instances.
- `states`: current DS state row.
- `state_idxs`: per-section starting indices into the state vector.
- `phi`: per-section inflow angles (radians).
- `y_ds`: environment vector (`[W, Ẇ, α, α̇]` per section).
- `p_ds`: parameter vector (`[chord, xcp]` per section).
- `Cx`, `Cy`, `Cm`: output vectors written in place.
"""
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