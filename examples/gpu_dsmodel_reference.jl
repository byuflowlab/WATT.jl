#=
Golden CPU trace generator for the GPU dynamic-stall prototype.

Builds the NREL 5MW blade and the prescribed (U, aoa) batch (see
`gpu_dsmodel_common.jl`), marches the CPU Beddoes-Leishman v3 (AeroDyn +
Gonzalez) model over every (section, sim) pair, and serializes the full state
history + loads to `examples/data/gpu_dsmodel_reference.jld2`. That file is the
versioned reference the GPU port validates against in `gpu_dsmodel_validate.jl`.

The march uses the same public DS routines WATT's `take_aero_step!` uses:
`DS.initialize` for the quasi-steady row-1 state, `DS.update_states_ADG!` for
each subsequent step (reads old states, writes new), and `DS.get_loads` for
(Cl, Cd, Cm). Only y[1]=U and y[3]=aoa are read by the ADG model, so Udot and
alphadot are held at zero — matching what the device kernel will consume.

Run:
    julia --project=examples examples/gpu_dsmodel_reference.jl

Adam Cardoza
=#

include(joinpath(@__DIR__, "gpu_dsmodel_common.jl"))

using JLD2
using Printf

const NSTATES = 32

# March a single 2D airfoil through its prescribed (U, aoa) history.
# Writes the (NSTATES × nt) state history and (3 × nt) load history into the
# supplied views. Returns nothing.
function march_section!(xds_hist_col, loads_col, airfoil, c, xcp, Uts, aoats, dt)
    nt = length(Uts)
    p  = [c, xcp]
    y  = [Uts[1], 0.0, aoats[1], 0.0]

    x0, l0 = _DS.initialize(airfoil, [0.0], y, p)
    xds_hist_col[:, 1] .= x0
    loads_col[:, 1]    .= l0

    xold = collect(x0)
    xnew = zero(xold)
    for i in 2:nt
        y[1] = Uts[i]; y[2] = 0.0; y[3] = aoats[i]; y[4] = 0.0
        _DS.update_states_ADG!(airfoil, xold, xnew, y, p, dt)
        Cl, Cd, Cm = _DS.get_loads(airfoil.model, airfoil, xnew, y, p)
        xds_hist_col[:, i] .= xnew
        loads_col[1, i] = Cl
        loads_col[2, i] = Cd
        loads_col[3, i] = Cm
        xold, xnew = xnew, xold
    end
    return nothing
end

function generate_reference()
    blade, _ = build_nrel5mw_blade()
    U, aoa, tvec, dt = make_ds_batch(blade)

    n_sections, n_sims, nt = size(U)
    @printf("Generating CPU DS reference: %d sections × %d sims × %d steps (dt=%.4g)\n",
            n_sections, n_sims, nt, dt)

    xds_hist = Array{Float64}(undef, NSTATES, n_sections, n_sims, nt)
    Cl = Array{Float64}(undef, n_sections, n_sims, nt)
    Cd = Array{Float64}(undef, n_sections, n_sims, nt)
    Cm = Array{Float64}(undef, n_sections, n_sims, nt)

    loads_col = Array{Float64}(undef, 3, nt)
    for s in 1:n_sims, j in 1:n_sections
        march_section!(view(xds_hist, :, j, s, :), loads_col,
                       blade.airfoils[j], blade.c[j], blade.xcp[j],
                       view(U, j, s, :), view(aoa, j, s, :), dt)
        Cl[j, s, :] .= loads_col[1, :]
        Cd[j, s, :] .= loads_col[2, :]
        Cm[j, s, :] .= loads_col[3, :]
    end

    n_nan = count(isnan, xds_hist) + count(isnan, Cl) + count(isnan, Cd) + count(isnan, Cm)
    n_nan > 0 && @warn "Reference contains $n_nan NaN entries — inspect scenarios/airfoils."

    scenario_names = [String(sc.name) for sc in DS_SCENARIOS]
    mkpath(dirname(REF_FILE))
    @save REF_FILE U aoa tvec dt xds_hist Cl Cd Cm scenario_names n_sections n_sims nt
    @printf("Saved reference to %s (%.1f MB)\n", REF_FILE, filesize(REF_FILE) / 1e6)
    return nothing
end

generate_reference()
