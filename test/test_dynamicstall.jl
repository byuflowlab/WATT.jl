#=
Tests for the dynamic-stall wrapper layer in src/dynamicstallmodels.jl.

The underlying DS math lives in DynamicStallModels.jl; what WATT owns is the
translation between WATT's Blade and the per-airfoil DS state/input/parameter
buffers, plus the post-step extraction of (Cx, Cy) from (Cl, Cd, phi).

Scope:
  * initialize_ds_model          — buffer shapes, eltype, p_ds layout
  * dsmodel_initial_condition!   — non-NaN, in-bounds, turbine sign-flip
  * update_ds_inputs!            — inflow/aoa update + first-difference rates
  * update_ds_states!            — determinism across repeated calls
  * extract_ds_loads!            — Cx/Cy reproduce BEM rotation identity

Adam Cardoza
=#

using Test
using WATT, DynamicStallModels

const DS_ = DynamicStallModels

localpath = @__DIR__
cd(localpath)

include("fixtures/nrel5mw.jl")


@testset "Dynamic stall wrapper" begin

    nt = 5
    TF = Float64

    @testset "initialize_ds_model: shapes and layout" begin
        states, stateidx, y, p = WATT.initialize_ds_model(blade, nt, TF)

        # `states` is (nt × total-DS-states). `stateidx[i]` is the first state
        # index for airfoil i; consecutive entries differ by numberofstates(i).
        ns_total = DS_.numberofstates_total(blade.airfoils)
        @test size(states) == (nt, ns_total)
        @test eltype(states) == TF

        @test length(stateidx) == n
        @test stateidx[1] == 1
        for i in 1:(n-1)
            @test stateidx[i+1] - stateidx[i] == DS_.numberofstates(blade.airfoils[i].model)
        end

        # y_ds holds 4 inputs per airfoil (W, Udot, alpha, alphadot).
        @test length(y) == 4n
        @test eltype(y) == TF

        # p_ds holds [c[i], xcp[i]] interleaved.
        @test length(p) == 2n
        for i in 1:n
            @test p[2i - 1] == blade.c[i]
            @test p[2i]     == blade.xcp[i]
        end
    end


    @testset "dsmodel_initial_condition!: finite, turbine sign-flip" begin
        states, stateidx, y, p = WATT.initialize_ds_model(blade, nt, TF)

        mesh = (xds_idxs = stateidx, y_ds = y, p_ds = p)

        # Use a realistic operating point at every node so DS init has
        # something nontrivial to do.
        phi = fill(0.1, n)          # rad
        W   = fill(50.0, n)         # m/s
        pitch = 0.0
        t0 = 0.0

        # turbine = true path
        xds_t = view(states, 1, :)
        WATT.dsmodel_initial_condition!(xds_t, phi, W, mesh, blade, true, t0, pitch)
        @test all(isfinite, xds_t)

        # The function also writes into mesh.y_ds; confirm those entries are
        # finite and that the alpha entry (idx+3) carries the turbine sign.
        # alpha_turbine = -(twist + pitch - phi)  (turbine flips the sign)
        i_mid = max(2, n ÷ 2)
        alpha_t = y[4*(i_mid-1) + 3]
        @test isapprox(alpha_t, -((blade.twist[i_mid] + pitch) - phi[i_mid]))

        # turbine = false path (fresh buffers so we don't carry-over y_ds).
        states2, stateidx2, y2, p2 = WATT.initialize_ds_model(blade, nt, TF)
        mesh2 = (xds_idxs = stateidx2, y_ds = y2, p_ds = p2)
        xds_p = view(states2, 1, :)
        WATT.dsmodel_initial_condition!(xds_p, phi, W, mesh2, blade, false, t0, pitch)
        alpha_p = y2[4*(i_mid-1) + 3]
        @test isapprox(alpha_p, (blade.twist[i_mid] + pitch) - phi[i_mid])
        @test isapprox(alpha_t, -alpha_p)
    end


    @testset "update_ds_inputs!: rates are first differences" begin
        states, stateidx, y, p = WATT.initialize_ds_model(blade, nt, TF)

        # Prime y_ds with a known previous-step state so the first-difference
        # rates have something to subtract against.
        W_prev = fill(40.0, n)
        alpha_prev = fill(0.05, n)
        for j in 1:n
            idx = 4*(j-1)
            y[idx + 1] = W_prev[j]
            y[idx + 2] = 0.0
            y[idx + 3] = alpha_prev[j]
            y[idx + 4] = 0.0
        end

        W   = fill(50.0, n)
        phi = fill(0.1, n)
        pitch = 0.0
        dt = 0.01
        turbine = true

        WATT.update_ds_inputs!(blade.airfoils, view(y, :), W, phi, blade.twist, pitch, dt, turbine, blade)

        for j in 1:n
            idx = 4*(j-1)
            alpha_new = -((blade.twist[j] + pitch) - phi[j])   # turbine sign-flip

            @test y[idx + 1] == W[j]
            @test isapprox(y[idx + 2], (W[j] - W_prev[j]) / dt)
            @test isapprox(y[idx + 3], alpha_new)
            @test isapprox(y[idx + 4], (alpha_new - alpha_prev[j]) / dt)
        end
    end


    @testset "update_ds_states!: deterministic" begin
        states, stateidx, y, p = WATT.initialize_ds_model(blade, nt, TF)

        # Seed a plausible old state by running dsmodel_initial_condition! into row 1.
        mesh = (xds_idxs = stateidx, y_ds = y, p_ds = p)
        xds_old_view = view(states, 1, :)
        WATT.dsmodel_initial_condition!(xds_old_view, fill(0.1, n), fill(50.0, n), mesh, blade, true, 0.0, 0.0)

        # Snapshot inputs so two calls see identical (states_old, y_ds, p_ds).
        xds_old = copy(xds_old_view)
        y_snap  = copy(y)
        p_snap  = copy(p)

        states_a = similar(xds_old)
        states_b = similar(xds_old)

        WATT.update_ds_states!(WATT.RK4(), blade.airfoils, xds_old, states_a, stateidx,
                               copy(y_snap), copy(p_snap), 0.0, 0.01)
        WATT.update_ds_states!(WATT.RK4(), blade.airfoils, xds_old, states_b, stateidx,
                               copy(y_snap), copy(p_snap), 0.0, 0.01)

        @test states_a == states_b
        @test all(isfinite, states_a)
    end


    @testset "extract_ds_loads!: Cx/Cy match BEM rotation identity" begin
        states, stateidx, y, p = WATT.initialize_ds_model(blade, nt, TF)

        # Use the DS initial condition as a self-consistent state so get_loads
        # returns sensible (Cl, Cd, Cm) for each airfoil.
        mesh = (xds_idxs = stateidx, y_ds = y, p_ds = p)
        xds = view(states, 1, :)
        phi = fill(0.1, n)
        W   = fill(50.0, n)
        WATT.dsmodel_initial_condition!(xds, phi, W, mesh, blade, true, 0.0, 0.0)

        Cx = zeros(n); Cy = zeros(n); Cm = zeros(n)
        WATT.extract_ds_loads!(blade.airfoils, xds, stateidx, phi, y, p, Cx, Cy, Cm)

        # Independently compute (Cl, Cd) from DS.get_loads and apply the BEM
        # rotation. WATT.extract_ds_loads! must reproduce these to machine ε.
        for j in 1:n
            nsi1, nsi2 = DS_.state_indices(blade.airfoils[j].model, stateidx[j])
            xs = view(xds, nsi1:nsi2)
            ys = view(y, 4*(j-1)+1:4*j)
            ps = view(p, 2*(j-1)+1:2*j)
            Cl, Cd, _ = DS_.get_loads(blade.airfoils[j].model, blade.airfoils[j], xs, ys, ps)

            sphi, cphi = sincos(phi[j])
            Cx_ref =  Cl*cphi + Cd*sphi
            Cy_ref = -(Cl*sphi - Cd*cphi)

            @test isapprox(Cx[j], Cx_ref; rtol=1e-12, atol=1e-14)
            @test isapprox(Cy[j], Cy_ref; rtol=1e-12, atol=1e-14)
        end

        @test all(isfinite, Cx)
        @test all(isfinite, Cy)
        @test all(isfinite, Cm)
    end

end #End dynamic stall wrapper

nothing
