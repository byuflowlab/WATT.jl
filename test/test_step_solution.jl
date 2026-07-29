#=
Equivalence gate for the single-step coupled primitive (src/aerostructural.jl):
`step_solution!` / `initialize_from_state` / `run_from_state!`.

The make-or-break contract: restarting from a snapshot of a `run_sim!` march at step
`s` and marching the window `s → s+k` with `run_from_state!` must reproduce
`run_sim!`'s Newmark solution step-for-step. Because both paths use the identical
Newmark-β update (the `i<=2` seeding of `step_system!` reads the rates the `i>2`
recursion would reconstruct), agreement is to floating point (~1e-13 on the states).

AD-through-the-step tests (frozen-start ∂u/∂p) live in test_ad.jl, gated behind
WATT_AD_TESTS, reusing that file's compliance/mass parameter machinery.

Adam Cardoza
=#

using Test
using WATT, GXBeam, LinearAlgebra, StaticArrays

localpath = @__DIR__
cd(localpath)

include("fixtures/nrel5mw.jl")


@testset "step_solution! / run_from_state! equivalence with run_sim!" begin

    nt    = 10
    tvec  = collect(range(0.0, 0.05, length=nt))    # 50 ms, ~3° azimuth
    pitch = 0.0
    s     = 4                                        # snapshot step (well into the i>2 regime)

    # --- Reference: full march from rest. ---
    aref, gref, mref = WATT.initialize_sim(blade, assembly, tvec; verbose=false)
    WATT.run_sim!(rotor, blade, mref, env_rated, tvec, aref, gref; pitch, verbose=false)

    # --- Snapshot the native state at s (the objects a downstream caller keeps). ---
    state_s   = gref[s]
    xds_s     = aref.xds[s, :]
    azimuth_s = aref.azimuth[s]

    # --- Restart from s on a fresh mesh over the window tvec[s:end]. ---
    tvec_win = tvec[s:end]
    nwin     = length(tvec_win)
    _, _, mwin = WATT.initialize_sim(blade, assembly, tvec_win; verbose=false)   # allocation-only
    init = WATT.initialize_from_state(state_s, xds_s, azimuth_s, mwin, blade,
                                      env_rated, tvec_win, tvec[s])

    # Capture the windowed structural states (local index j ⇔ global step s+j-1).
    wstates    = Vector{Any}(undef, nwin)
    wstates[1] = state_s
    WATT.run_from_state!(init..., mwin, rotor, blade, env_rated, tvec_win;
                         pitch, out=(st, j) -> (wstates[j] = st))

    @testset "point states match run_sim! (u/θ/V/Ω)" begin
        for j in 2:nwin
            g   = s + j - 1
            wst = wstates[j]
            ref = gref[g]
            @test length(wst.points) == length(ref.points)
            for ip in eachindex(ref.points)
                @test isapprox(wst.points[ip].u,     ref.points[ip].u;     atol=1e-10, rtol=1e-8)
                @test isapprox(wst.points[ip].theta, ref.points[ip].theta; atol=1e-10, rtol=1e-8)
                @test isapprox(wst.points[ip].V,     ref.points[ip].V;     atol=1e-10, rtol=1e-8)
                @test isapprox(wst.points[ip].Omega, ref.points[ip].Omega; atol=1e-10, rtol=1e-8)
            end
        end
    end

    @testset "element internal loads match run_sim! (Fi/Mi)" begin
        for j in 2:nwin
            g = s + j - 1
            for ie in eachindex(gref[g].elements)
                @test isapprox(wstates[j].elements[ie].Fi, gref[g].elements[ie].Fi; atol=1e-4, rtol=1e-8)
                @test isapprox(wstates[j].elements[ie].Mi, gref[g].elements[ie].Mi; atol=1e-4, rtol=1e-8)
            end
        end
    end

    @testset "run_from_state! reaches the same final state as run_sim!" begin
        final = wstates[nwin]
        @test isapprox(final.points[end].u, gref[nt].points[end].u; atol=1e-10, rtol=1e-8)
    end

end #End equivalence gate

nothing
