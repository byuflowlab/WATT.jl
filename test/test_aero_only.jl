#=
Tests for the aero-only simulation path in src/aero_only.jl.

Current status (Phase 2 baseline): aero-only `initialize` is broken — it calls
an outdated `initialize_ds_model(blade.airfoils, nt; inittype=...)` signature
that no longer exists, and the mesh NamedTuple it builds is missing `p_ds`
which `take_aero_step!` → `update_ds_states!` requires. Phase 1 only fixed the
`.c` field access on line 71; the deeper signature/mesh mismatch is deferred
to Phase 3 (API consolidation) / Phase 5 (static solver integration).

This file is therefore mostly @test_broken: it records the *intended*
assertions so that when the path is repaired, the maintainer flips
@test_broken → @test and immediately knows the contract.

Adam Cardoza
=#

using Test
using WATT

localpath = @__DIR__
cd(localpath)

include("fixtures/nrel5mw.jl")


@testset "Aero-only simulation" begin

    tvec = collect(range(0.0, 0.05, length=5))   # 5 steps, 12.5 ms dt
    pitch = 0.0

    # ---------------------------------------------------------------------
    # Documented broken state: `initialize` crashes on a MethodError because
    # `initialize_ds_model` was changed to (blade, nt, TF) but aero_only.jl
    # still calls the old (blade.airfoils, nt; inittype=...) form. When this
    # is fixed in Phase 3/5, the @test_throws below will flip to a failure —
    # at that point swap it for the @test_skip block that follows.
    # ---------------------------------------------------------------------
    @testset "initialize: currently broken (Phase 3/5 target)" begin
        @test_throws MethodError WATT.initialize(blade, tvec)
    end


    # ---------------------------------------------------------------------
    # Intended contract — once initialize works, enable these.
    # ---------------------------------------------------------------------
    @testset "initialize: shape & eltype contract (skipped until fix)" begin
        @test_skip begin
            aerostates, mesh = WATT.initialize(blade, tvec)
            @test size(aerostates.azimuth) == (length(tvec),)
            @test size(aerostates.phi)     == (length(tvec), n)
            @test size(aerostates.Fx)      == (length(tvec), n)
            @test eltype(aerostates.phi) == Float64
            for k in (:delta, :def_theta, :aerov, :xcc, :xds_idxs, :y_ds, :p_ds)
                @test haskey(mesh, k)
            end
            true
        end
    end


    @testset "simulate!: end-to-end (deferred — needs p_ds in mesh)" begin
        @test_broken begin
            aerostates, mesh = WATT.initialize(blade, tvec)
            WATT.simulate!(aerostates, mesh, rotor, blade, env_rated, tvec; pitch)
            all(isfinite, aerostates.Fx) && all(isfinite, aerostates.Fy)
        end
    end

end #End aero-only simulation

nothing
