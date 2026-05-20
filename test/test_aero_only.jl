#=
Tests for the aero-only simulation path in src/aero_only.jl.

Phase 3 fixed the `initialize` → `initialize_aero` rename together with the
two carry-over bugs (DS-init signature, missing `p_ds` in the mesh).
The contract assertions below are now active `@test` rather than `@test_broken`.

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

    @testset "initialize_aero: shape & eltype contract" begin
        aerostates, mesh = WATT.initialize_aero(blade, tvec)
        @test size(aerostates.azimuth) == (length(tvec),)
        @test size(aerostates.phi)     == (length(tvec), n)
        @test size(aerostates.Fx)      == (length(tvec), n)
        @test eltype(aerostates.phi) == Float64
        for k in (:delta, :def_theta, :aerov, :xcc, :xds_idxs, :y_ds, :p_ds)
            @test k in propertynames(mesh)
        end
    end


    @testset "simulate!: end-to-end finite, no NaN/Inf" begin
        aerostates, mesh = WATT.initialize_aero(blade, tvec)
        WATT.simulate!(aerostates, mesh, rotor, blade, env_rated, tvec; pitch)
        @test all(isfinite, aerostates.Fx)
        @test all(isfinite, aerostates.Fy)
    end

end #End aero-only simulation

nothing
