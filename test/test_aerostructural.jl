#=
Fast-feedback tests for the coupled aerostructural path
(src/aerostructural.jl). Short tvec (≤10 steps) so this fails quickly when
the contract drifts; the long ~100 s regression lives in
simple_NREL5MW.jl and gates physics drift.

Scope:
  * initialize_sim — aerostates / mesh shape & eltype contract
  * run_sim!       — runs without NaN/Inf, structural deflection bounded
  * idempotency   — two fresh runs produce bitwise-identical outputs

Adam Cardoza
=#

using Test
using WATT, GXBeam, LinearAlgebra, StaticArrays

localpath = @__DIR__
cd(localpath)

include("fixtures/nrel5mw.jl")


@testset "Aerostructural simulation" begin

    nt = 8
    tvec = collect(range(0.0, 0.05, length=nt))   # 50 ms total, ~3° azimuth
    pitch = 0.0


    @testset "initialize_sim: shape & eltype contract" begin
        aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=false)

        for f in (:azimuth,)
            @test size(getfield(aerostates, f)) == (nt,)
        end
        for f in (:phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx)
            @test size(getfield(aerostates, f)) == (nt, n)
            @test eltype(getfield(aerostates, f)) == Float64
        end
        @test size(aerostates.xds, 1) == nt

        @test length(gxhistory) == nt

        # Mesh must carry every field run_sim! / take_aero_step! references.
        @test mesh isa WATT.SimMesh
        @test aerostates isa WATT.AeroStates{Float64}
        for k in (:interpolationpoints, :delta, :def_theta, :aerov, :xcc,
                  :xds_idxs, :y_ds, :p_ds, :assembly, :system,
                  :prescribed_conditions, :distributed_loads)
            @test k in propertynames(mesh)
        end
    end


    @testset "run_sim!: no NaN/Inf, bounded deflection" begin
        aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=false)
        WATT.run_sim!(rotor, blade, mesh, env_rated, tvec, aerostates, gxhistory; pitch, verbose=false)

        for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
            @test all(isfinite, getfield(aerostates, f))
        end

        # Tip deflection should be much smaller than blade length — sanity
        # bound, not a precision check.
        tip_u = gxhistory[end].points[end].u
        @test all(isfinite, tip_u)
        @test norm(tip_u) < rtip

        # Azimuth is monotonically non-decreasing
        for i in 2:nt
            @test aerostates.azimuth[i] >= aerostates.azimuth[i-1]
        end
    end


    @testset "idempotency: two fresh runs identical" begin
        a1, g1, m1 = WATT.initialize_sim(blade, assembly, tvec; verbose=false)
        WATT.run_sim!(rotor, blade, m1, env_rated, tvec, a1, g1; pitch, verbose=false)

        a2, g2, m2 = WATT.initialize_sim(blade, assembly, tvec; verbose=false)
        WATT.run_sim!(rotor, blade, m2, env_rated, tvec, a2, g2; pitch, verbose=false)

        for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
            @test getfield(a1, f) == getfield(a2, f)
        end
    end

    @testset "run_sim (non-mutating) matches initialize_sim + run_sim!" begin
        a1, g1, m1 = WATT.initialize_sim(blade, assembly, tvec; verbose=false)
        WATT.run_sim!(rotor, blade, m1, env_rated, tvec, a1, g1; pitch, verbose=false)

        a2, g2, _ = WATT.run_sim(rotor, blade, assembly, env_rated, tvec; pitch, verbose=false)

        for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
            @test getfield(a1, f) == getfield(a2, f)
        end
        @test length(g1) == length(g2)
    end

end #End aerostructural simulation

nothing
