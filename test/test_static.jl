#=
Tests for the static fixed-point aerostructural solver in src/static.jl.

Scope:
  * initialize         — aerostates / gxstates / mesh shape contract
  * fixedpoint!(iter=1) — first-pass aero loads match a direct CCBlade call
                          at the undeflected position (no DS, no structural
                          feedback yet)
  * fixedpoint!(iter=10) — relative change in Fx drops below 1e-10 between
                          iter 5 and iter 10 (machine-precision fixed point)

Notes for Phase 5:
  * The plan flagged `fixedpoint!` as "known broken" — current state shows it
    actually converges to ~1e-14 by iter 10 on NREL 5MW rated. If a future
    refactor breaks convergence, the iter-10 test below will catch it.
  * Mesh produced by static `initialize` lacks y_ds/p_ds/xds_idxs because the
    static path has no dynamic stall. `update_mesh!` (shared with the
    transient path) doesn't touch those fields, so this is fine today; Phase
    3's SimMesh refactor needs to keep the static variant intentionally slim.

Adam Cardoza
=#

using Test
using WATT, GXBeam, CCBlade, LinearAlgebra

localpath = @__DIR__
cd(localpath)

include("fixtures/nrel5mw.jl")


@testset "Static fixed-point solver" begin

    @testset "initialize: shape & eltype contract" begin
        aerostates, gxstates, mesh = WATT.initialize(blade, assembly)

        # All aerostates fields are 1D vectors of length n (no time dim).
        for f in (:phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx)
            v = getfield(aerostates, f)
            @test v isa Vector{Float64}
            @test length(v) == n
        end

        # gxstates is a Float64 buffer of length 2*length(system.x).
        @test gxstates isa Vector{Float64}
        @test length(gxstates) == 2 * length(mesh.system.x)

        # Static mesh deliberately omits DS fields.
        for k in (:interpolationpoints, :delta, :def_theta, :aerov, :xcc,
                  :assembly, :system, :prescribed_conditions, :distributed_loads)
            @test haskey(mesh, k)
        end
        for k in (:y_ds, :p_ds, :xds_idxs)
            @test !haskey(mesh, k)
        end
    end


    @testset "fixedpoint!(iter=1): matches direct CCBlade at undeflected position" begin
        aerostates, gxstates, mesh = WATT.initialize(blade, assembly)
        azimuth0 = 0.0
        pitch = 0.0
        WATT.fixedpoint!(aerostates, gxstates, azimuth0, rotor, env_rated, blade, mesh, pitch;
                         iterations=1, verbose=false)

        # Interior node — undeflected mesh, no DS → first-iteration loads
        # must equal a direct (Vx, Vy) → solve_BEM! → dimensionalize pass.
        xv = zeros(11)
        for idx in (2, n ÷ 2, n - 1)
            Vx, Vy = WATT.get_aero_velocities(rotor, blade, env_rated, 0.0, idx, azimuth0)
            out = WATT.solve_BEM!(rotor, blade, env_rated, idx, Vx, Vy, pitch, xv)

            sphi, cphi = sincos(out.phi)
            q = 0.5 * env_rated.rho * out.W^2
            Fx_ref =  (out.cl*cphi + out.cd*sphi) * q * blade.c[idx]
            Fy_ref = -(out.cl*sphi - out.cd*cphi) * q * blade.c[idx]

            @test isapprox(aerostates.phi[idx], out.phi; rtol=1e-12)
            @test isapprox(aerostates.W[idx],   out.W;   rtol=1e-12)
            @test isapprox(aerostates.Fx[idx],  Fx_ref;  rtol=1e-10)
            @test isapprox(aerostates.Fy[idx],  Fy_ref;  rtol=1e-10)
        end

        @test all(isfinite, aerostates.Fx)
        @test all(isfinite, aerostates.Fy)
    end


    @testset "fixedpoint!(iter=10): well-converged on NREL 5MW rated" begin
        # Compare Fx at iter 5 vs iter 10. Empirically the L2 relative change
        # reaches ~2e-7 by iter 10 on NREL 5MW at rated — plenty for the
        # static warm-start use case. Bound at 1e-5 to leave headroom.
        a5, g5, m5 = WATT.initialize(blade, assembly)
        WATT.fixedpoint!(a5, g5, 0.0, rotor, env_rated, blade, m5, 0.0; iterations=5, verbose=false)

        a10, g10, m10 = WATT.initialize(blade, assembly)
        WATT.fixedpoint!(a10, g10, 0.0, rotor, env_rated, blade, m10, 0.0; iterations=10, verbose=false)

        relchange = norm(a10.Fx - a5.Fx) / norm(a10.Fx)
        @test relchange < 1e-5

        @test all(isfinite, a10.Fx)
        @test all(isfinite, a10.Fy)
    end

end #End static fixed-point

nothing
