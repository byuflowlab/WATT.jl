#=
Full end-to-end regression test for WATT.jl.

Loads a reference simulation (NREL 5MW, ~100 s) saved when WATT was known to be
working correctly, recomputes it with the current source, and asserts that the
outputs match the saved reference within tight tolerance. This is the safety
net for the Phase 1+ rework — any drift here means cleanup changed physics.

Reference data: simpleNREL5_100s.jld2

Adam Cardoza 2026-05-15
=#

using Test, WATT, JLD2, GXBeam, FLOWMath, DelimitedFiles, StaticArrays

localpath = @__DIR__
cd(localpath)

const RTOL = 1e-8
const ATOL = 1e-10

# Reference fields:
#   blade, assembly, tvec_r, rotor_r, aerostates, gxhistory, mesh
# Environment scalars (used to reconstruct `env` in-test — closure-based
# SimpleEnvironment cannot be JLD2-serialized directly):
#   rho, mu, a, U, Omega, shearexp
# @load "simpleNREL5_100s.jld2" blade assembly tvec_r aerostates gxhistory mesh rotor_r \
#                               rho mu a U Omega shearexp
@load "simpleNREL5_100s.jld2" blade assembly tvec_r aerostates gxhistory mesh rotor_r env

env_ = env # Because env is a closure-based SimpleEnvironment, we can't directly compare it to a freshly constructed one. So we load it as env_, then reconstruct env from the scalars and functions in the file. ->  #Todo: This whole thing is a closure... and I don't know if I like that. It's providing some issues. 

RPM = 11.44
omega = RPM * 2pi / 60

rho = env_.rho
mu = env_.mu
a = env_.a
shearexp = env_.shearexp



turbfile = "../data/openfast/TurbSim.dat"

turb = readdlm(turbfile, skipstart=4)
n, m = size(turb)
tvec_turb = range(turb[1, 1], turb[n, 1], length=n) #Because the file doesn't get the time vector correctly. 
Ufit = Akima(tvec_turb, turb[:, 2])
Vfit = Akima(tvec_turb, turb[:, 5])
Wfit = Akima(tvec_turb, turb[:, 6])

ufun(t) = SVector(Ufit(t), Vfit(t), Wfit(t))
omegafun(t) = SVector(0.0, 0.0, 0.0)
udotfun(t) = SVector(0.0, 0.0, 0.0) #Note: Could probably use the gradient function. 
omegadotfun(t) = SVector(0.0, 0.0, 0.0)
# Vinf(t) = Ufit(t) 
Vinf(t) = sqrt(Ufit(t)^2 + Vfit(t)^2)
RS(t) = omega
Vinfdot(t) = 0.0
RSdot(t) = 0.0
env = WATT.SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)

aerostates_ref = aerostates
gxhistory_ref = gxhistory


"""
Compare two GXBeam AssemblyStates field-by-field at a single time index.
Returns true if all per-point and per-element kinematic + load fields agree.
"""
function assembly_states_match(a, b; rtol=RTOL, atol=ATOL)
    length(a.points) == length(b.points) || return false
    length(a.elements) == length(b.elements) || return false
    for (pa, pb) in zip(a.points, b.points)
        isapprox(pa.u,     pb.u;     rtol=rtol, atol=atol) || return false
        isapprox(pa.theta, pb.theta; rtol=rtol, atol=atol) || return false
    end
    for (ea, eb) in zip(a.elements, b.elements)
        isapprox(ea.u,     eb.u;     rtol=rtol, atol=atol) || return false
        isapprox(ea.theta, eb.theta; rtol=rtol, atol=atol) || return false
        isapprox(ea.Fi,    eb.Fi;    rtol=rtol, atol=atol) || return false
        isapprox(ea.Mi,    eb.Mi;    rtol=rtol, atol=atol) || return false
    end
    return true
end


@testset "Full WATT simulation regression" begin

    @testset "initialize_sim: shape & eltype contract" begin
        aerostates_new, gxhistory_new, mesh_new = WATT.initialize_sim(blade, assembly, tvec_r; verbose=false)

        for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
            @test size(getfield(aerostates_new, f)) == size(getfield(aerostates_ref, f))
            @test eltype(getfield(aerostates_new, f)) == eltype(getfield(aerostates_ref, f))
        end
        @test length(gxhistory_new) == length(gxhistory_ref)
        @test eltype(gxhistory_new) == eltype(gxhistory_ref)
    end


    @testset "run_sim!: golden-value regression" begin
        aerostates_new, gxhistory_new, mesh_new = WATT.initialize_sim(blade, assembly, tvec_r; verbose=false)
        WATT.run_sim!(rotor_r, blade, mesh_new, env, tvec_r, aerostates_new, gxhistory_new; verbose=false)

        @testset "aerostates fields" begin
            for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
                @test isapprox(getfield(aerostates_new, f), getfield(aerostates_ref, f); rtol=RTOL, atol=ATOL)
            end
        end

        @testset "gxhistory spot-checks" begin
            nt = length(gxhistory_new)
            for i in (1, max(1, div(nt, 2)), nt)
                @test assembly_states_match(gxhistory_new[i], gxhistory_ref[i])
            end
        end

        @testset "no NaN / no Inf" begin
            for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
                @test all(isfinite, getfield(aerostates_new, f))
            end
        end

        @testset "integrated thrust & torque (last step)" begin
            Fx_new = aerostates_new.Fx[end, :]
            Fy_new = aerostates_new.Fy[end, :]
            Fx_old = aerostates_ref.Fx[end, :]
            Fy_old = aerostates_ref.Fy[end, :]

            T_new, Q_new = WATT.thrusttorque(blade.r, Fx_new, Fy_new, rotor_r.B)
            T_old, Q_old = WATT.thrusttorque(blade.r, Fx_old, Fy_old, rotor_r.B)

            @test isapprox(T_new, T_old; rtol=RTOL, atol=ATOL)
            @test isapprox(Q_new, Q_old; rtol=RTOL, atol=ATOL)
        end
    end


    @testset "idempotency: two fresh runs give identical results" begin
        a1, g1, m1 = WATT.initialize_sim(blade, assembly, tvec_r; verbose=false)
        WATT.run_sim!(rotor_r, blade, m1, env, tvec_r, a1, g1; verbose=false)

        a2, g2, m2 = WATT.initialize_sim(blade, assembly, tvec_r; verbose=false)
        WATT.run_sim!(rotor_r, blade, m2, env, tvec_r, a2, g2; verbose=false)

        for f in (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
            @test getfield(a1, f) == getfield(a2, f)
        end
    end

end
