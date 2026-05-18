#=
Tests for the WATT-side functions that interface with GXBeam.jl as used by
WATT's aerostructural simulation (initialize_sim / run_sim! in
src/aerostructural.jl).

Scope: the wrapper layer WATT defines on top of GXBeam — rotation-parameter
conversions, beam arc-length computation, and the per-step force-update
function that maps aero loads onto GXBeam DistributedLoads. The actual
GXBeam math (system assembly, time stepping) is upstream and tested there.

The full end-to-end aerostructural regression lives in simple_NREL5MW.jl.

Adam Cardoza
=#

using Test
using WATT, GXBeam, OpenFASTTools, DynamicStallModels
using LinearAlgebra, StaticArrays, Statistics, StructArrays

of = OpenFASTTools
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)


# ---------- Shared setup: a real NREL 5MW blade + assembly ----------
# Matches the setup pattern in test_mesh.jl and test_environments.jl so
# the tests exercise the same data shapes run_sim! actually sees.

ofpath = "../data/openfast"
adblade = of.read_adblade("sn5_adblade.dat", ofpath)
edfile  = of.read_edfile("sn5_EDfile.dat", ofpath)
bdfile  = of.read_bdfile("sn5_BDfile.dat", ofpath)
bdblade = of.read_bdblade("sn5_BDblade.dat", ofpath)

aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder1.dat"))
aftypes[2] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder2.dat"))
aftypes[3] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU40_A17.dat"))
aftypes[4] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU35_A17.dat"))
aftypes[5] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU30_A17.dat"))
aftypes[6] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU25_A17.dat"))
aftypes[7] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU21_A17.dat"))
aftypes[8] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "NACA64_A17.dat"))

af_idx = Int.(adblade["BlAFID"])
afs    = aftypes[af_idx]

chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"]
rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
n    = length(rvec)

airfoils = StructArray{DS.Airfoil}(undef, n)
xcp = Vector{Float64}(undef, n)
for i = 1:n
    airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
end

blade    = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; rhub, rtip)
assembly = of.make_assembly(edfile, bdfile, bdblade)
ne       = length(assembly.elements)


@testset "GXBeam interface" begin

    @testset "retrieve_eulerangles" begin
        # WATT.retrieve_eulerangles(R) extracts (theta_x, theta_y, theta_z)
        # from a 3x3 rotation matrix. Used by WMPtoangle and (via
        # interpolate_angle) by update_mesh! every time step.

        # Identity → zero angles
        @test WATT.retrieve_eulerangles(Matrix{Float64}(I, 3, 3)) == SVector(0.0, 0.0, 0.0)

        # Pure x-axis rotation
        α = 0.3
        Rx = [1.0  0.0      0.0;
              0.0  cos(α)  -sin(α);
              0.0  sin(α)   cos(α)]
        θ = WATT.retrieve_eulerangles(Rx)
        @test isapprox(θ[1], α;   atol=1e-12)
        @test isapprox(θ[2], 0.0; atol=1e-12)
        @test isapprox(θ[3], 0.0; atol=1e-12)

        # Pure y-axis rotation
        β = 0.25
        Ry = [ cos(β)  0.0  sin(β);
                0.0    1.0   0.0;
              -sin(β)  0.0  cos(β)]
        θ = WATT.retrieve_eulerangles(Ry)
        @test isapprox(θ[1], 0.0; atol=1e-12)
        @test isapprox(θ[2], β;   atol=1e-12)
        @test isapprox(θ[3], 0.0; atol=1e-12)

        # Pure z-axis rotation
        γ = 0.4
        Rz = [cos(γ)  -sin(γ)  0.0;
              sin(γ)   cos(γ)  0.0;
              0.0      0.0     1.0]
        θ = WATT.retrieve_eulerangles(Rz)
        @test isapprox(θ[1], 0.0; atol=1e-12)
        @test isapprox(θ[2], 0.0; atol=1e-12)
        @test isapprox(θ[3], γ;   atol=1e-12)

        # Check a generic rotation - Todo: 
    end


    @testset "WMPtoangle" begin
        # WMPtoangle converts Wiener-Milenkovic rotation parameters (as
        # returned by GXBeam's PointState.theta) into Euler angles. Used by
        # interpolate_angle in mesh.jl, which feeds delta_theta to
        # get_aerostructural_velocities every time step.

        # Zero WMP → zero angles
        @test WATT.WMPtoangle(SVector(0.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)

        # Single-axis WMP must produce single-axis Euler angles. The exact
        # magnitude isn't checked here (the WMP↔angle relationship has a
        # nonlinear scaling), only that off-axis components stay zero.
        θx = WATT.WMPtoangle(SVector(0.5, 0.0, 0.0))
        @test isapprox(θx[2], 0.0; atol=1e-12)
        @test isapprox(θx[3], 0.0; atol=1e-12)
        @test θx[1] != 0.0 #Todo: We need to check that the angle is actually correct. 
        
        θy = WATT.WMPtoangle(SVector(0.0, 0.5, 0.0))
        @test isapprox(θy[1], 0.0; atol=1e-12)
        @test isapprox(θy[3], 0.0; atol=1e-12)
        @test θy[2] != 0.0

        θz = WATT.WMPtoangle(SVector(0.0, 0.0, 0.5))
        @test isapprox(θz[1], 0.0; atol=1e-12)
        @test isapprox(θz[2], 0.0; atol=1e-12)
        @test θz[3] != 0.0

        # No NaN / no Inf on a generic small rotation
        @test all(isfinite, WATT.WMPtoangle(SVector(0.1, -0.2, 0.05)))
    end


    @testset "get_bladelength_vector" begin
        # Returns the cumulative arc length from origin through every
        # structural point. Used by create_interpolationpoints to build
        # the aero↔structural coupling map.

        rgx = WATT.get_bladelength_vector(assembly)

        # One entry per structural point
        @test length(rgx) == length(assembly.points)

        # Monotonically increasing (root → tip)
        @test issorted(rgx)

        # First entry is the norm of the first point (distance from origin)
        @test isapprox(rgx[1], norm(assembly.points[1]))

        # Last entry equals the total arc length, computed independently
        # as the sum of segment norms (plus the initial offset from origin).
        np = length(assembly.points)
        L_independent = norm(assembly.points[1]) +
            sum(norm(assembly.points[i] - assembly.points[i-1]) for i in 2:np)
        @test isapprox(rgx[end], L_independent)

        # No NaN / no Inf
        @test all(isfinite, rgx)
    end


    @testset "update_forces!" begin
        # update_forces! maps per-aero-node loads (Fx axial, Fy tangential,
        # Mx torsion) onto a GXBeam DistributedLoads dict keyed by element
        # index. Called every time step inside run_sim! after take_aero_step!
        # before GXBeam.step_system!.
        #
        # Sign convention (from gxbeam.jl:406–407): GXBeam fz = -Fx (axial
        # thrust pushes the blade in -z of the beam frame), GXBeam fy = Fy.

        distributed_loads = Dict{Int, GXBeam.DistributedLoads{Float64}}()


        ### Zero loads: dict is fully populated; every element's
        ### distributed-load object is finite.
        Fx = zeros(n); Fy = zeros(n); Mx = zeros(n)
        WATT.update_forces!(distributed_loads, Fx, Fy, Mx, blade, assembly)

        @test length(distributed_loads) == ne
        for ielem in 1:ne
            @test haskey(distributed_loads, ielem)
            @test isa(distributed_loads[ielem], GXBeam.DistributedLoads{Float64})
            @test all(isfinite, distributed_loads[ielem].f1)
            @test all(isfinite, distributed_loads[ielem].f2)
        end


        ### Non-zero loads: distributed loads change and remain finite.
        # Capture the zero-load f1/f2 for comparison.
        f1_zero = [copy(distributed_loads[i].f1) for i in 1:ne]

        empty!(distributed_loads)
        Fx_const = 100.0 .* ones(n)   # 100 N/m axial thrust
        Fy_const = 50.0  .* ones(n)
        WATT.update_forces!(distributed_loads, Fx_const, Fy_const, Mx, blade, assembly)

        @test length(distributed_loads) == ne
        # At least one element's distributed-load integral must differ from
        # the zero-load case (otherwise the function isn't doing anything).
        nonzero_seen = false
        for ielem in 1:ne
            if !isapprox(distributed_loads[ielem].f1, f1_zero[ielem]; atol=1e-12)
                nonzero_seen = true
                break
            end
        end
        @test nonzero_seen

        # Everything still finite
        for ielem in 1:ne
            @test all(isfinite, distributed_loads[ielem].f1)
            @test all(isfinite, distributed_loads[ielem].f2)
        end


        ### Linearity: doubling the input loads doubles the output
        ### distributed-load integrals. The fit (DS.linear) is a linear
        ### operator, so f1/f2 of DistributedLoads must scale exactly.
        empty!(distributed_loads)
        WATT.update_forces!(distributed_loads, 2 .* Fx_const, 2 .* Fy_const, Mx, blade, assembly)
        f1_doubled = [copy(distributed_loads[i].f1) for i in 1:ne]

        empty!(distributed_loads)
        WATT.update_forces!(distributed_loads, Fx_const, Fy_const, Mx, blade, assembly)
        f1_single = [copy(distributed_loads[i].f1) for i in 1:ne]

        for ielem in 1:ne
            @test isapprox(f1_doubled[ielem], 2 .* f1_single[ielem]; rtol=1e-10)
        end
    end

end #End testing GXBeam interface
