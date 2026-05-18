#=
Test the meshing and mesh transfer functions. 

=#

using Test
using WATT, GXBeam, OpenFASTTools, DynamicStallModels
using LinearAlgebra, StaticArrays, StructArrays

of = OpenFASTTools
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)

@testset "Mesh" begin


    @testset "Mesh generation" begin

        ### Test find_paired_indices()
        #inbounds tests
        rvec = collect(0:10)

        r = .5
        pair = WATT.find_point_indices(rvec, r)
        pair_gold = (1, 2)
        @test pair==pair_gold

        r = pi
        pair = WATT.find_point_indices(rvec, r)
        pair_gold = (4, 5)
        @test pair==pair_gold

        #Boundary tests
        r = 0
        pair = WATT.find_point_indices(rvec, r)
        pair_gold = (1, 2)
        @test pair==pair_gold

        r = 10
        pair = WATT.find_point_indices(rvec, r)
        pair_gold = (10, 11)
        @test pair==pair_gold

        #Out of bounds tests. Current behavior: warn + clamp the returned
        #pair to the first or last valid bracket. We codify both the
        #fallback pair AND the warning here so a future flip to "error on
        #out of bounds" surfaces as an obvious test failure rather than a
        #silent behavior change.
        #Todo: Decide whether to keep warn-and-clamp or switch to throwing
        #an error. find_point_indices is only called at simulation setup
        #(via create_interpolationpoints in initialize_sim), so a hard
        #error there is cheap and catches misconfigurations loudly.
        r = -1
        local pair
        @test_logs (:warn, r"not found") (pair = WATT.find_point_indices(rvec, r))
        @test pair == (1, 2)

        r = 4*pi
        @test_logs (:warn, r"not found") (pair = WATT.find_point_indices(rvec, r))
        @test pair == (10, 11)

        ### Test find_interpolation_percent
        #Test a random spot
        r = 4.75
        pair = WATT.find_point_indices(rvec, r)
        percent = WATT.find_interpolation_percent(rvec, pair, r)
        percent_gold = 0.75
        @test percent==percent_gold

        #Test out of bounds. find_point_indices warns and clamps; the
        #subsequent find_interpolation_percent linearly extrapolates,
        #returning percent = r (since the bracket is (1,2) on rvec=0:10).
        r = -0.45
        @test_logs (:warn, r"not found") (pair = WATT.find_point_indices(rvec, r))
        percent = WATT.find_interpolation_percent(rvec, pair, r)
        percent_gold = -0.45
        @test percent == percent_gold

        #Test something that doesn't have an interval of 1
        rvec = collect(0:10:100)
        r = 12
        pair = WATT.find_point_indices(rvec, r)
        percent = WATT.find_interpolation_percent(rvec, pair, r)
        percent_gold = 0.2
        @test percent == percent_gold

        ### Test create_interpolationpoints(assembly, blade)
        # Prep the ASD rotor and operating conditions 
        ofpath = "../data/openfast"
        adblade = of.read_adblade("sn5_adblade.dat", ofpath)
        edfile = of.read_edfile("sn5_EDfile.dat", ofpath)
        bdfile = of.read_bdfile("sn5_BDfile.dat", ofpath)
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

        # indices correspond to which airfoil is used at which station
        af_idx = Int.(adblade["BlAFID"])

        # create airfoil array
        afs = aftypes[af_idx]

        chordvec = adblade["BlChord"]
        twistvec = adblade["BlTwist"]
        rhub = edfile["HubRad"]
        rvec = adblade["BlSpn"] .+ rhub
        rtip = rvec[end]
        hubht = 80.0
        n = length(rvec)

        rvec = collect(range(rvec[1], rvec[end], length=n))

        airfoils = StructArray{DS.Airfoil}(undef, n)
        xcp = Vector{Float64}(undef, n)
        for i = 1:n
            airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
        end

        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; rhub=rhub, rtip=rtip)
        assembly = of.make_assembly(edfile, bdfile, bdblade)

        ips = WATT.create_interpolationpoints(assembly, blade)

        rgx = [norm(assembly.points[i]) for i in eachindex(assembly.points)]

        ### Structural sanity: pairs are consecutive, in range, percent ∈ [0,1]
        for ip in ips #todo: I bet we could rewrite this test to be more condensed, rather than like 75 tests.
            @test 1 <= ip.pair[1] < ip.pair[2] <= length(rgx)
            @test ip.pair[2] == ip.pair[1] + 1
            @test 0.0 <= ip.percent <= 1.0
        end

        ### Identity check: interpolating the (undeflected) structural radii
        ### at ip.percent between ip.pair must recover the aero node radius.
        for i in eachindex(ips)
            r_lo = rgx[ips[i].pair[1]]
            r_hi = rgx[ips[i].pair[2]]
            r_interp = (1 - ips[i].percent) * r_lo + ips[i].percent * r_hi
            @test isapprox(r_interp, blade.r[i]; atol=1e-10)
        end

        ### Between-node check: for every aero node, derive the expected
        ### bracket from rgx and compare to what create_interpolationpoints found.
        for i in eachindex(ips)
            roi = blade.r[i]
            j = findlast(r -> r <= roi, rgx)
            j === nothing && continue            # aero node below first structural point
            j == length(rgx) && continue         # aero node at or beyond tip
            @test ips[i].pair == (j, j + 1)
            @test isapprox(ips[i].percent, (roi - rgx[j]) / (rgx[j+1] - rgx[j]); atol=1e-12)
        end

    end #End mesh generation tests

    @testset "Interpolation Functions" begin

        # Tests the per-aero-node interpolation helpers used by update_mesh!
        # every time step in run_sim!:
        #
        #   - interpolate_deflection(ip, assembly, state)  → SVector{3} translation
        #   - interpolate_angle(ip, assembly, state)       → SVector{3} Euler rotation
        #   - interpolate_velocity(ip, assembly, state)    → SVector{3} structural velocity
        #   - convert_velocities(blade, env, assembly, state, ips, t, idx)
        #         → SVector{3} velocity in aero frame, with rotational component removed
        #
        # We drive structural deflection via a direct GXBeam.time_domain_analysis
        # with a known follower load, then check that the interpolation at root
        # and tip aero nodes (which coincide with structural endpoints) recovers
        # the corresponding state.points[].u / .theta / .V exactly, and that
        # an interior node matches a hand-derived linear blend.

        ofpath = "../data/openfast"
        adblade = of.read_adblade("sn5_adblade.dat", ofpath)
        edfile = of.read_edfile("sn5_EDfile.dat", ofpath)
        bdfile = of.read_bdfile("sn5_BDfile.dat", ofpath)
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
        afs = aftypes[af_idx]

        chordvec = adblade["BlChord"]
        twistvec = adblade["BlTwist"]
        rhub = edfile["HubRad"]
        rvec = adblade["BlSpn"] .+ rhub
        rtip = rvec[end]
        n = length(rvec)

        rvec = collect(range(rvec[1], rvec[end], length=n))

        airfoils = StructArray{DS.Airfoil}(undef, n)
        xcp = Vector{Float64}(undef, n)
        for i = 1:n
            airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
        end

        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; rhub=rhub, rtip=rtip)
        assembly = of.make_assembly(edfile, bdfile, bdblade)

        ips = WATT.create_interpolationpoints(assembly, blade)

        distributed_loads = Dict{Int64, GXBeam.DistributedLoads{Float64}}()
        nelem = length(assembly.elements)
        for ielem in 1:nelem
            distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem;
                fz_follower = (s) -> 10)
        end

        # Cantilever root.
        prescribed_conditions = Dict(
            1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

        rpm = 10.0
        omega = rpm*2*pi/60
        env = environment(1, 1, 1, 1, omega, 0)
        angular_velocity = SVector(0.0, 0.0, -omega)

        azimuth = 90*pi/180
        g = 9.81
        gravity = SVector(g*sin(azimuth), -g*cos(azimuth), 0.0)

        tvec = 0:0.1:0.3

        system, history, converged = GXBeam.time_domain_analysis(assembly, tvec;
            prescribed_conditions,
            distributed_loads,
            angular_velocity,
            gravity)

        @test converged

        state = history[end]
        na = length(blade.airfoils)

        delta  = [WATT.interpolate_deflection(ips[i], assembly, state) for i in 1:na]
        dtheta = [WATT.interpolate_angle(ips[i], assembly, state)      for i in 1:na]
        V      = [WATT.interpolate_velocity(ips[i], assembly, state)   for i in 1:na]
        aerov  = [norm(WATT.convert_velocities(blade, env, assembly, state, ips, tvec[end], i)) for i in 1:na]

        theta_root = WATT.WMPtoangle(state.points[1].theta)
        theta_tip  = WATT.WMPtoangle(state.points[end].theta)

        ### Root: the first aero node coincides with the first structural
        ### point (cantilevered, zero deflection by prescribed_conditions).
        @test isapprox(blade.r[1], assembly.points[1][1])
        @test isapprox(delta[1], state.points[1].u)
        @test isapprox(dtheta[1], SVector(0.0, 0.0, 0.0); atol=1e-12)
        @test isapprox(dtheta[1][1], theta_root[1]; atol=1e-12)
        @test isapprox(V[1], state.points[1].V)

        ### Tip: aero node is ~0.0001 m off from the structural tip point.
        @test isapprox(blade.r[end], assembly.points[end][1], rtol=1e-4)
        @test isapprox(delta[end], state.points[end].u, rtol=1e-4)
        @test isapprox(dtheta[end][1], theta_tip[1], rtol=1e-4)
        @test isapprox(V[end], state.points[end].V, rtol=1e-4)

        ### Intermediate point: a hand-derived linear blend at idx=2 should
        ### match interpolate_velocity's output.
        idx = 2
        r = blade.r[idx]
        rgx = WATT.get_bladelength_vector(assembly)
        pair = WATT.find_point_indices(rgx, r)
        percent = WATT.find_interpolation_percent(rgx, pair, r)
        V_gold = (1 - percent) * state.points[pair[1]].V + percent * state.points[pair[2]].V
        @test isapprox(V[idx], V_gold)

        ### Aero-frame velocity magnitudes should be non-negative (norms).
        @test all(>=(0.0), aerov[2:end])
    end

end #End Mesh Tests