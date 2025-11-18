#=
Test the meshing and mesh transfer functions. 

=#

using Test
using WATT, GXBeam, OpenFASTsr, DynamicStallModels
using LinearAlgebra, StaticArrays

of = OpenFASTsr
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

        #Out of bounds tests
        r = -1
        pair = WATT.find_point_indices(rvec, r)
        pair_gold = (1,2)
        @test pair==pair_gold

        r = 4*pi
        pair = WATT.find_point_indices(rvec, r)
        pair_gold = (10,11)
        @test pair==pair_gold

        ### Test find_interpolation_percent
        #Test a random spot
        r = 4.75
        pair = WATT.find_point_indices(rvec, r)
        percent = WATT.find_interpolation_percent(rvec, pair, r)
        percent_gold = 0.75
        @test percent==percent_gold

        #Test out of bounds
        r = -0.45
        pair = WATT.find_point_indices(rvec, r)
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
        ofpath = "../testing/OpenFAST_NREL5MW"
        adblade = of.read_adblade("NREL5MW_adblade.dat", ofpath)
        edfile = of.read_edfile("NREL5MW_EDfile.dat", ofpath)
        bdfile = of.read_bdfile("NREL5MW_BDfile.dat", ofpath)
        bdblade = of.read_bdblade("NREL5MW_BDblade.dat", ofpath)

        aftypes = Array{of.AirfoilInput}(undef, 8)
        aftypes[1] = of.read_airfoilinput(ofpath*"/Airfoils/Cylinder1.dat") 
        aftypes[2] = of.read_airfoilinput(ofpath*"/Airfoils/Cylinder2.dat") 
        aftypes[3] = of.read_airfoilinput(ofpath*"/Airfoils/DU40_A17.dat") 
        aftypes[4] = of.read_airfoilinput(ofpath*"/Airfoils/DU35_A17.dat") 
        aftypes[5] = of.read_airfoilinput(ofpath*"/Airfoils/DU30_A17.dat") 
        aftypes[6] = of.read_airfoilinput(ofpath*"/Airfoils/DU25_A17.dat") 
        aftypes[7] = of.read_airfoilinput(ofpath*"/Airfoils/DU21_A17.dat") 
        aftypes[8] = of.read_airfoilinput(ofpath*"/Airfoils/NACA64_A17.dat") 

        # indices correspond to which airfoil is used at which station
        af_idx = Int.(adblade["BlAFID"])

        # create airfoil array
        afs = aftypes[af_idx]

        chordvec = adblade["BlChord"]
        twistvec = adblade["BlTwist"]
        rhub = edfile["HubRad"]
        rvec = adblade["BlSpn"] .+ rhub
        hubht = 80.0
        n = length(rvec)

        rvec = collect(range(rvec[1], rvec[end], length=n))

        airfoils = Vector{DS.Airfoil}(undef, n)
        for i = 1:n
            airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
        end

        blade = Blade(rvec, twistvec, airfoils)
        assembly = of.make_assembly(edfile, bdfile, bdblade)

        ips = WATT.create_interpolationpoints(assembly, blade)

        rgx = [norm(assembly.points[i]) for i in eachindex(assembly.points)]
        
        idx = 2
        roi = blade.r[idx]
        pair_gold = (5, 6)
        a = roi - rgx[5]
        L = rgx[6]-rgx[5]
        percent_gold = a/L

        @test ips[idx].pair == pair_gold
        @test ips[idx].percent == percent_gold

    end #End mesh generation tests

    @testset "Interpolation Functions" begin

        # Prep the ASD rotor and operating conditions 
        ofpath = "../testing/OpenFAST_NREL5MW"
        adblade = of.read_adblade("NREL5MW_adblade.dat", ofpath)
        edfile = of.read_edfile("NREL5MW_EDfile.dat", ofpath)
        bdfile = of.read_bdfile("NREL5MW_BDfile.dat", ofpath)
        bdblade = of.read_bdblade("NREL5MW_BDblade.dat", ofpath)

        aftypes = Array{of.AirfoilInput}(undef, 8)
        aftypes[1] = of.read_airfoilinput(ofpath*"/Airfoils/Cylinder1.dat") 
        aftypes[2] = of.read_airfoilinput(ofpath*"/Airfoils/Cylinder2.dat") 
        aftypes[3] = of.read_airfoilinput(ofpath*"/Airfoils/DU40_A17.dat") 
        aftypes[4] = of.read_airfoilinput(ofpath*"/Airfoils/DU35_A17.dat") 
        aftypes[5] = of.read_airfoilinput(ofpath*"/Airfoils/DU30_A17.dat") 
        aftypes[6] = of.read_airfoilinput(ofpath*"/Airfoils/DU25_A17.dat") 
        aftypes[7] = of.read_airfoilinput(ofpath*"/Airfoils/DU21_A17.dat") 
        aftypes[8] = of.read_airfoilinput(ofpath*"/Airfoils/NACA64_A17.dat") 

        # indices correspond to which airfoil is used at which station
        af_idx = Int.(adblade["BlAFID"])

        # create airfoil array
        afs = aftypes[af_idx]

        chordvec = adblade["BlChord"]
        twistvec = adblade["BlTwist"]
        rhub = edfile["HubRad"]
        rvec = adblade["BlSpn"] .+ rhub
        hubht = 80.0
        n = length(rvec)

        rvec = collect(range(rvec[1], rvec[end], length=n))

        airfoils = Vector{DS.Airfoil}(undef, n)
        for i = 1:n
            airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
        end

        blade = Blade(rvec, twistvec, airfoils)
        assembly = of.make_assembly(edfile, bdfile, bdblade)

        ips = WATT.create_interpolationpoints(assembly, blade)

        distributed_loads = Dict{Int64, DistributedLoads{Float64}}()

        nelem = length(assembly.elements)

        for ielem in 1:nelem
            distributed_loads[ielem] = DistributedLoads(assembly, ielem;
                fz_follower = (s) -> 10)
        end

        # root section is fixed
        prescribed_conditions = Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) 

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
            gravity,
            steady=false) 

        state = history[end]

        na = length(blade.airfoils)

        delta = [WATT.interpolate_deflection(ips[i], assembly, state) for i in 1:na]
        dtheta = [WATT.interpolate_angle(ips[i], assembly, state) for i in 1:na]
        V = [WATT.interpolate_velocity(ips[i], assembly, state) for i in 1:na]
        aerov = [norm(WATT.convert_velocities(blade, env, assembly, state, ips, tvec[end], i)) for i in 1:na]


        theta_root = WATT.WMPtoangle(state.points[1].theta)
        theta_tip = WATT.WMPtoangle(state.points[end].theta)

        ## Check the endpoints where the nodes match up.
        #Root 
        @test isapprox(blade.r[1], assembly.points[1][1])
        @test isapprox(delta[1], state.points[1].u)
        @test isapprox(dtheta[1], 0.0)
        @test isapprox(dtheta[1], theta_root[1])
        @test isapprox(V[1], state.points[1].V)


        #Tip (There is a tolerance because the tip aero node is 0.0001 m off of the structural node).
        @test isapprox(blade.r[end], assembly.points[end][1], rtol=0.0001)
        @test isapprox(delta[end], state.points[end].u, rtol=0.0001)
        @test isapprox(dtheta[end], theta_tip[1], rtol=0.0001)
        @test isapprox(V[end], state.points[end].V, rtol=0.0001)

        
        ## Check an intermediate point 
        idx = 2
        r = blade.r[idx]

        rgx = WATT.get_bladelength_vector(assembly)
        pair = WATT.find_point_indices(rgx, r)
        percent = WATT.find_interpolation_percent(rgx, pair, r)

        V_gold = state.points[pair[1]].V*(1-percent) + percent*state.points[pair[2]].V

        @test V[2]==V_gold

        

        ## All point checks
        @test !any(i->i<0, aerov[2:end])




        ### Test that the the structural velocity goes to zero at the steady state solution.
        system, history, converged = GXBeam.time_domain_analysis(assembly, tvec;
            prescribed_conditions,
            distributed_loads,
            angular_velocity,
            gravity,
            steady=true) 

        state = history[end]
        aerov = [WATT.convert_velocities(blade, env, assembly, state, ips, tvec[end], i) for i in 1:na]


        steadyflag = true
        for i in 1:na
            if any(i->!isapprox(i, 0.0, atol=1e-12), aerov[i])
                steadyflag = false
            end
        end
        @test steadyflag
    end

end #End Mesh Tests