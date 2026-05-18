#=
Test all of the environment constructors and functions

#Todo: Minor errors on the multiple rotations. 
=#

using Test
using WATT, CCBlade, OpenFASTTools, DynamicStallModels, GXBeam, StaticArrays, StructArrays

DS = DynamicStallModels
of = OpenFASTTools

localpath = @__DIR__
cd(localpath)



@testset "Environments" begin
    ### Prep the ASD rotor and operating conditions 
    ofpath = "../data/openfast"
    adblade = of.read_adblade("sn5_adblade.dat", ofpath)
    edfile = of.read_edfile("sn5_EDfile.dat", ofpath)

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
    hubht = 80.0
    n = length(rvec)


    airfoils = StructArray{DS.Airfoil}(undef, n)
    xcp = Vector{Float64}(undef, n)
    for i = 1:n
        airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    end 

    vinf = 10.0
    tsr = 7.55
    rotorR = rvec[end]
    omega = vinf*tsr/rotorR

    rho = 1.225
    mu = 1.464e-5 
    a = 343.0
    shearexp = 0.0

    B = 3
    hubht = 80.0
    turbine = true

    @testset "Simple Environments" begin

        # Each block below tests WATT.get_aero_velocities against CCBlade's
        # reference windturbine_op for one geometric configuration. The aero
        # velocities (Vx, Vy) feeding the BEM must match a known-correct
        # implementation across every combination of precone / tilt / yaw /
        # azimuth — a coordinate-transform regression suite.
        #
        # Reference: CCBlade.windturbine_op assembles (Vx, Vy) from the wind
        # vector by transforming the global wind frame into the local airfoil
        # frame. WATT's get_aero_velocities must produce identical values.
        #
        # All blocks evaluate at idx=10 (midspan-ish). The common setup —
        # blade chord/twist tables, airfoils, environmental constants — is
        # built once above this testset.

        idx = 10
        r = rvec[idx]
        pitch = 0.0

        ### --- Case 1: undeflected, axisymmetric, t=0 ---
        ### No precone, tilt, or yaw. Reduces to pure axial inflow / rigid
        ### rotation — sanity check on the basic transformation chain.
        precone = yaw = tilt = 0.0
        t = 0.0
        azimuth = omega*t   # constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils)
        env = WATT.environment(rho, mu, a, vinf, omega, shearexp)

        @test isa(env, WATT.Environment)   # constructor returns the right abstract subtype

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)


        ### --- Case 2: undeflected, t=4.5 s ---
        ### Same geometry, nonzero azimuth. Verifies the rotation transform
        ### is azimuth-invariant for the axisymmetric case.
        t = 4.5
        azimuth = omega*t

        op = WATT.windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)


        ### --- Case 3: precone only, t=0 ---
        ### Adds 25° precone. Vx should drop by cos(precone) since the rotor
        ### disk tilts away from the wind. Vy picks up a wind component too.
        t = 0.0
        precone = 25*pi/180
        azimuth = omega*t

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = WATT.environment(rho, mu, a, vinf, omega, shearexp)

        op = WATT.windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)


        ### --- Case 4: precone + azimuth, t=3.2 s ---
        ### Precone shouldn't depend on azimuth (rotational symmetry about
        ### the shaft axis), so Vx/Vy should match Case 3 at this azimuth.
        t = 3.2
        precone = 25*pi/180
        azimuth = omega*t

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = WATT.environment(rho, mu, a, vinf, omega, shearexp)

        op = WATT.windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)





        ### --- Case 5: tilt only, t=0 ---
        ### 25° shaft tilt. Rotor disk tips relative to wind; both Vx and Vy
        ### pick up azimuth-dependent projections of the wind vector.
        t = 0.0
        precone = 0.0
        tilt = 25*pi/180
        azimuth = omega*t

        rotor = WATT.Rotor(B, hubht, turbine; tilt)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = WATT.environment(rho, mu, a, vinf, omega, shearexp)

        op = CCBlade.windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)


        ### --- Case 6: tilt + azimuth, t=3.2 s ---
        ### Same tilt at a different blade position. Tilted disk + nonzero
        ### azimuth is the regime where most coordinate-transform bugs
        ### historically lived (sign-flip in transform_BC_HR etc).
        t = 3.2
        precone = 0.0
        tilt = 25*pi/180
        azimuth = omega*t

        rotor = WATT.Rotor(B, hubht, turbine; tilt)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)






        ### --- Case 7: yaw only, t=0 ---
        ### 25° nacelle yaw. Wind hits the disk at an angle; expect skewed
        ### inflow, with Vx reduced and a cross-flow component appearing.
        t = 0.0
        precone = 0.0
        tilt = 0.0
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)




        ### --- Case 8: yaw + azimuth, t=5.6 s ---
        t = 5.6
        precone = 0.0
        tilt = 0.0
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)





        #Todo: add a sweep test case once sweep is supported.




        ### --- Case 9: tilt + yaw, t=0 ---
        ### Two non-axisymmetric rotations composed. Catches order-of-rotation
        ### bugs that pure-tilt or pure-yaw cases miss.
        t = 0.0
        precone = 0.0
        tilt = 34*pi/180
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)



        ### --- Case 10: tilt + yaw + azimuth, t=2.2 s ---
        t = 2.2
        precone = 0.0
        tilt = 34*pi/180
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)



        ### --- Case 11: tilt + yaw + precone, t=0 ---
        ### Full rotation chain. If the order of transforms (precone → hub →
        ### tilt → yaw → global) is wrong anywhere, this case fails.
        t = 0.0
        precone = 25*pi/180
        tilt = 38*pi/180
        yaw = 49*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)




        ### --- Case 12: tilt + yaw + precone + azimuth, t=π/2 ---
        ### Same full chain at a 90° azimuth — picks a specific blade
        ### orientation (horizontal) where many transforms reduce to simple
        ### swaps; useful for catching subtle sign errors.
        t = pi/2
        precone = 25*pi/180
        tilt = 38*pi/180
        yaw = 49*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)

    end #End testing simple environments


    @testset "Aerostructural Velocities" begin

        # Tests WATT.get_aerostructural_velocities, which maps structural
        # deflections + velocities at an aero node into the airfoil-frame
        # velocities (Vx, Vy) consumed by the BEM.
        #
        # Design: drive structural deflection by applying a known external
        # load through GXBeam directly (constant 10 N/m radial follower load
        # on every element, plus gravity at a fixed azimuth). This isolates
        # the velocity transform from WATT's aero-structural simulation loop
        # — the full-pipeline regression lives in simple_NREL5MW.jl.
        #
        # Then interpolate the structural state onto each aero node using the
        # same helpers update_mesh! uses (mesh.jl:202–206), and assert:
        #   - GXBeam solve converged
        #   - Vy increases monotonically with r (dominated by Ω·r)
        #   - No NaN / no Inf
        #   - At the prescribed (fixed) root, get_aerostructural_velocities
        #     reduces exactly to get_aero_velocities

        bdfile  = of.read_bdfile("sn5_BDfile.dat", ofpath)
        bdblade = of.read_bdblade("sn5_BDblade.dat", ofpath)

        precone = 0.0
        tilt = 0.0
        yaw = 0.0
        pitch = 0.0

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; precone)
        env   = environment(rho, mu, a, vinf, omega, shearexp)

        assembly = of.make_assembly(edfile, bdfile, bdblade)
        ips      = WATT.create_interpolationpoints(assembly, blade)

        # Constant radial follower load on each element to produce a
        # non-zero structural state. Magnitude is arbitrary.
        distributed_loads = Dict{Int, GXBeam.DistributedLoads{Float64}}()
        for ielem in 1:length(assembly.elements)
            distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem;
                fz_follower = s -> 10.0)
        end

        # Cantilever root.
        prescribed_conditions = Dict(
            1 => GXBeam.PrescribedConditions(
                ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

        t = 0.2
        azimuth = omega*t
        g = 9.81
        gravity = SVector(g*sin(azimuth), -g*cos(azimuth), 0.0)
        angular_velocity = SVector(0.0, 0.0, -omega)

        tvec = 0:0.1:t

        system, history, converged = GXBeam.time_domain_analysis(
            assembly, tvec;
            prescribed_conditions, distributed_loads,
            angular_velocity, gravity)

        @test converged
        state = history[end]
        na    = length(blade.airfoils)

        # Interpolate structural quantities onto each aero node.
        delta       = [WATT.interpolate_deflection(ips[i], assembly, state)            for i in 1:na]
        delta_theta = [WATT.interpolate_angle(ips[i], assembly, state)                 for i in 1:na]
        aerov       = [WATT.convert_velocities(blade, env, assembly, state, ips, t, i) for i in 1:na]

        Vxvec = zeros(na)
        Vyvec = zeros(na)
        for i in 1:na
            Vxvec[i], Vyvec[i] = WATT.get_aerostructural_velocities(
                rotor, blade, env, t, i, azimuth,
                delta[i], delta_theta[i], aerov[i])
        end

        # Physical sanity: Vy is dominated by Ω·r so it must increase with r.
        @test issorted(Vyvec)

        # No NaN / no Inf
        @test all(isfinite, Vxvec)
        @test all(isfinite, Vyvec)

        # Root identity: at the prescribed (fixed) root, structural
        # quantities are zero, so get_aerostructural_velocities must reduce
        # exactly to get_aero_velocities.
        Vx_aero, Vy_aero = WATT.get_aero_velocities(rotor, blade, env, t, 1, azimuth)
        @test isapprox(Vxvec[1], Vx_aero; atol=1e-10)
        @test isapprox(Vyvec[1], Vy_aero; atol=1e-10)

    end #End testing aerostructural velocities

end #End testing environments