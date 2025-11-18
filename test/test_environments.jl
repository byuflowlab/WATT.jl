#=
Test all of the environment constructors and functions

#Todo: Minor errors on the multiple rotations. 
=#

using Test
using WATT, CCBlade, OpenFASTsr, DynamicStallModels, GXBeam, StaticArrays

DS = DynamicStallModels
of = OpenFASTsr

localpath = @__DIR__
cd(localpath)



@testset "Environments" begin
    ### Prep the ASD rotor and operating conditions 
    ofpath = "../testing/OpenFAST_NREL5MW"
    addriver = of.read_addriver("NREL5MW_ADdriver.dvr", ofpath)
    adblade = of.read_adblade("NREL5MW_adblade.dat", ofpath)
    edfile = of.read_edfile("NREL5MW_EDfile.dat", ofpath)

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

    # @show af_idx[10]

    # create airfoil array
    afs = aftypes[af_idx]

    chordvec = adblade["BlChord"]
    twistvec = adblade["BlTwist"]
    rhub = edfile["HubRad"]
    rvec = adblade["BlSpn"] .+ rhub
    hubht = 80.0
    n = length(rvec)


    airfoils = Vector{DS.Airfoil}(undef, n)
    for i = 1:n
        airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
    end 

    vinf = addriver["HWndSpeed_mat"][1] #10.0
    rpm = addriver["RotSpd_mat"][1]
    omega = rpm*(2*pi)/60 #vinf*tsr/rotorR

    rho = addriver["FldDens"] #1.225
    mu = addriver["KinVisc"] #1.464e-5 #18.13e-6
    a = addriver["SpdSound"] #343.0
    shearexp = addriver["PLExp"][1] #0.0

    B = 3
    hubht = 80.0
    turbine = true

    @testset "Simple Environments" begin

        ### No tilt, yaw, precone
        pitch = precone = yaw = tilt = 0.0
        dx = dy = dz = 0.0
        t = 0.0

        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        @test isa(env, WATT.Environment) #Test the typing. 

        idx = 10
        r = rvec[idx]
        # @show r

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("Nothing: ")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)
        # println("")





        ### Test at a different azimuthal angle
        t = 4.5
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("Azimuthal")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)





        ### Test with precone
        t = 0.0
        precone = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("precone")
        # println("My code:")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)

        


        ### Test with precone at a different azimuthal angle.
        t = 3.2
        precone = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("precone azimuth:")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)





        ### Test with tilt
        t = 0.0
        precone = 0.0
        tilt = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("Tilted flow:")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)



        ### Test with tilt at a different azimuth
        t = 3.2
        precone = 0.0
        tilt = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        # println("")
        # @show (azimuth)*180/pi

        rotor = WATT.Rotor(B, hubht, turbine; tilt)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("tilt with azimuth")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)






        ### Test with yaw
        t = 0.0
        precone = 0.0
        tilt = 0.0
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("yaw")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)




        ### Test with yaw at a different azimuth
        t = 5.6
        precone = 0.0
        tilt = 0.0
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("yaw with azimuth")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)





        # ### Test with sweep #Todo:




        ### Test with tilt, and yaw
        t = 0.0
        precone = 0.0
        tilt = 34*pi/180
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("Tilt and yaw")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        # @show op.Vx, op.Vy

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy) #rtol=0.02 #Todo: This increased error occurs when I rotate the freestream velocities into the aerodynamic frame... They appear to be off by 2 mm/s in the z direction (comparing z to y, where y is CCBlade's equivalent velocity)




        ## Test tilt, yaw, and a different azimuth
        t = 2.2
        precone = 0.0
        tilt = 34*pi/180
        yaw = 25*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("Tilt, yaw, and azimuth")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        # @show op.Vx, op.Vy

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy) #, rtol=0.02 #Todo: This increased error occurs when I rotate the freestream velocities into the aerodynamic frame... They appear to be off by 2 mm/s in the z direction (comparing z to y, where y is CCBlade's equivalent velocity)




        ### Test with tilt, yaw, and precone
        t = 0.0
        precone = 25*pi/180
        tilt = 38*pi/180
        yaw = 49*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("tilt, yaw, and precone")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo
        # # println("")

        @test isapprox(Vx, op.Vx) #, rtol=0.0001
        @test isapprox(Vy, op.Vy) #, rtol=0.02





        ### Test with tilt, yaw, and precone
        t = pi/2
        precone = 25*pi/180
        tilt = 38*pi/180
        yaw = 49*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("tilt, yaw, precone, and azimuth")
        # println("My code: ")
        
        Vx, Vy = WATT.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = WATT.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo
        # # println("")

        @test isapprox(Vx, op.Vx) #, rtol=0.0001
        @test isapprox(Vy, op.Vy) #, rtol=0.02



        #----------------------------------------#


        ### Test aerostructural velocities
        precone = 0.0 # 3*pi/180
        tilt = 0.0 #4*pi/180
        yaw = 0.0 #5*pi/180

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

        rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)
        
        # @show tilt, yaw, precone

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

        
        angular_velocity = SVector(0.0, 0.0, -omega)

        t = 0.2
        azimuth = omega*t
        g = 9.81
        gravity = SVector(g*sin(azimuth), -g*cos(azimuth), 0.0)

        tvec = 0:0.1:t

        system, history, converged = GXBeam.time_domain_analysis(assembly, tvec;
            prescribed_conditions,
            distributed_loads,
            angular_velocity,
            gravity,
            steady=false) 

        state = history[end]

        na = length(blade.airfoils)

        delta = [WATT.interpolate_deflection(ips[i], assembly, state) for i in 1:na]
        aerov = [WATT.convert_velocities(blade, env, assembly, state, ips, t, i) for i in 1:na]

        # @show state.points[1].u
        # @show state.points[1].V
        # @show env.RS(0)
        # @show env.RS(t)
        # @show aerov[1]

        Vxvec = zeros(na)
        Vyvec = zeros(na)

        for i = 1:na
            Vxvec[i], Vyvec[i] = WATT.get_aerostructural_velocities(rotor, blade, env, t, i, azimuth, delta[i], aerov[i])
        end

        # @show blade.r

        increasingflag = true
        for i = 2:na
            if Vyvec[i-1]>=Vyvec[i]
                increasingflag = false
            end
        end
        @test increasingflag


        ### #Todo: I don't know what other tests to do here because I don't know what the values should be... they look pretty good. I mean, I guess I could pull velocities out of OpenFAST and compare those. 



    end #End testing simple environments
end #End testing environments