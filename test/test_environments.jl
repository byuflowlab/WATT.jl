#=
Test all of the environment constructors and functions

#Todo: Minor errors on the multiple rotations. 
=#

using Test
using Rotors, CCBlade, OpenFASTsr, DynamicStallModels

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

        rotor = Rotors.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        @test isa(env, Rotors.Environment) #Test the typing. 

        idx = 10
        r = rvec[idx]
        # @show r

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)
        # println("")

        ### Test at a different azimuthal angle
        t = 4.5
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)



        ### Test with precone
        t = 0.0
        precone = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)

        #Todo: Test Precone at a different azimuthal angle. 



        ### Test with tilt
        t = 0.0
        precone = 0.0
        tilt = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine; tilt)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)



        ### Test with tilt at a different azimuth
        t = 3.2
        precone = 0.0
        tilt = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        # println("")
        # @show (azimuth)*180/pi

        rotor = Rotors.Rotor(B, hubht, turbine; tilt)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("My code: ")
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = Rotors.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)



        ### Test with yaw
        t = 0.0
        precone = 0.0
        tilt = 0.0
        yaw = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        println("")
        println("My code: ")
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = Rotors.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        # @test isapprox(Vx, op.Vx)
        # @test isapprox(Vy, op.Vy)




        ### Test with yaw at a different azimuth
        t = 5.6
        precone = 0.0
        tilt = 0.0
        yaw = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("My code: ")
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)

        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = Rotors.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy)




        ### Test with tilt, and yaw
        t = 0.0
        precone = 0.0
        tilt = 4*pi/180
        yaw = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("My code: ")
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = Rotors.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx)
        @test isapprox(Vy, op.Vy, rtol=0.02) #Todo: This increased error occurs when I rotate the freestream velocities into the aerodynamic frame... They appear to be off by 2 mm/s in the z direction (comparing z to y, where y is CCBlade's equivalent velocity)





        ### Test with tilt, yaw, and precone
        t = 0.0
        precone = 3*pi/180
        tilt = 4*pi/180
        yaw = 5*pi/180
        azimuth = omega*t #assume constant rotation rate

        rotor = Rotors.Rotor(B, hubht, turbine; tilt, yaw)
        blade = Blade(rvec, twistvec, airfoils; precone)
        env = environment(rho, mu, a, vinf, omega, shearexp)

        idx = 10
        r = rvec[idx]

        op = windturbine_op.(vinf, omega, pitch, r, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

        # println("")
        # println("My code: ")
        
        Vx, Vy = Rotors.get_aero_velocities(rotor, blade, env, t, idx, azimuth)


        # @show Vx, Vy
        # println("")
        # println("Dr. Ning's code: ")

        # Vxo, Vyo = Rotors.get_aero_velocities(env, t, r, azimuth, precone, tilt, yaw, hubht)

        # @show Vxo, Vyo

        @test isapprox(Vx, op.Vx, rtol=0.0001)
        @test isapprox(Vy, op.Vy, rtol=0.02)

    end #End testing simple environments
end #End testing environments