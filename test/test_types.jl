#=
Test the Blade constructors, functions and methods. 

=#

using Test, Statistics
using Rotors, FLOWMath, DynamicStallModels, OpenFASTsr

DS = DynamicStallModels
of = OpenFASTsr

localpath = @__DIR__
cd(localpath)


@testset "Blades" begin
    ### Prep the ASD rotor and operating conditions 
    ofpath = "../testing/OpenFAST_NREL5MW"
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

    # create airfoil array
    afs = aftypes[af_idx]

    chordvec = adblade["BlChord"]
    twistvec = adblade["BlTwist"]
    rhub = edfile["HubRad"]
    rvec = adblade["BlSpn"] .+ rhub
    hubht = 80.0
    n = length(rvec)

    rvec = collect(range(rvec[1], rvec[end], length=n))

    ### Test that rvec is linearly spaced (some of the tests depend on that). 
    dr = rvec[2]-rvec[1]
    drvec = [rvec[i]-rvec[i-1] for i = 2:n]
    @test any(i->isapprox(i,dr), drvec)

    
    airfoils = Vector{DS.Airfoil}(undef, n)
    for i = 1:n
        airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
    end

    @testset "Constructors" begin

        ### Test blade with no curvature
        blade = Blade(rvec, twistvec, airfoils)
        @test isapprox(mean(blade.rx), 0.0)
        @test isapprox(mean(blade.ry), 0.0)
        @test mean(blade.rz)>0.0




        ### Test blade with precone
        precone=2*pi/180
        blade = Blade(rvec, twistvec, airfoils; precone)
        @test mean(blade.rx)>0
        @test mean(blade.ry)==0

        preconeflag = true
        for i in 2:n
            if blade.rx[i]<=blade.rx[i-1]
                preconeflag = false
            end
        end
        @test preconeflag
        # @test isapprox(blade.rx[end], (blade.rtip-blade.rhub)*sin(precone))
        @test isapprox(blade.rx[end], (blade.rtip)*sin(precone))




        ### Test blade with single curve value
        curve = 2*pi/180
        blade = Blade(rvec, twistvec, airfoils; curve)

        @test mean(blade.rx)>0
        @test mean(blade.ry)==0

        curveflag = true
        for i in 2:n
            if blade.rx[i]<=blade.rx[i-1]
                curveflag = false
            end
        end
        @test curveflag
        # @test isapprox(blade.rx[end], (blade.rtip-blade.rhub)*sin(curve))
        @test isapprox(blade.rx[end], (blade.rtip)*sin(curve))





        ### Test blade with increasing curve
        curve = collect(range(0, 10*pi/180, length=n))
        blade = Blade(rvec, twistvec, airfoils; curve)

        @test mean(blade.rx)>0
        @test mean(blade.ry)==0

        #Test that the dx is indeed increasing. 
        curveflag = true
        dx_old = [0.0]
        for i in 3:n
            dx = blade.rx[i] - blade.rx[i-1]
            if dx<dx_old[1]
                curveflag = false
            end
            dx_old[1] = dx
        end
        @test curveflag #Note: Test depends on linearly spaced rvec





        ### Test blade with single sweep
        sweep = 2*pi/180
        blade = Blade(rvec, twistvec, airfoils; sweep)

        @test mean(blade.ry)>0
        @test mean(blade.rx)==0


        sweepflag = true
        for i in 2:n
            if blade.ry[i]<=blade.ry[i-1]
                sweepflag = false
            end
        end
        @test sweepflag
        @test isapprox(blade.ry[end], (blade.rtip)*sin(sweep))





        ### Test blade with increasing sweep
        sweep = collect(range(0, 10*pi/180, length=n))
        blade = Blade(rvec, twistvec, airfoils; sweep)

        @test mean(blade.ry)>0  
        @test mean(blade.rx)==0
        bladelength = @. sqrt(blade.rx^2 + blade.ry^2 + blade.rz^2)



        #Test that the dy is indeed increasing. 
        sweepflag = true
        dy_old = [0.0]
        for i in 3:n
            dy = abs(blade.ry[i] - blade.ry[i-1])
            if dy<dy_old[1]
                sweepflag = false
            end
            dy_old[1] = dy
        end
        @test sweepflag #Note: Test depends on linearly spaced rvec

        # @show curve
        # @show blade.rx
        # @show blade.ry
        # @show blade.rz

        # using Plots
        # curve = collect(range(0, 10*pi/180, length=n))
        # sweep = collect(range(0, 10*pi/180, length=n))
        # precone = 2*pi/180

        # blade1 = Blade(rvec, twistvec, airfoils; curve)
        # blade2 = Blade(rvec, twistvec, airfoils; precone)
        # blade3 = Blade(rvec, twistvec, airfoils; sweep)

        # plt = plot(rvec, blade1.rx, lab="x", leg=:topleft)
        # # plot!(rvec, blade1.ry, lab="y")
        # plot!(rvec, blade1.rz, lab="z")
        # plot!(rvec, blade2.rx, lab="x", linestyle=:dash)
        # # plot!(rvec, blade2.ry, lab="y", linestyle=:dash)
        # plot!(rvec, blade2.rz, lab="z", linestyle=:dash)
        # plot!(rvec, blade3.rx, lab="x", linestyle=:dash)
        # # plot!(rvec, blade2.ry, lab="y", linestyle=:dash)
        # plot!(rvec, blade3.rz, lab="z", linestyle=:dash)
        # display(plt)
    end
end #End Blade testset