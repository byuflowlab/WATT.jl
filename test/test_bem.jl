#=
Tests to test the four different BEM implementations: 1) single residual (inflow angle), 2) single fixed point iteration, 3) dual residual (induction factors), 4) dual fixed point iteration. 

Adam Cardoza
=#
using Test, WATT, DelimitedFiles, FLOWMath, CCBlade, OpenFASTTools, DynamicStallModels
# using Test, DelimitedFiles, FLOWMath, CCBlade

of = OpenFASTTools
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)

# include("/Users/adamcardoza/.julia/dev/WATT/src/bem.jl")

errfun(x, xt) = 100*(x - xt)/xt

# @testset "Blade Element Momentum Theory" begin

#     @testset "CCBlade" begin
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
        rtip = rvec[end]
        precone = 0.0
        hubht = 80.0
        n = length(rvec)


        airfoils = StructArray{DS.Airfoil}(undef, n)
        xcp = Vector{Float64}(undef, n)
        for i = 1:n
            airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
        end 

        blade = Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils)

        B = 3
        hubht = 80.0
        turbine = true
        rotor = WATT.Rotor(B, hubht, turbine)

        vinf = 10.0
        tsr = 7.55
        rotorR = rtip*cos(precone)
        omega = vinf*tsr/rotorR

        rho = 1.225
        mu = 1.464e-5 
        a = 343.0
        shearexp = 0.0
        env = environment(rho, mu, a, vinf, omega, shearexp)

        pitch = 0.0
        Vx = vinf

        xv = zeros(11)   # scratch buffer required by solve_BEM!

        for idx = 1:n
            Vy = omega*rvec[idx]

            rotorout = WATT.solve_BEM!(rotor, blade, env, idx, Vx, Vy, pitch, xv; npts=10)

            polar = blade.airfoils[idx].polar
            af = CCBlade.AlphaAF(polar[:,1], polar[:,2], polar[:,3])
            section = Section(rvec[idx], chordvec[idx], twistvec[idx], af)
            rotor_cc = CCBlade.Rotor(blade.rhub, blade.rtip, B; turbine=true, tip=nothing)
            op = CCBlade.OperatingPoint(Vx, Vy, rho, pitch, mu, a)
            ccout = solve(rotor_cc, section, op)

            ### Wrapper-fidelity: same concrete return type as CCBlade.
            @test isa(rotorout, typeof(ccout))

            ### Numerical agreement, field-by-field, against a direct CCBlade call.
            @testset "idx=$idx field=$f" for f in fieldnames(typeof(ccout))
                @test isapprox(getfield(rotorout, f), getfield(ccout, f); rtol=1e-8, atol=1e-10)
            end

            ### Physical sanity: induction factor should not exceed 1.
            ### (a > 1 means BEMT has entered the turbulent-wake regime where standard BEMT breaks down.)
            @test ccout.a < 1
        end


#     end #End test CCBlade
# end #End test BEMT

nothing

#=
Things I inserted into CCBlade that I might need: 

- Just after the zero Outputs initializer
Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G) = Outputs(promote(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)...)

=#