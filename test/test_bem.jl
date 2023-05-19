#=
Tests to test the four different BEM implementations: 1) single residual (inflow angle), 2) single fixed point iteration, 3) dual residual (induction factors), 4) dual fixed point iteration. 

Adam Cardoza
=#
using Revise
using Test, Rotors, DelimitedFiles, FLOWMath, CCBlade, OpenFASTsr, DynamicStallModels
# using Test, DelimitedFiles, FLOWMath, CCBlade

of = OpenFASTsr
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)

# include("/Users/adamcardoza/.julia/dev/Rotors/src/bem.jl")

errfun(x, xt) = 100*(x - xt)/xt

@testset "Blade Element Momentum Theory" begin

    ### Adam Cardoza 3/17/23
    @testset "Dual fixed point iteration" begin

        U = 10.0
        rpm = 50.0
        omega = rpm*2*pi/60
        B = 3
        c = 1.0
        r = 42.431438
        rtip = 63.0
        rhub = 1.5
        twist = 0.0
        pitch = 0.0
        beta = twist + pitch
        polar = readdlm("../data/polars/NACA64_A17.dat", skipstart=3)
        alpha = polar[:,1].*(pi/180)
        clfit = Akima(alpha, polar[:,2])
        cdfit = Akima(alpha, polar[:,3])
        maxiters = 25
        tolerance = 2e-12

        phi0 = 1.301071424322167E-004
        alpha0 = phi0 - beta #Angle of attack, in the paragraph before EQ 20
        Cl = clfit(alpha0) #Static section coefficient of lift
        Cd = cdfit(alpha0) #Static section coefficient of drag
        cphi = cos(phi0)
        sphi = sin(phi0)

        F = 1.0
        lambda = omega*r/U

        #cn:   0.442850496085517     
        #ct:  -5.142382031116966E-003
        cx = Cl*cphi + Cd*sphi #x
        cy = Cl*sphi - Cd*cphi #x
        a, ap, R = Rotors.inductionfactors(phi0, U, omega, c, r, B, cx, cy, F) 

        # OpenFAST's first iteration values
        a_gold = 0.997398000367958  
        ap_gold = -0.100062053901359 


        #I'm not sure that this function is going to work. 
        # phi, a, ap, converged = solve_inflowangle_aerodyn(U, omega, r, rtip, rhub, B, c, beta, clfit, cdfit, maxiters, tolerance)
        # @show phi, a, ap, converged

        rho = 1.225
        mu = 1.4639e-5
        asound = 335.0

        af = AlphaAF("../data/polars/NACA64_A17.dat", radians=false)
        section = Section(r, c, beta, af)
        rotor = CCBlade.Rotor(rhub, rtip, B; turbine=true, tip=nothing)
        op = CCBlade.OperatingPoint(U, omega*r, rho, pitch, mu, asound)
        out = solve(rotor, section, op)

        function residual(phir)
            alphar = phir - beta
            clr = clfit(alphar)
            cdr = cdfit(alphar)
            cxr = clr*cphi + cdr*sphi #x
            cyr = clr*sphi - cdr*cphi #x
            _, _, Ri = Rotors.inductionfactors(phir, U, omega, c, r, B, cxr, cyr, F)
            return Ri
        end

        # package up variables and parameters for residual
        xv = [section.r, section.chord, section.theta, rotor.Rhub, rotor.Rtip, op.Vx, op.Vy, op.rho, op.pitch, op.mu, op.asound]
        pv = (section.af, rotor.B, rotor.turbine, rotor.re, rotor.mach, rotor.rotation, rotor.tip)
        
        function resid(phir)
            Ri, _ = CCBlade.residual_and_outputs(phir, xv, pv)
            return Ri
        end

        phistar, info = brent(residual, eps(), pi/2 - eps(); maxiter=5000)

        phi2, info2 = Rotors.sub_brent(residual, eps(), pi/2 - eps(), 0.0; maxiter=5000)
        phi3, info3 = Rotors.sub_brent(resid, eps(), pi/2 - eps(), 0.0; maxiter=5000)

        phi4, info4 = brent(resid, eps(), pi/2 - eps(); maxiter=5000, atol=2e-15, rtol=1e-10) #Tightening these tolerances, or loosening them makes it difficult to tie down if it matches the OpenFAST solutions. 



        

        ### OpenFAST's Converged values
        # phi_gold = 0.00742743941*pi/180 #Tightening the convergence criteria by altering the source code you get 0.0001304168554254917 (radians)... 
        # a_gold = 0.997409395 
        # ap_gold = -0.100506008
        ### OpenFAST's converged values after artificially tightening the convergence criteria (by altering the source code). 
        phi_gold = 0.0001304168554254917
        a = 0.997391823156366     
        ap = -9.984575201549731E-002

        println("")
        # @show phistar #AD residual, flowmath brent
        # @show phi2 #AD residual, sub_brent
        # @show phi3 #CC residual, sub_brent
        # @show phi4 #CC Residual, flowmath brent
        # # @show phi, a, ap, converged
        # @show out.phi, out.a, out.ap #CC residual, flowmath brent
        # @show phi_gold, a_gold, ap_gold #OpenFAST

        phiAfbC = errfun(phistar, out.phi)
        phiAsbC = errfun(phi2, out.phi)
        phiCsbC = errfun(phi3, out.phi)
        phiOFC = errfun(phi_gold, out.phi)

        phiAfbO = errfun(phistar, phi_gold)
        phiAsbO = errfun(phi2, phi_gold)
        phiCsbO = errfun(phi3, phi_gold)
        phiCOF = errfun(out.phi, phi_gold)

        # println("")
        # println("")

        # println("AD fb: ", phistar, ", ", phiAfbC, ", ", phiAfbO)
        # println("AD sb: ", phi2, ", ", phiAsbC, ", ", phiAsbO)
        # println("CC sb: ", phi3, ", ", phiCsbC, ", ", phiCsbO)
        # println("OF   : ", phi_gold, ", ", phiOFC, ", ", 0.0)
        # println("CCBla: ", out.phi, ", ", 0.0, ", ", phiCOF)
        

        # @show phierr, aerr, aperr #Phi converges within a half of a percent, but a and ap are like more than 10% off. Which is less fun. So maybe I should be comparing to OpenFAST, not CCBlade... because they are two different methods and I'm getting two different results. 

    end #End Dual fixed point iteration

    @testset "CCBlade" begin
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

        blade = Blade(rvec, twistvec, airfoils)

        B = 3
        hubht = 80.0
        turbine = true
        rotor = Rotors.Rotor(B, hubht, turbine)

        vinf = addriver["HWndSpeed_mat"][1] #10.0
        # tsr = 7.55
        # rotorR = rtip*cos(precone)
        rpm = addriver["RotSpd_mat"][1]
        omega = rpm*(2*pi)/60 #vinf*tsr/rotorR

        rho = addriver["FldDens"] #1.225
        mu = addriver["KinVisc"] #1.464e-5 #18.13e-6
        a = addriver["SpdSound"] #343.0
        shearexp = addriver["PLExp"][1] #0.0
        env = environment(rho, mu, a, vinf, omega, shearexp)

        pitch = 0.0
        Vx = vinf

        for idx = 1:n
            Vy = omega*rvec[idx]

            rotorout = Rotors.solve_BEM(rotor, blade, env, idx, Vx, Vy, pitch; npts=10)

            # af = AlphaAF("../data/polars/DU25_A17.dat", radians=false)
            polar = blade.airfoils[idx].polar
            af = CCBlade.AlphaAF(polar[:,1], polar[:,2], polar[:,3])
            section = Section(rvec[idx], chordvec[idx], twistvec[idx], af)
            rotor_cc = CCBlade.Rotor(blade.rhub, blade.rtip, B; turbine=true, tip=nothing)
            op = CCBlade.OperatingPoint(Vx, Vy, rho, pitch, mu, a)
            ccout = solve(rotor_cc, section, op)

            @test isa(rotorout, typeof(ccout))

            flags = Rotors.compare_fieldnames(rotorout, ccout)
            @test any(i -> i, flags)

            if ccout.a>1
                @show idx, ccout.a, ccout.phi
            end
        end


    end #End test CCBlade
end #End test BEMT

nothing

#=
Things I inserted into CCBlade that I might need: 

- Just after the zero Outputs initializer
Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G) = Outputs(promote(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)...)

=#