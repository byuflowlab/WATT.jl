#=
Tests to test the four different BEM implementations: 1) single residual (inflow angle), 2) single fixed point iteration, 3) dual residual (induction factors), 4) dual fixed point iteration. 

Adam Cardoza
=#
using Revise
using Test, Rotors, DelimitedFiles, FLOWMath, CCBlade
# using Test, DelimitedFiles, FLOWMath, CCBlade

localpath = @__DIR__
cd(localpath)

include("/Users/adamcardoza/.julia/dev/Rotors/src/bem.jl")

errfun(x, xt) = 100*(x - xt)/xt

# @testset "Blade Element Momentum Theory" begin

    # Adam Cardoza 3/17/23
    # @testset "Dual fixed point iteration" begin

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
        a, ap, R = inductionfactors(phi0, U, omega, c, r, B, cx, cy, F) 

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
        rotor = Rotor(rhub, rtip, B; turbine=true, tip=nothing)
        op = CCBlade.OperatingPoint(U, omega*r, rho, pitch, mu, asound)
        out = solve(rotor, section, op)

        function residual(phir)
            alphar = phir - beta
            clr = clfit(alphar)
            cdr = cdfit(alphar)
            cxr = clr*cphi + cdr*sphi #x
            cyr = clr*sphi - cdr*cphi #x
            _, _, Ri = inductionfactors(phir, U, omega, c, r, B, cxr, cyr, F)
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



        

        #OpenFAST's Converged values
        # phi_gold = 0.00742743941*pi/180 #Tightening the convergence criteria by altering the source code you get 0.0001304168554254917 (radians)... 
        # a_gold = 0.997409395 
        # ap_gold = -0.100506008
        phi_gold = 0.0001304168554254917
        a = 0.997391823156366     
        ap = -9.984575201549731E-002

        println("")
        @show phistar #AD residual, flowmath brent
        @show phi2 #AD residual, sub_brent
        @show phi3 #CC residual, sub_brent
        @show phi4 #CC Residual, flowmath brent
        # @show phi, a, ap, converged
        @show out.phi, out.a, out.ap #CC residual, flowmath brent
        @show phi_gold, a_gold, ap_gold #OpenFAST

        phiAfbC = errfun(phistar, out.phi)
        phiAsbC = errfun(phi2, out.phi)
        phiCsbC = errfun(phi3, out.phi)
        phiOFC = errfun(phi_gold, out.phi)

        phiAfbO = errfun(phistar, phi_gold)
        phiAsbO = errfun(phi2, phi_gold)
        phiCsbO = errfun(phi3, phi_gold)
        phiCOF = errfun(out.phi, phi_gold)

        println("")
        println("")

        println("AD fb: ", phistar, ", ", phiAfbC, ", ", phiAfbO)
        println("AD sb: ", phi2, ", ", phiAsbC, ", ", phiAsbO)
        println("CC sb: ", phi3, ", ", phiCsbC, ", ", phiCsbO)
        println("OF   : ", phi_gold, ", ", phiOFC, ", ", 0.0)
        println("CCBla: ", out.phi, ", ", 0.0, ", ", phiCOF)
        

        # @show phierr, aerr, aperr #Phi converges within a half of a percent, but a and ap are like more than 10% off. Which is less fun. So maybe I should be comparing to OpenFAST, not CCBlade... because they are two different methods and I'm getting two different results. 

#     end #End Dual fixed point iteration
# end #End test BEMT

nothing

#=
Things I inserted into CCBlade that I might need: 

- Just before the Mach corrections: 
struct AFfun <: AFType
    clfun
    cdfun
end

function afeval(af::AFfun, alpha, Re, Mach) #I don't know if CCBlade will see this.
    return af.clfun(alpha), af.cdfun(alpha)
end


- Just after the zero Outputs initializer
Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G) = Outputs(promote(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)...)

=#