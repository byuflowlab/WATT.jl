using Test, Rotors



@testset "Utility Functions" begin

    @testset "Solvers" begin
        examplefun(x) = x^3
        a = -pi/2
        b = 1

        x0, info = Rotors.sub_brent(examplefun, a, b, 0.0)
        xgold = 0.0

        @test isapprox(x0, xgold, atol=1e-6)
    end #End testing solvers
end
