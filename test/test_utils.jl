using Test, Rotors



@testset "Utility Functions" begin

    @testset "Rotations" begin
        ### Test rotation about same axis. 
        r = [0, 1, 0]
        x, y, z = Rotors.rotate_vector(r..., 0, 90*pi/180, 0)
        @test isapprox([x, y, z], r)

        ### Test simple rotation
        r = [1, 0, 0]
        x, y, z = Rotors.rotate_vector(r..., 0, 90*pi/180, 0)
        @test isapprox([x, y, z], [0, 0, -1])

        ### Test more complex rotation
        r = [0, 0, 1]
        x, y, z = Rotors.rotate_vector(r..., 45*pi/180, 0, 0)
        @test isapprox([x, y, z], [0, -sqrt(2)/2, sqrt(2)/2])

        ### Test multiple rotations at once.
        r = [0, 0, 1]
        x, y, z = Rotors.rotate_vector(r..., -pi/2, 0, pi/2)
        @test isapprox([x, y, z], [-1, 0, 0])
        
        ### Test multiple rotations strung together
        r = [1, 0, 0]
        x, y, z = Rotors.rotate_vector(Rotors.rotate_vector(r..., 0, pi/2, 0)..., pi/2, 0, 0)
        @test isapprox([x, y, z], [0, 1, 0])

        ### Test simple reverse rotations
        r = [0, 0, 1]
        x, y, z = Rotors.rotate_vector(r..., -pi/2, 0, -pi/2)
        @test isapprox([x, y, z], [1, 0, 0])
        x, y, z = Rotors.rotate_vector(x, y, z, -pi/2, 0, -pi/2; forward=false)
        @test isapprox([x, y, z], r)

        ### Test complex reverse rotations
        r = [0, 0, 1]
        x, y, z = Rotors.rotate_vector(r..., -pi/3, pi/47, pi/2)
        x, y, z = Rotors.rotate_vector(x, y, z, -pi/3, pi/47, pi/2; forward=false)
        @test isapprox([x, y, z], r)
    end

    @testset "Solvers" begin
        examplefun(x) = x^3
        a = -pi/2
        b = 1

        x0, info = Rotors.sub_brent(examplefun, a, b, 0.0)
        xgold = 0.0

        @test isapprox(x0, xgold, atol=1e-6)
    end #End testing solvers
end
