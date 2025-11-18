using Test, WATT



@testset "Utility Functions" begin

    @testset "Rotations" begin
        ### Test rotation about same axis. 
        r = [0, 1, 0]
        theta = 90*pi/180
        x, y, z = WATT.rotate_vector(r..., 0, theta, 0)
        x2, y2, z2 = WATT.rotate_y(r..., theta)
        @test isapprox([x, y, z], r)
        @test isapprox([x2, y2, z2], r)

        ### Test simple rotation
        r = [1, 0, 0]
        x, y, z = WATT.rotate_vector(r..., 0, theta, 0)
        x2, y2, z2 = WATT.rotate_y(r..., theta)
        r_gold = [0, 0, -1]
        @test isapprox([x, y, z], r_gold)
        @test isapprox([x2, y2, z2], r_gold)

        ### Test more complex rotation
        r = [0, 0, 1]
        theta = 45*pi/180
        x, y, z = WATT.rotate_vector(r..., theta, 0, 0)
        x2, y2, z2 = WATT.rotate_x(r..., theta)
        r_gold = [0, -sqrt(2)/2, sqrt(2)/2]
        @test isapprox([x, y, z], r_gold)
        @test isapprox([x2, y2, z2], r_gold)

        ### Test multiple rotations at once.
        r = [0, 0, 1]
        theta = pi/2
        x, y, z = WATT.rotate_vector(r..., -theta, 0, theta)
        x2, y2, z2 = WATT.rotate_z(WATT.rotate_x(r..., -theta)..., theta)
        r_gold = [-1, 0, 0]
        @test isapprox([x, y, z], r_gold)
        @test isapprox([x2, y2, z2], r_gold)
        
        ### Test multiple rotations strung together
        r = [1, 0, 0]
        x, y, z = WATT.rotate_vector(WATT.rotate_vector(r..., 0, pi/2, 0)..., pi/2, 0, 0)
        x2, y2, z2 = WATT.rotate_x(WATT.rotate_y(r..., theta)..., theta)
        r_gold = [0, 1, 0]
        @test isapprox([x, y, z], r_gold)
        @test isapprox([x2, y2, z2], r_gold)

        ### Test simple reverse rotations
        r = [0, 0, 1]
        x, y, z = WATT.rotate_vector(r..., -pi/2, 0, -pi/2)
        @test isapprox([x, y, z], [1, 0, 0])
        x, y, z = WATT.rotate_vector(x, y, z, -pi/2, 0, -pi/2; forward=false)
        @test isapprox([x, y, z], r)

        ### Test complex reverse rotations
        r = [0, 0, 1]
        x, y, z = WATT.rotate_vector(r..., -pi/3, pi/47, pi/2)
        x, y, z = WATT.rotate_vector(x, y, z, -pi/3, pi/47, pi/2; forward=false)
        @test isapprox([x, y, z], r)

        ### Test by comparison to individual, slower rotation matrices
        V = [1, 2, 3]
        theta = [10, 38, 73.4].*(pi/180)

        R = WATT.rotate_z(theta[3])*WATT.rotate_y(theta[2])*WATT.rotate_x(theta[1])

        #Rotate the 
        Vt = R*V
        Vt2 = WATT.rotate_vector(V..., theta...; forward=true)

        Vb = R'*Vt
        Vb2 = WATT.rotate_vector(Vt2..., theta...; forward=false)

        @test isapprox(R', inv(R))
        forwardflag = true
        reverseflag = true

        for i = 1:3
            if !isapprox(Vt[i], Vt2[i]) 
                forwardflag = false
            end

            if !isapprox(Vb[i], Vb2[i])
                reverseflag = false
            end
        end
        @test forwardflag
        @test reverseflag
    
    end

    @testset "Solvers" begin
        examplefun(x) = x^3
        a = -pi/2
        b = 1

        x0, info = WATT.sub_brent(examplefun, a, b, 0.0)
        xgold = 0.0

        @test isapprox(x0, xgold, atol=1e-6)
    end #End testing solvers
end
