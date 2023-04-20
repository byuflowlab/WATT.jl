#=

=#

using Test
using Rotors, GXBeam, OpenFASTsr, DynamicStallModels
using LinearAlgebra, StaticArrays, Statistics


of = OpenFASTsr
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)

@testset "GXBeam" begin
    @testset "helper functions" begin

        ### Test converting Wiener-Milenkovic parameters to angles
        wmp = SVector(0.0, 0.0, 0.0)
        theta = Rotors.WMPtoangle(wmp)
        @test !any(i -> i!=0, theta)

        #Todo: Come up with more tests for this. 
    end
end #End testing GXBeam functions