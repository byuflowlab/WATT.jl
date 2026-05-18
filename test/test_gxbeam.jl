#=

=#

using Test
using WATT, GXBeam, OpenFASTTools, DynamicStallModels
using LinearAlgebra, StaticArrays, Statistics


of = OpenFASTTools
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)

@testset "GXBeam" begin
    @testset "helper functions" begin

        ### Test converting Wiener-Milenkovic parameters to angles
        wmp = SVector(0.0, 0.0, 0.0)
        theta = WATT.WMPtoangle(wmp)
        @test !any(i -> i!=0, theta)

        #Todo: Come up with more tests for this. 
    end
end #End testing GXBeam functions