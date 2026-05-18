using Test

localpath = @__DIR__
cd(localpath)

@testset "WATT" begin
    include("test_utils.jl") #Checked 5/15/26
    include("test_types.jl") #Checked 5/18/26 - Todo: there are tests to be completed. 
    include("test_mesh.jl")  #checked 5/18/26 
    include("test_bem.jl") #Improved 5/18/26
    include("test_environments.jl") #Improved 5/18/26 
    include("test_gxbeam.jl") #Improved 5/18/26

    # Slow full-simulation regression test against a saved reference (~100 s sim).
    # Gates the Phase 1+ rework — any drift here means cleanup changed physics. #Implemented 5/15/26
    include("simple_NREL5MW.jl")

end