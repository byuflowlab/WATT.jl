using Test

localpath = @__DIR__
cd(localpath)

@testset "WATT" begin
    include("test_utils.jl") #Checked 5/15/26
    include("test_types.jl")
    include("test_mesh.jl") 
    #Note: It is normal for warnings to show up, they are testing the boundaries of functions. 
    include("test_bem.jl")
    include("test_environments.jl")
    include("test_gxbeam.jl")

    # Slow full-simulation regression test against a saved reference (~100 s sim).
    # Gates the Phase 1+ rework — any drift here means cleanup changed physics.
    include("simple_NREL5MW.jl")

end