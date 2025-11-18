using Test

localpath = @__DIR__
cd(localpath)

@testset "WATT" begin
    include("test_utils.jl")
    include("test_types.jl")
    include("test_mesh.jl") 
    #Note: It is normal for warnings to show up, they are testing the boundaries of functions. 
    include("test_bem.jl")
    include("test_environments.jl")
    include("test_gxbeam.jl")

end