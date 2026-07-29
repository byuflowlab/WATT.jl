using Test

localpath = @__DIR__
cd(localpath)

@testset "WATT" begin
    include("test_utils.jl") #Checked 5/15/26
    include("test_types.jl") #Checked 5/18/26 - Todo: there are tests to be completed. 
    include("test_mesh.jl")  #checked 5/18/26 
    include("test_bemt.jl") #Improved 5/18/26
    include("test_environments.jl") #Improved 5/18/26 
    include("test_gxbeam.jl") #Improved 5/18/26
    include("test_dynamicstall.jl") #Added 2026-05-19 (Phase 2)
    include("test_aero_only.jl") #Added 2026-05-19 (Phase 2) — mostly @test_broken pending Phase 3/5
    include("test_aerostructural.jl") #Added 2026-05-19 (Phase 2)
    include("test_step_solution.jl") #Added 2026-07-29 — single-step primitive equivalence gate
    include("test_static.jl") #Added 2026-05-19 (Phase 2)

    # AD tests are slow — gate behind WATT_AD_TESTS so CI can opt in.
    if get(ENV, "WATT_AD_TESTS", "false") == "true"
        include("test_ad.jl") #Added 2026-05-19 (Phase 2)
    end

    # Slow full-simulation regression test against a saved reference (~100 s sim).
    # Gates the Phase 1+ rework — any drift here means cleanup changed physics. #Implemented 5/15/26
    include("simple_NREL5MW.jl")

end