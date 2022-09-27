


abstract type ModelInit end #TODO: I'd like to generalize this to both the ds model and gxbeam (and anything else for that matter. )

struct Hansen <: ModelInit
end

struct BeddoesLeishman <: ModelInit #Todo: This is a repeated struct from DSM. I should change the name. 
end

struct Steady <: ModelInit
end