


abstract type ModelInit end #TODO: I'd like to generalize this to both the ds model and gxbeam (and anything else for that matter. )

struct Hansen <: ModelInit
end

struct Steady <: ModelInit
end