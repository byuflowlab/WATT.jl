import CCBlade.afeval

export Blade

struct Blade{TF} #TODO: Do I want structural information in this? Or a different struct? or any struct at all. 
    rhub::TF
    rtip::TF
    rR::Array{TF,1}
    airfoils::AbstractVector{<:DS.Airfoil}
end



# function (af::ds.Airfoil)(alpha, Re, Mach)
#     return af.cl(alpha), af.cd(alpha)
# end

# function afeval(af::ds.Airfoil, alpha, Re, Mach)
#     return af(alpha, Re, Mach)
# end

function afeval(af::DS.Airfoil, alpha, Re, Mach)
    return af.cl(alpha), af.cd(alpha)
end



