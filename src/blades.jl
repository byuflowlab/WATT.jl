struct Airfoil{TF, Tfit}
    polar::Array{TF, 2}
    cl::Tfit
    cd::Tfit
    cm::Tfit
    A::Array{TF,1}
    b::Array{TF,1}
    T::Array{TF,1}
end

function simpleairfoil(polar)
    cl = Akima(polar[:,1], polar[:,2])
    cd = Akima(polar[:,1], polar[:,3])
    cm = Akima(polar[:,1], zeros(length(polar[:,1])))
    A = zeros(2)
    b = zeros(2)
    T = zeros(2)
    return Airfoil(polar, cl, cd, cm, A, b, T)
end

struct Blade
    airfoils
end