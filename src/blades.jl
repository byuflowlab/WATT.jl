struct Blade #TODO: Do I want structural information in this? Or a different struct? or any struct at all. 
    airfoils
end

struct Airfoil{TF, Tfit}
    polar::Array{TF, 2}
    cl::Tfit
    cd::Tfit
    cm::Tfit
    dcldalpha::TF
    alpha0::TF
    A::Array{TF,1}
    b::Array{TF,1}
    T::Array{TF,1}
end

function simpleairfoil(polar)
    cl = Akima(polar[:,1], polar[:,2])
    cd = Akima(polar[:,1], polar[:,3])
    cm = Akima(polar[:,1], zeros(length(polar[:,1])))
    dcldalpha = 2*pi
    alpha0 = 0.0
    A = [0.3, 0.7]
    b = [0.14, 0.53]
    T = [1.7, 3.0]
    return Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, A, b, T)
end

function nearestto(xvec, x) #TODO: Move to a utilities file. 
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function complexairfoil(polar; A = [0.3, 0.7], b = [0.14, 0.53], T = [1.7, 3.0])
    cl = Akima(polar[:,1], polar[:,2])
    cd = Akima(polar[:,1], polar[:,3])
    cm = Akima(polar[:,1], zeros(length(polar[:,1])))
    alpha0, _ = brent(cl, -0.25, 0.25)
    _, minclidx = findmin(polar[:,2])
    _, maxclidx = findmax(polar[:,2])
    middlepolar = polar[minclidx:maxclidx,:]
    _, cl0idx = nearestto(middlepolar[:,2], 0.0)
    alpha50 = middlepolar[end,1]*0.25
    _, alf50idx = nearestto(middlepolar[:,1], alpha50)
    _, dcldalpha = linear_fit(middlepolar[cl0idx:alf50idx,1], middlepolar[cl0idx:alf50idx,2]) #TODO: Create my own linear fit function so I don't have to pull in a package. 
    return Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, A, b, T)
end

