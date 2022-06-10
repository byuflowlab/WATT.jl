

function derivative(sol, tvec)
    solatt = sol(tvec)
    x = Array(solatt)'
    m, n = size(x)
    du = zero(x)

    for i= 1:n
        spline = Akima(tvec, x[:,i])
        du[:,i] = FLOWMath.gradient(spline, tvec)
    end
    return du
end

function mat_derivative(data, tvec)
    m,n = size(data)

    du = zero(data)

    for i = 1:n
        spline = Akima(tvec, data[:,i])
        du[:,i] = FLOWMath.gradient(spline, tvec)
    end
    return du
end

function prepareextra(tspan, filename, numitems)
    mat = [tspan[1] repeat(["Initial"], numitems)...]
    writedlm(filename, mat, ',')
end

function saveextra(t, filename, items...; verbose=false, saveiteration=false) #Todo: There is a saving Callback apparently. 
    file = readdlm(filename,',')

    # newitems = Array{eltype(items)}(undef, length(items)) #Todo: Needs to be adapted to take any base type. 
    newitems = zeros(length(items))

    for i = 1:length(items)
        if isa(items[i], ForwardDiff.Dual)
            newitems[i] = items[i].value
        else
            newitems[i] = items[i]
        end
    end

    tvec = file[:,1]
    line = [t newitems...]

    if tvec[end] != t #Save a new line if the solver has moved forward a timestep
        newfile = vcat(file, line)
        writedlm(filename, newfile, ',')

    else #Save over the last line if the solver hasn't moved forward a timestep
        newfile = vcat(file[1:end-1,:], line)
        writedlm(filename, newfile, ',')

    end
end

function readextra(filename)
    file = readdlm(filename, ',')
    tvec = file[:,1]

    m, n = size(file)
    keepidxs = [1]
    for i = 2:m
        if !(tvec[i]<=tvec[1])
            push!(keepidxs, i)
        end
    end
    return file[keepidxs, :]
end