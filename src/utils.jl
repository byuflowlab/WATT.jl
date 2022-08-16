function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end

function isnanvec(vec)
    for i=1:length(vec)
        if isnan(vec[i])
            return true
        end
    end
    return false
end


function derivative_me(sol, tvec)
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

function plotdshistory(dshistory, tvec, index; legloc=:topright, titletext=nothing)
    nt = length(tvec)
    states = [dshistory[it][index] for it = 1:nt]
    x = zeros(nt, 4)
    for i = 1:nt
        x[i,:] = states[i].x
    end

    plt = plot(tvec, x[:,1], lab="State 1", xaxis="Time (s)", leg=legloc, title=titletext)
    plot!(tvec, x[:,2], lab="State 2")
    plot!(tvec, x[:,3], lab="State 3")
    plot!(tvec, x[:,4], lab="State 4")
    return plt
end