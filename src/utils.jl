

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end

function compare_fieldnames(obj1, obj2)
    
    names = getfieldnames(obj1)
    flags = Vector{Bool}(undef, length(names))

    if !isa(obj1, typeof(obj2))
        @warn("compare_fieldnames(): the objects are not of the same type.")
        flags .= false
        return flags
    end

    for i in eachindex(names)
        i1 = getproperty(obj1, names[i])
        i2 = getproperty(obj2, names[i])
        if isapprox(i1, i2)
            flags[i] = true
        else
            flags[i] = false
        end
    end
    return flags
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

"""
    sub_brent()

Translated and adapted from OpenFAST
"""
function sub_brent(fun, a, b, toler; maxiter::Int = 100, xtoler=1e-6, epsilon=eps())

    # Set the user chosen tolerance t to xtoler
    if (xtoler<0.0)
        @warn("WARNING: xtoler must be positive. Resetting xtoler.")
        xtoler = 0.0
    end

    fa = fun(a)
    fb = fun(b)

    if fa==0
        # println("Initial boundary is a zero.")
        return a, fa
    elseif fb == 0
        # println("Initial boundary is a zero.")
        return b, (fb, 0)
    end

    # Test whether root is bracketed
    if !(fa*fb<0)
        if abs(fa)<abs(fb)
            @warn("brent: WARNING: root is not bracketed, returning best endpoint a = $a")
            return a, (fa, 0)
        else
            @warn("brent: WARNING: root is not bracketed, returning best endpoint b = $b")
            return b, (fb, 0)
        end
    end

    # step = 'init'

    c = a
    fc = fa
    e = b - a
    d = e

    # At any point in time, b is the best guess of the root, a is the previous value of b, and the root is bracketed by b and c.
    for iter = 1:maxiter

        if (fb>0.0 && fc>0.0) || (fb<=0.0 && fc<=0.0)
            c = a
            fc = fa
            e = b - a
            d = e
        end

        # If c is strictly better than b, swap b and c so b is the best guess. 
        if abs(fc)<abs(fb)
            a = b
            b = c
            c  = a
            fa = fb
            fb = fc
            fc = fa
        end

        # Set the tolerance. Note: brent is very careful with these things, so don't deviate from this.
        # tol = 2.0*epsilon*abs(b) + xtoler
        tol = 2e-12

        # Determine what half the length of the bracket [b,c] is
        m = 0.5*(c-b)

        # If taking a bisection step would move the guess of the root less than tol, then return b the best guess.
        # if (abs(m)<=tol)
        #     # println("bracket size below tolerance")
        #     return b, (fb, iter)
        # elseif (fb==0.0)
        #     # println("Residual converged mid step.")
        #     return b, (fb, iter)
        # end

        if (fb==0.0)
            println("Residual converged mid step.")
            return b, (fb, iter)
        end


        # If still here, then check whether need to do bisection or can do interpolation
        if ((abs(e)>=tol) && (abs(fa)>abs(fb)))
            s = fb/fa
            if a != c
                # Inverse quadratic interpolation
                q = fa/fc
                r = fb/fc
                p = s*(2.0*m*q*(q-r) - (b-a)*(r-1.0))
                q = (q-1.0)*(r-1.0)*(s-1.0)
                
                # step = 'quad'
            else
                # Linear interpolation
                p = 2.0*m*s
                q = 1.0-s

                # step = 'linear'
            end

            # Ensure p is positive
            if p<=0.0
                p = -p
            else
                q = -q
            end

            s = e
            e = d

            if (2.0*p>=3.0*m*q-abs(tol*q)) || (p>=abs(0.5*s*q))
                # Interpolation step failed to produce good step, bisect instead
                e = m
                d = m # m is half the distance between b and c
                # step = 'bisect'
            else
                # Do interpolation step (either quadratic or linear)
                d = p/q
            end
        else

            # Do bisection step
            e = m 
            d = m
            
        end 

        # Get new points. 
        ## Replace a (the old b) with b.
        a = b
        fa = fb

        ### Increment b by d if that is greater than the tolerance. O/w, increment by tol.
        if abs(d)<=tol
            # m is .5*(c-b) with the bracket either [b,c] or [c,b]. 
            if m > 0.0
                # If m>0d0, then bracket is [b,c] so move towards c by tol
                b = b + tol
            else
                # If m<=0d0, then bracket is [c,b] so move towards c by tol
                b = b - tol
            end
        else
            b = b + d
        end

        ### Evaluate at the new point
        fb = fun(b)

        # Check my custom tolerance 
        if abs(fb)<toler
            println("Residual Converged below tolerance.")
            return b, (fb, iter)
        end
            
    end #End convergence iterations
    return b, (fb, maxiter)
end

function rotate_vector(x, y, z, theta_x, theta_y, theta_z; forward::Bool=true)
    cx = cos(theta_x)
    cy = cos(theta_y)
    cz = cos(theta_z)
    sx = sin(theta_x)
    sy = sin(theta_y)
    sz = sin(theta_z)

    if forward
        x_new = x*(cz*cy) + y*(cz*sy*sx - sz*cx) + z*(cz*sy*cx + sz*sx)
        y_new = x*(cy*sz) + y*(sx*sy*sz + cx*cz) + z*(cx*sy*sz - sx*cz)
        z_new = x*(-sy) + y*(sx*cy) + z*(cx*cy)
        return x_new, y_new, z_new
    else
        # error("rotate_vector: Reverse transformation not implemented yet.")
        x_new = x*(cz*cy) + y*(cy*sz) + z*(-sy)
        y_new = x*(cz*sy*sx - sz*cx) + y*(sx*sy*sz + cx*cz) + z*(sx*cy)
        z_new = x*(cz*sy*cx + sz*sx) + y*(cx*sy*sz - sx*cz) + z*(cx*cy)
        return x_new, y_new, z_new
    end

end