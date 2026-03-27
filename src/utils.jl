

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function dualcopy(x)
    if isa(x[1], ReverseDiff.TrackedReal) || isa(x[1], ForwardDiff.Dual)
        TF = typeof(x[1])
        ns = length(x)
        xnew = Vector{TF}(undef, ns)
        for i = 1:ns
            xnew[i] = x[i]
        end

        return xnew
    else
        # If not, then just return the original vector.
        return deepcopy(x)
    end
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


function linear_interp(xnew, x0, x1, y0, y1)
    top = y0*(x1-xnew) + y1*(xnew-x0)
    bot = x1-x0
    return top/bot
end


# """
#     brent_init(f, a, b, x0; args=(), atol=2e-12, rtol=4*eps(), maxiter=100)

# As Brent's method from FLOWMath, but with an initial guess. 

# **Inputs**
# """
# function brent_init(f, a, b, x0; args=(), atol=2e-12, rtol=4*eps(), maxiter=100)

#     fx0 = f(x0)
#     if isapprox(fx0, 0; atol, rtol)
#         return x0, (iter=0, fcalls=1, flag="CONVERGED")
#     end

#     fa = f(a)
#     fb = f(b)

#     #Todo: I need to come up with a way to select out of the three spots which bound to use. -> Instead of creating a new function... I could just check if x0 is inbetween a and b outside of the loop. If it is, then I can check if it is a zero. If not, then replace one of the bounds into Brent's method. 
#     if a<x0<b && fa*fx0<0
#         xpre = a; xcur = b
#     elseif x0<a
#         xpre = a; xcur = b
#     elseif b<x0
#         xpre = a; xcur = b
#     else #Original 
#         xpre = a; xcur = b
#     end
    
#     # xblk = 0.0; fblk = 0.0; spre = 0.0; scur = 0.0
#     error_num = "INPROGRESS"

#     fpre = f(xpre, args...)
#     fcur = f(xcur, args...)
#     xblk = zero(fpre); fblk = zero(fpre); spre = zero(fpre); scur = zero(fpre)
#     funcalls = 3
#     iterations = 0
    
#     if fpre*fcur > 0
#         error_num = "SIGNERR"
#         return 0.0, (iter=iterations, fcalls=funcalls, flag=error_num)
#     end
#     if fpre == zero(fpre)
#         error_num = "CONVERGED"
#         return xpre, (iter=iterations, fcalls=funcalls, flag=error_num)
#     end
#     if fcur == zero(fcur)
#         error_num = "CONVERGED"
#         return xcur, (iter=iterations, fcalls=funcalls, flag=error_num)
#     end

#     for i = 1:maxiter
#         iterations += 1
#         if fpre*fcur < 0
#             xblk = xpre
#             fblk = fpre
#             spre = scur = xcur - xpre
#         end
#         if abs(fblk) < abs(fcur)
#             xpre = xcur
#             xcur = xblk
#             xblk = xpre

#             fpre = fcur
#             fcur = fblk
#             fblk = fpre
#         end

#         delta = (atol + rtol*abs(xcur))/2
#         sbis = (xblk - xcur)/2
#         if fcur == zero(fcur) || abs(sbis) < delta
#             error_num = "CONVERGED"
#             return xcur, (iter=iterations, fcalls=funcalls, flag=error_num)
#         end

#         if abs(spre) > delta && abs(fcur) < abs(fpre)
#             if xpre == xblk
#                 # interpolate
#                 stry = -fcur*(xcur - xpre)/(fcur - fpre)
#             else
#                 # extrapolate
#                 dpre = (fpre - fcur)/(xpre - xcur)
#                 dblk = (fblk - fcur)/(xblk - xcur)
#                 stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre))
#             end
#             if 2*abs(stry) < min(abs(spre), 3*abs(sbis) - delta)
#                 # good short step
#                 spre = scur
#                 scur = stry
#             else
#                 # bisect
#                 spre = sbis
#                 scur = sbis
#             end
#         else 
#             # bisect
#             spre = sbis
#             scur = sbis
#         end

#         xpre = xcur; fpre = fcur
#         if abs(scur) > delta
#             xcur += scur
#         else
#             xcur += (sbis > 0 ? delta : -delta)
#         end

#         fcur = f(xcur, args...)
#         funcalls += 1
#     end
#     error_num = "CONVERR"
#     return xcur, (iter=iterations, fcalls=funcalls, flag=error_num)
# end



function rotate_x(alpha_x)
    return [
        1.0     0.0             0.0;
        0.0     cos(alpha_x)   -sin(alpha_x);
        0.0     sin(alpha_x)    cos(alpha_x)]
end

"""
    rotate_x(x, y, z, theta; T=false) -> xnew, ynew, znew

Rotate the vector x, y, z by theta about the X axis.
Use T to calculate the transpose (inverse). 

**Arguments**
- `x`, `y`, `z::Number`: The x, y, and z components of the vector
- `theta::Number`: The rotation angle (radians). 
- `T::Bool`: A flag of whether to calculate the rotation or inverse rotation. 
"""
function rotate_x(x, y, z, theta; T::Bool=false)

    st, ct = sincos(theta)

    if T
        xnew = 1*x + 0*y + 0*z
        ynew = 0*x + ct*y + st*z
        znew = 0*x - st*y + ct*z
        return xnew, ynew, znew
    else
        xnew = 1*x + 0*y + 0*z
        ynew = 0*x + ct*y - st*z
        znew = 0*x + st*y + ct*z
        return xnew, ynew, znew
    end
end

function rotate_y(alpha_y)
    return [cos(alpha_y) 0 sin(alpha_y);
            0.0 1.0 0.0;
            -sin(alpha_y) 0 cos(alpha_y)]
end

"""
    rotate_y(x, y, z, theta; T=false) -> xnew, ynew, znew

Rotate the vector x, y, z by theta about the Y axis.
Use T to calculate the transpose (inverse). 

**Arguments**
- `x`, `y`, `z::Number`: The x, y, and z components of the vector
- `theta::Number`: The rotation angle (radians). 
- `T::Bool`: A flag of whether to calculate the rotation or inverse rotation. 
"""
function rotate_y(x, y, z, theta; T::Bool=false)

    st, ct = sincos(theta)

    if T
        xnew = ct*x + 0*y - st*z
        ynew = 0*x + 1*y + 0*z
        znew = st*x + 0*y + ct*z
        return xnew, ynew, znew
    else
        xnew = ct*x + 0*y + st*z
        ynew = 0*x + 1*y + 0*z
        znew = -st*x + 0*y + ct*z
        return xnew, ynew, znew
    end
end

function rotate_z(alpha_z)
    return [cos(alpha_z) -sin(alpha_z) 0.0;
            sin(alpha_z) cos(alpha_z) 0.0;
            0.0 0.0 1.0]
end

"""
    rotate_z(x, y, z, theta; T=false) -> xnew, ynew, znew

Rotate the vector x, y, z by theta about the Z axis.
Use T to calculate the transpose (inverse). 

**Arguments**
- `x`, `y`, `z::Number`: The x, y, and z components of the vector
- `theta::Number`: The rotation angle (radians). 
- `T::Bool`: A flag of whether to calculate the rotation or inverse rotation. 
"""
function rotate_z(x, y, z, theta; T::Bool=false)

    st, ct = sincos(theta)

    if T #Calculate the transposed rotation. 
        xnew = ct*x + st*y + 0*z
        ynew = -st*x + ct*y + 0*z
        znew = 0*x + 0*y + 1*z
        return xnew, ynew, znew
    else
        xnew = ct*x - st*y + 0*z
        ynew = st*x + ct*y + 0*z
        znew = 0*x + 0*y + 1*z
        return xnew, ynew, znew
    end
end



function rotate_vector(x, y, z, theta_x, theta_y, theta_z; forward::Bool=true)
    sx, cx = sincos(theta_x)
    sy, cy = sincos(theta_y)
    sz, cz = sincos(theta_z)

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

function cross(a, b) #Todo: Test this
    i = a[2]*b[3] - a[3]*b[2]
    j = a[3]*b[1] - a[1]*b[3]
    k = a[1]*b[2] - a[2]*b[1]
    return i, j, k
end