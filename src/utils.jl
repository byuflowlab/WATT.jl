

"""
    nearestto(xvec, x) -> (minval, minidx)

Find the element of `xvec` closest in value to `x` and return both the
value and its index.

**Arguments**
- `xvec::AbstractVector`: vector to search.
- `x::Number`: target value.

**Returns**
- `minval`: the entry of `xvec` nearest to `x`.
- `minidx::Int`: its index in `xvec`.
"""
function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

"""
    dualcopy(x) -> Vector

Copy `x` while preserving any AD tracking attached to its elements. When
`x` carries `ForwardDiff.Dual` or `ReverseDiff.TrackedReal` entries, the
result is a new `Vector` with the same dual element type so the tape /
seed survives. Falls back to `deepcopy` otherwise.

**Arguments**
- `x::AbstractVector`: vector whose element type may carry AD tracking.

**Returns**
- A copy of `x` with AD tracking preserved when present.
"""
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

"""
    getfieldnames(obj) -> Tuple{Symbol}

Convenience wrapper around `fieldnames(typeof(obj))`. Lets callers obtain
the field-name tuple of a value without first reaching for its type.
"""
function getfieldnames(obj)
    return fieldnames(typeof(obj))
end

"""
    compare_fieldnames(obj1, obj2) -> Vector{Bool}

Element-wise `isapprox` comparison of the fields of two same-typed
structs. Returns a `Vector{Bool}` aligned with `getfieldnames(obj1)` —
`true` for fields that match, `false` otherwise. If the two objects are
not of the same type, warns and returns all-`false`.

**Arguments**
- `obj1`, `obj2`: values to compare. Expected to share a concrete type.

**Returns**
- `Vector{Bool}` of per-field match flags.
"""
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


"""
    isnanvec(vec) -> Bool

Short-circuiting check for any `NaN` entry in `vec`. Returns `true` as
soon as one is found.

**Arguments**
- `vec::AbstractVector`: vector to scan.

**Returns**
- `true` if any element is `NaN`, `false` otherwise.
"""
function isnanvec(vec)
    for i=1:length(vec)
        if isnan(vec[i])
            return true
        end
    end
    return false
end


"""
    derivative_me(sol, tvec) -> Matrix

Evaluate an ODE solution `sol` at `tvec`, then differentiate each
component column-wise using an Akima spline through the sampled values.

**Arguments**
- `sol`: callable solution object (`sol(tvec)` must return a matrix-like result).
- `tvec::AbstractVector`: time points to evaluate / differentiate at.

**Returns**
- `Matrix` with the same `size(sol(tvec)')` whose columns are the time
  derivatives of the corresponding components of the solution.
"""
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

"""
    mat_derivative(data, tvec) -> Matrix

Column-wise time derivative of a sampled-signal matrix `data` (rows =
time, columns = signals) using an Akima spline through each column.

**Arguments**
- `data::AbstractMatrix`: `length(tvec) × n` matrix of samples.
- `tvec::AbstractVector`: time points corresponding to the rows of `data`.

**Returns**
- `Matrix` of the same shape as `data` holding the per-column derivatives.
"""
function mat_derivative(data, tvec)
    m,n = size(data)

    du = zero(data)

    for i = 1:n
        spline = Akima(tvec, data[:,i])
        du[:,i] = FLOWMath.gradient(spline, tvec)
    end
    return du
end


"""
    linear_interp(xnew, x0, x1, y0, y1) -> Number

Single-point linear interpolation between `(x0, y0)` and `(x1, y1)`
evaluated at `xnew`. No bounds check — extrapolates if `xnew` lies
outside `[x0, x1]`.

**Arguments**
- `xnew`: location to evaluate at.
- `x0`, `x1`: bracketing x-coordinates.
- `y0`, `y1`: corresponding y-values.

**Returns**
- Interpolated value at `xnew`.
"""
function linear_interp(xnew, x0, x1, y0, y1)
    top = y0*(x1-xnew) + y1*(xnew-x0)
    bot = x1-x0
    return top/bot
end


"""
    rotate_x(alpha_x) -> Matrix{Float64}

Build the 3×3 right-handed rotation matrix for an angle `alpha_x` about
the X axis. See [`rotate_x(x, y, z, theta; T=false)`](@ref) for the
allocation-free per-component version.
"""
function rotate_x(alpha_x)
    return [
        1.0     0.0             0.0;
        0.0     cos(alpha_x)   -sin(alpha_x);
        0.0     sin(alpha_x)    cos(alpha_x)]
end

"""
    rotate_x(x, y, z, theta; T=false) -> (xnew, ynew, znew)

Rotate the vector `(x, y, z)` by `theta` about the X axis without
allocating an intermediate matrix.

**Arguments**
- `x`, `y`, `z::Number`: components of the input vector.
- `theta::Number`: rotation angle (rad).

**Keyword Arguments**
- `T::Bool = false`: when `true`, apply the transposed (inverse) rotation.

**Returns**
- `(xnew, ynew, znew)`: components of the rotated vector.
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

"""
    rotate_y(alpha_y) -> Matrix{Float64}

Build the 3×3 right-handed rotation matrix for an angle `alpha_y` about
the Y axis. See [`rotate_y(x, y, z, theta; T=false)`](@ref) for the
allocation-free per-component version.
"""
function rotate_y(alpha_y)
    return [cos(alpha_y) 0 sin(alpha_y);
            0.0 1.0 0.0;
            -sin(alpha_y) 0 cos(alpha_y)]
end

"""
    rotate_y(x, y, z, theta; T=false) -> (xnew, ynew, znew)

Rotate the vector `(x, y, z)` by `theta` about the Y axis without
allocating an intermediate matrix.

**Arguments**
- `x`, `y`, `z::Number`: components of the input vector.
- `theta::Number`: rotation angle (rad).

**Keyword Arguments**
- `T::Bool = false`: when `true`, apply the transposed (inverse) rotation.

**Returns**
- `(xnew, ynew, znew)`: components of the rotated vector.
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

"""
    rotate_z(alpha_z) -> Matrix{Float64}

Build the 3×3 right-handed rotation matrix for an angle `alpha_z` about
the Z axis. See [`rotate_z(x, y, z, theta; T=false)`](@ref) for the
allocation-free per-component version.
"""
function rotate_z(alpha_z)
    return [cos(alpha_z) -sin(alpha_z) 0.0;
            sin(alpha_z) cos(alpha_z) 0.0;
            0.0 0.0 1.0]
end

"""
    rotate_z(x, y, z, theta; T=false) -> (xnew, ynew, znew)

Rotate the vector `(x, y, z)` by `theta` about the Z axis without
allocating an intermediate matrix.

**Arguments**
- `x`, `y`, `z::Number`: components of the input vector.
- `theta::Number`: rotation angle (rad).

**Keyword Arguments**
- `T::Bool = false`: when `true`, apply the transposed (inverse) rotation.

**Returns**
- `(xnew, ynew, znew)`: components of the rotated vector.
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



"""
    rotate_vector(x, y, z, theta_x, theta_y, theta_z; forward=true) -> NTuple{3}

Apply the composite Z–Y–X intrinsic rotation
`Rz(theta_z) · Ry(theta_y) · Rx(theta_x)` to the input vector
`(x, y, z)`. Pass `forward=false` to apply the transpose (inverse)
rotation instead.

**Arguments**
- `x`, `y`, `z`: vector components.
- `theta_x`, `theta_y`, `theta_z`: rotation angles about the X, Y, Z axes (rad).

**Keyword Arguments**
- `forward::Bool = true`: when `false`, applies the transposed rotation.

**Returns**
- `(x_new, y_new, z_new)`: rotated vector components.
"""
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

"""
    cross(a, b) -> NTuple{3}

Tuple-returning 3D cross product `a × b`. Avoids the array allocation of
`LinearAlgebra.cross` when only the scalar components are needed.

**Arguments**
- `a`, `b`: indexable 3-component vectors.

**Returns**
- `(i, j, k)`: the three components of `a × b`.
"""
function cross(a, b) #Todo: Test this
    i = a[2]*b[3] - a[3]*b[2]
    j = a[3]*b[1] - a[1]*b[3]
    k = a[1]*b[2] - a[2]*b[1]
    return i, j, k
end