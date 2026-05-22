#=
Code used to interact with GXBeam.jl. 
Also some code to run GXBeam in the way that I'm doing it in the aerostructural simulation. 

=#
"""
    retrieve_eulerangles(R) -> SVector{3}

Extract the intrinsic 3-2-1 (roll, pitch, yaw) Euler angles from a 3×3
rotation matrix `R`.

**Arguments**
- `R`: 3×3 rotation matrix.

**Returns**
- `SVector{3}` of `(roll, pitch, yaw)` in radians.
"""
function retrieve_eulerangles(R)
    return SVector{3}([atan(R[3,2], R[3,3]), asin(-R[3,1]), atan(R[2,1], R[1,1])])
end

"""
    WMPtoangle(c) -> SVector{3}

Convert a Wiener–Milenkovic rotation parameter vector `c` into the
equivalent 3-2-1 Euler angles. Applies GXBeam's parameter scaling before
building the rotation matrix.

**Arguments**
- `c`: Wiener–Milenkovic 3-vector.

**Returns**
- `SVector{3}` of `(roll, pitch, yaw)` in radians.
"""
function WMPtoangle(c)
    scaling = GXBeam.rotation_parameter_scaling(c)
    R = GXBeam.wiener_milenkovic(scaling*c)'
    return retrieve_eulerangles(R)
end

"""
    get_bladelength_vector(assembly) -> Vector

Cumulative arc-length along the assembly's point list, measured from the
origin to each point. The first entry is `‖p₁‖`, and each subsequent
entry adds the segment length `‖pᵢ − pᵢ₋₁‖`. Useful for mapping
GXBeam points back to a 1D radial coordinate.

**Arguments**
- `assembly::GXBeam.Assembly`: assembly to measure.

**Returns**
- `Vector` of arc-lengths, one per point.
"""
function get_bladelength_vector(assembly::GXBeam.Assembly)
    inittype = eltype(assembly)
    ns = length(assembly.points)
    rgx = zeros(inittype, ns)
    rgx[1] = norm(assembly.points[1])
    for i = 2:ns
        delr = assembly.points[i] - assembly.points[i-1]
        rgx[i] = rgx[i-1] + norm(delr)
    end
    return rgx
end


"""
    update_assembly(assembly; compliance=nothing, stiffness=nothing) -> GXBeam.Assembly

Return a copy of `assembly` with the per-element compliance matrices
replaced. Either `compliance` or `stiffness` may be supplied — if
`stiffness` is given it is inverted element-by-element. A single matrix
broadcasts to all elements; a vector of matrices must be the same length
as `assembly.elements`.

**Arguments**
- `assembly::GXBeam.Assembly`: assembly to copy.
- `compliance`: a 6×6 matrix or a `length(elements)`-long vector of 6×6 matrices.
- `stiffness`: a 6×6 matrix or a `length(elements)`-long vector of 6×6 matrices; inverted to produce compliance.

**Returns**
- A new `GXBeam.Assembly` with the same points/connectivity but updated element compliances.
"""
function update_assembly(assembly; compliance=nothing, stiffness=nothing)

    points = assembly.points

    start = assembly.start
    stop = assembly.stop

    ### Update Elements
    if !isnothing(compliance) || !isnothing(stiffness)
        if !isnothing(stiffness)
            if isa(stiffness, Vector) #Todo: This doesn't allow for unique replacement of the compliance matrix. (As stands this just gives a stiffness to an entire matrix... yeah? Unless... it's a vector of matrices... which might work as coded.)
                if length(stiffness) != length(assembly.elements)
                    error("The number of elements in the stiffness vector does not match the number of elements in the assembly.")
                end
                compliance = [SMatrix{6,6}(inv(stiffness[i])) for i = eachindex(assembly.elements)]
            else
                compliance = [inv(stiffness) for _ in eachindex(assembly.elements)]
            end
        else
            if isa(compliance, Vector)
                if length(compliance) != length(assembly.elements)
                    error("The number of elements in the compliance vector does not match the number of elements in the assembly.")
                end
            else
                compliance = [SMatrix{6,6}(compliance) for i = eachindex(assembly.elements)]
            end
        end
    end

    elements = [GXBeam.Element(assembly.elements[i].L, assembly.elements[i].x, compliance[i], assembly.elements[i].mass, assembly.elements[i].Cab, assembly.elements[i].mu) for i = eachindex(assembly.elements)]

    

    return GXBeam.Assembly(points, start, stop, elements)

end

"""
    interpolate_matrix_symmetric(f1, f2, Kmat; fit=Linear)

Interpolate a set of symmetric 6×6 matrices from one vector of fractional
locations to any number of new fractions. Only the upper triangle is
interpolated; the lower triangle is filled by symmetry at each output
location.

**Arguments**
- `f1::Vector`: fractional locations of the input matrices.
- `f2::Vector`: fractional locations to interpolate to.
- `Kmat::Array{Float64, 3}`: input matrices, sized `6×6×length(f1)`.
- `fit`: interpolant constructor (default `Linear`).

**Returns**
- `Kfit::Array{Float64, 3}`: interpolated matrices, sized `6×6×length(f2)`.
"""
function interpolate_matrix_symmetric(f1, f2, Kmat; fit=Linear)
    if length(f1)!=size(Kmat, 3)
        error("interpolate_stiffness(): The number of f1 points and stiffness matrices must be the same.")
    end

    K11fit = fit(f1, Kmat[1,1,:])
    K12fit = fit(f1, Kmat[1,2,:])
    K13fit = fit(f1, Kmat[1,3,:])
    K14fit = fit(f1, Kmat[1,4,:])
    K15fit = fit(f1, Kmat[1,5,:])
    K16fit = fit(f1, Kmat[1,6,:])

    K22fit = fit(f1, Kmat[2,2,:])
    K23fit = fit(f1, Kmat[2,3,:])
    K24fit = fit(f1, Kmat[2,4,:])
    K25fit = fit(f1, Kmat[2,5,:])
    K26fit = fit(f1, Kmat[2,6,:])

    K33fit = fit(f1, Kmat[3,3,:])
    K34fit = fit(f1, Kmat[3,4,:])
    K35fit = fit(f1, Kmat[3,5,:])
    K36fit = fit(f1, Kmat[3,6,:])

    K44fit = fit(f1, Kmat[4,4,:])
    K45fit = fit(f1, Kmat[4,5,:])
    K46fit = fit(f1, Kmat[4,6,:])

    K55fit = fit(f1, Kmat[5,5,:])
    K56fit = fit(f1, Kmat[5,6,:])

    K66fit = fit(f1, Kmat[6,6,:])

    Kfit = zeros(6,6,length(f2))
    for i in eachindex(f2)
        Kfit[1,1,i] = K11fit(f2[i])
        Kfit[1,2,i] = K12fit(f2[i])
        Kfit[1,3,i] = K13fit(f2[i])
        Kfit[1,4,i] = K14fit(f2[i])
        Kfit[1,5,i] = K15fit(f2[i])
        Kfit[1,6,i] = K16fit(f2[i])

        Kfit[2,1,i] = K12fit(f2[i])
        Kfit[2,2,i] = K22fit(f2[i])
        Kfit[2,3,i] = K23fit(f2[i])
        Kfit[2,4,i] = K24fit(f2[i])
        Kfit[2,5,i] = K25fit(f2[i])
        Kfit[2,6,i] = K26fit(f2[i])

        Kfit[3,1,i] = K13fit(f2[i])
        Kfit[3,2,i] = K23fit(f2[i])
        Kfit[3,3,i] = K33fit(f2[i])
        Kfit[3,4,i] = K34fit(f2[i])
        Kfit[3,5,i] = K35fit(f2[i])
        Kfit[3,6,i] = K36fit(f2[i])

        Kfit[4,1,i] = K14fit(f2[i])
        Kfit[4,2,i] = K24fit(f2[i])
        Kfit[4,3,i] = K34fit(f2[i])
        Kfit[4,4,i] = K44fit(f2[i])
        Kfit[4,5,i] = K45fit(f2[i])
        Kfit[4,6,i] = K46fit(f2[i])

        Kfit[5,1,i] = K15fit(f2[i])
        Kfit[5,2,i] = K25fit(f2[i])
        Kfit[5,3,i] = K35fit(f2[i])
        Kfit[5,4,i] = K45fit(f2[i])
        Kfit[5,5,i] = K55fit(f2[i])
        Kfit[5,6,i] = K56fit(f2[i])

        Kfit[6,1,i] = K16fit(f2[i])
        Kfit[6,2,i] = K26fit(f2[i])
        Kfit[6,3,i] = K36fit(f2[i])
        Kfit[6,4,i] = K46fit(f2[i])
        Kfit[6,5,i] = K56fit(f2[i])
        Kfit[6,6,i] = K66fit(f2[i])
    end

    return Kfit
end


"""
    pane_assembly(assembly; ne=nothing, verbose=false, fit=Linear)

Resample the GXBeam point distribution of `assembly` to `ne` equal-length
elements along the arc-length of the original point list. Coordinates of
the new points are interpolated from the original ones using the provided
`fit` constructor.

!!! warning
    Currently returns the per-element non-dimensional arc-length vector
    `sevec_new`, not a rebuilt `GXBeam.Assembly`. Treat the return value
    as a meshing helper, not a full assembly replacement.

**Arguments**
- `assembly::GXBeam.Assembly`: assembly to resample.
- `ne::Int`: target number of elements. Pass `nothing` to no-op and return `assembly`.
- `verbose::Bool`: print a warning when `ne` is clamped to 1.
- `fit`: interpolant constructor used to resample point coordinates.

**Returns**
- `sevec_new::Vector`: non-dimensional midpoint locations of the new elements (see warning above).
"""
function pane_assembly(assembly; ne=nothing, verbose::Bool=false, fit=Linear)
    if isnothing(ne)
        return assembly
    elseif ne<1
        if verbose
            @warn("The number of elements must be at least one. Increasing number of elements to 1.")
        end
        ne = 1
    end


    ### Extract the current assembly points
    np_cur = length(assembly.points)
    ne_cur = np_cur - 1

    points = zeros(np_cur, 3)
    for i = 1:np_cur
        points[i,:] = collect(assembly.points[i].data)
    end

    ### Calculate the element lengths (so I can calculate the assembly length)
    lvec = zeros(ne_cur)
    for i = 2:np_cur
        p1 = view(points, i, :)
        p2 = view(points, i-1, :)
        delp = p2.-p1
        
        lvec[i-1] = sqrt(sum(delp.^2))
    end

    L = sum(lvec) #Calculate the assembly length

    ### Calculate non-dimensional position of points along beam
    svec = zeros(np_cur)
    for i = 2:np_cur
        svec[i] = sum(lvec[1:i-1])/L
    end

    Xfit = fit(svec, points[:,1])
    Yfit = fit(svec, points[:,2])
    Zfit = fit(svec, points[:,3])

    # elvec = [assembly.elements[i].L for i = 1:ne_cur]
    sevec = [(sum(lvec[1:i-1])+(lvec[i]/2))/L for i in 1:ne_cur]


    ### New Beam
    np = ne + 1
    svec_new = collect(range(0.0, 1.0, np))
    points_new = zeros(np, 3)

    points_new[:,1] = Xfit.(svec_new)
    points_new[:,2] = Yfit.(svec_new)
    points_new[:,3] = Zfit.(svec_new)

    pointsvec = [SVector{3}(points_new[i,:]) for i = 1:np]
    newlvec = (L/ne).*ones(ne)

    x_elements = [SVector{3}((points_new[i,:].+points_new[i+1,:])./2) for i in 1:ne]

    elvec_new = [(L/(2*ne)) + (i-1)*(L/ne) for i = 1:ne]

    sevec_new = [elvec_new[i]/L for i in 1:ne]

    return sevec_new #Todo: Shouldn't this return a new assembly? Where is this used? 
end

"""
    BladePoints(points; xdim=true, ydim=true, zdim=true)

Plot-recipe wrapper around a vector of 3-component points. Use as
`plot(BladePoints(points))` after `using Plots`. Drops any axis whose
coordinates are all zero, producing a 2D plot when the points are planar.
"""
struct BladePoints{T}
    points::T
    xdim::Bool
    ydim::Bool
    zdim::Bool
end
BladePoints(points; xdim=true, ydim=true, zdim=true) = BladePoints(points, xdim, ydim, zdim)

@recipe function f(bp::BladePoints)
    n = length(bp.points)
    x = [bp.points[i][1] for i in 1:n]
    y = [bp.points[i][2] for i in 1:n]
    z = [bp.points[i][3] for i in 1:n]

    xdim, ydim, zdim = bp.xdim, bp.ydim, bp.zdim
    if iszero(x); xdim = false; end
    if iszero(y); ydim = false; end
    if iszero(z); zdim = false; end

    seriestype := :scatter
    legend --> false
    aspect_ratio --> :equal

    if xdim && ydim && zdim
        xguide --> "X distance"
        yguide --> "Y distance"
        zguide --> "Z distance"
        return x, y, z
    elseif !xdim
        xguide --> "Y distance"
        yguide --> "Z distance"
        return y, z
    elseif !ydim
        xguide --> "X distance"
        yguide --> "Z distance"
        return x, z
    else
        xguide --> "X distance"
        yguide --> "Y distance"
        return x, y
    end
end

"""
    AssemblyPlot(assembly; xdim=true, ydim=true, zdim=true)

Plot-recipe wrapper around a `GXBeam.Assembly`. Use as
`plot(AssemblyPlot(assembly))` after `using Plots`. Draws the assembly
points connected by a black line with the element midpoints overlaid;
collapses to 2D if one coordinate is all zero.
"""
struct AssemblyPlot{A}
    assembly::A
    xdim::Bool
    ydim::Bool
    zdim::Bool
end
AssemblyPlot(assembly; xdim=true, ydim=true, zdim=true) = AssemblyPlot(assembly, xdim, ydim, zdim)

@recipe function f(ap::AssemblyPlot)
    points = cat(ap.assembly.points...; dims=2)'
    nel = length(ap.assembly.elements)
    elements = zeros(nel, 3)
    for i in 1:nel
        elements[i, :] = ap.assembly.elements[i].x
    end

    xdim, ydim, zdim = ap.xdim, ap.ydim, ap.zdim
    if iszero(points[:, 1]); xdim = false; end
    if iszero(points[:, 2]); ydim = false; end
    if iszero(points[:, 3]); zdim = false; end

    legend --> false
    aspect_ratio --> :equal

    if xdim && ydim && zdim
        xguide --> "X distance"
        yguide --> "Y distance"
        zguide --> "Z distance"

        @series begin
            seriestype := :path
            linewidth --> 4
            linecolor --> :black
            points[:, 1], points[:, 2], points[:, 3]
        end
        @series begin
            seriestype := :scatter
            points[:, 1], points[:, 2], points[:, 3]
        end
        @series begin
            seriestype := :scatter
            markershape --> :cross
            elements[:, 1], elements[:, 2], elements[:, 3]
        end
    else
        if !xdim
            xplt, yplt = points[:, 2], points[:, 3]
            xmplt, ymplt = elements[:, 2], elements[:, 3]
            xguide --> "Y distance"
            yguide --> "Z distance"
        elseif !ydim
            xplt, yplt = points[:, 1], points[:, 3]
            xmplt, ymplt = elements[:, 1], elements[:, 3]
            xguide --> "X distance"
            yguide --> "Z distance"
        else
            xplt, yplt = points[:, 1], points[:, 2]
            xmplt, ymplt = elements[:, 1], elements[:, 2]
            xguide --> "X distance"
            yguide --> "Y distance"
        end

        @series begin
            seriestype := :path
            linewidth --> 4
            linecolor --> :black
            xplt, yplt
        end
        @series begin
            seriestype := :scatter
            markersize --> 6
            xplt, yplt
        end
        @series begin
            seriestype := :scatter
            markershape --> :diamond
            markersize --> 4
            xmplt, ymplt
        end
    end
end




"""
    update_forces!(distributed_loads, Fx, Fy, Mx, blade, assembly; fit=DS.linear)

Populate `distributed_loads` element-by-element from the current aero
loads at the blade's radial stations. For each GXBeam element, builds a
`GXBeam.DistributedLoads` whose `fy(s)` and `fz(s)` are 1D interpolants of
`Fy` and `-Fx` over `blade.r`, evaluated on the element's arc-length
sub-range `[s1, s2]` defined by the bracketing assembly points.

When `Fx` carries dual numbers (`ForwardDiff.Dual` or
`ReverseDiff.TrackedReal`), the loads are first run through
`GXBeam.dual_safe_copy!` to strip tracking before interpolation —
preserving the dual seed without confusing GXBeam's load builder.

Mutates `distributed_loads` in place. `Mx` is currently accepted but not
applied (bending moment coupling is not yet wired in).

**Arguments**
- `distributed_loads::Dict{Int, GXBeam.DistributedLoads}`: per-element load dict to fill.
- `Fx`, `Fy`: per-station aerodynamic force components along the blade.
- `Mx`: per-station moment (accepted but unused — see note above).
- `blade::Blade`: blade providing the radial-station vector `blade.r`.
- `assembly::GXBeam.Assembly`: assembly whose elements will receive loads.
- `fit`: 1D interpolant constructor taking `(snew, x, y)` signature (default `DS.linear`).
"""
function update_forces!(distributed_loads, Fx, Fy, Mx, blade, assembly; fit=DS.linear)

    #Todo: The problem here is that I'm creating a new interpolation struct every iteration... which is taking a significant amount of time. So I need to rethink how this is being done. If I'm going to use an arbitrary interpolation like this, then I need to pass the interpolation from timestep to time step and update it. Otherwise, I need to manually interpolate (without a functor). -> I think the function approach should work well here. 

    #todo: I think that this is a bit of a problem, because what if the rvec already includes rhub? (problem from before that needs to be resolved, see next Todo statement. ) #todo: I need to nail down behavior outside of the aero node regions. 
    #todo: I might need to extract the value out of Fx, Fy, and Mx if there is a tracked real present. 
    # @show eltype(Fx)
    if isa(Fx[1], ReverseDiff.TrackedReal)
        fz = zeros(length(Fx)) #TODO: Allocation every time step!!!
        fy = zeros(length(Fx))
        # mx = zeros(length(Fx))

        GXBeam.dual_safe_copy!(fz, -Fx)
        GXBeam.dual_safe_copy!(fy, Fy)
        # GXBeam.dual_safe_copy!(fx, Fx)

    elseif isa(Fx[1], ForwardDiff.Dual)
        fz = zeros(length(Fx))
        fy = zeros(length(Fx))
        # mx = zeros(length(Fx))

        GXBeam.dual_safe_copy!(fz, -Fx)
        GXBeam.dual_safe_copy!(fy, Fy)
        # GXBeam.dual_safe_copy!(fx, Fx)

    else
        fz = -Fx
        fy = Fy
        # Mxfit = fit(blade.r, Mx)
    end


    for ielem = eachindex(assembly.elements)
        # r1 = assembly.points[ielem][1] #Todo,: I want a vector of just lengths, of the points. Not just the X distance. 
        # r2 = assembly.points[ielem+1][1]
        r1 = norm(assembly.points[ielem])
        r2 = norm(assembly.points[ielem+1])
        # @show typeof(r1), typeof(r2) #Correct types

        # linear(xnew, x0, x1, y0, y1)
        # distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem; fy_follower = (s) -> fit(s, blade.r, fy), fz_follower = (s) -> fit(s, blade.r, fz), s1=r1, s2=r2) #, mx = (s) -> Mxfit(s) #Todo: Bending moment isn't coupled in!!!
        distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem; fy = (s) -> fit(s, blade.r, fy), fz = (s) -> fit(s, blade.r, fz), s1=r1, s2=r2) #, mx = (s) -> Mxfit(s) #Todo: Bending moment isn't coupled in!!!
        # distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(s), fz = (s) -> Fzfit(s), s1=r1, s2=r2) #, mx = (s) -> Mxfit(s)
        #todo: There is a slight problem here, if changing from follower loads to dead loads does absolutely nothing... then I'm not sure that what Taylor says they are doing is what they are actually doing. I need to look into that behavior. -> He applies the rotation matrix to the follower loads... And it looks like he does it correctly, or rather 

        # if ielem==length(assembly.elements)
        #     @show distributed_loads[ielem]
        # end
    end

end



"""
    get_blade_weight(assembly) -> Float64

Estimate the blade mass per unit length integrated along the assembly's
radial coordinate via the trapezoidal rule. Uses each element's
`mass[1,1]` entry (translational inertia along the first axis) as the
linear mass density at the element's midpoint.

**Arguments**
- `assembly::GXBeam.Assembly`: assembly with populated element mass matrices.

**Returns**
- Integrated blade mass (units follow the element mass-matrix convention).
"""
function get_blade_weight(assembly::GXBeam.Assembly)

    n = length(assembly.elements)

    rvec = zeros(n)
    masses = zeros(n)

    for i = 1:n
        rvec[i] = norm(assembly.elements[i].x)
        masses[i] = assembly.elements[i].mass[1,1]
    end

    return trapz(rvec, masses)
end