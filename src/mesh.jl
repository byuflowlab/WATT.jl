




"""
    find_point_indices(rvec, r)

Find the indices of the two nearest neighboring nodes in a 1D mesh of a given position. Note that this function assumes that the mesh is monatonically increasing.

**Arguments**
- `rvec::Vector{Number}`: The 1D mesh of interest.
- r - the position of the point of interest

**Returns**
- pair - a tuple containing the indices of the point before and after the aerodynamic point (before, after)
"""
function find_point_indices(rvec, r)

    n1 = length(rvec)

    if r<rvec[1]
        @warn("find_point_indices(): The point was not found in the 1D mesh.")
        return 1, 2
    elseif r>rvec[end]
        @warn("find_point_indices(): The point was not found in the 1D mesh.")
        return (n1-1, n1)
    end

    ### Iterate through the points 
    for i = 1:n1-1
        if rvec[i]>r
            return (i-1, i)
        end
    end

    return (n1-1, n1)
end


"""
    find_interpolation_percent(rvec, pair, r)

Find the percentage location of the point of interest between the neighboring nodes.

**Inputs**
- `rvec::Vector{Number}`: The 1D mesh of interest. 
- `pair::Tuple{Int, Int}`: The indicial pair of the surrounding neighbors.
- `r::Number`: The point of interest (on the second mesh... not the mesh of interest)

**Outputs**
- `percent::Number`: The percentage location between the two neighboring nodes.
"""
function find_interpolation_percent(rvec, pair, r)

    a = r - rvec[pair[1]]
    L = rvec[pair[2]]- rvec[pair[1]]
    return a/L
end

"""
    InterpolationPoint(pair, percent)

A struct to quickly interpolate from one mesh to another. . 
"""
struct InterpolationPoint{TF}
    pair::Tuple{Int64, Int64}
    percent::TF
end

#=
Note: Maybe down the line I can switch InterpolationPoint to a type, and make different structs of linear, quadratic, cubic, or maybe even more general of polynomial, Lagrange, Legendre and interpolate that way.

=#

"""
    create_interpolationpoints(assembly, rvec)

**Inputs**
- assembly - A GXBeam assembly
- rvec - The radial location of the aerodynamic nodes

**Outputs**
- `interpolationpoints::Vector{InterpolationPoint}`

### Notes
- This assumes that the structural points are in order and go from root to tip. Additionally, the aerodynamic points follow a similar order. 
"""
function create_interpolationpoints(assembly::GXBeam.Assembly, blade::Blade)

    ### Allocate memory
    na = length(blade.r)
    points = Vector{InterpolationPoint}(undef, na)

    ### Find the beam length 
    rgx = get_bladelength_vector(assembly)

    for i = 1:na
        pair = find_point_indices(rgx, blade.r[i])
        percent = find_interpolation_percent(rgx, pair, blade.r[i])
        points[i] = InterpolationPoint(pair, percent)
    end

    return points
end


"""
    interpolate_position(ip, assembly, state)

Get the position of aerodynamic nodes from the interpolation struct 
and the structural structs. 
"""
function interpolate_position(ip, assembly, state)  

    p1 = assembly.points[ip.pair[1]] + state.points[ip.pair[1]].u
    p2 = assembly.points[ip.pair[2]] + state.points[ip.pair[2]].u

    return (1-ip.percent)*p1 + ip.percent*p2
end

"""
    interpolate_deflection(ip, assembly, state)

Interpolate the structural deflection onto the aerodynamic mesh. 

**Arguments**
- `ip::InterpolationPoint`: The interpolation point of interest.
"""
function interpolate_deflection(ip, assembly, state)
    p1 = state.points[ip.pair[1]].u
    p2 = state.points[ip.pair[2]].u

    return (1-ip.percent)*p1 + ip.percent*p2
end

"""
    interpolate_velocity(ip, assembly, state, env)

Get the relative velocity of an aerodynamic point. 
"""
function interpolate_velocity(ip, assembly, state) 

    v1 = state.points[ip.pair[1]].V
    v2 = state.points[ip.pair[2]].V

    return (1-ip.percent)*v1 + ip.percent*v2
end

function interpolate_angle(ip, assembly, state) 
    
    theta1 = WMPtoangle(state.points[ip.pair[1]].theta)
    theta2 = WMPtoangle(state.points[ip.pair[2]].theta) 

    return (1-ip.percent)*theta1 + ip.percent*theta2
end

"""
    convert_velocities(blade, env, assembly, state, interpolationpoints, t, idx)

Interpolate the aerostructural velocities, transform them into the aerodynamic reference frame, and remove the rotational portion.

**Arguments**

"""
function convert_velocities(blade::Blade, env::Environment, assembly, state, interpolationpoints, t, idx)

    ### Interpolate the structural quantities to the aerodynamic mesh. 
    delta = interpolate_deflection(interpolationpoints[idx], assembly, state)
    Vs = interpolate_velocity(interpolationpoints[idx], assembly, state)
    

    ### Transform the quantities into the aerodynamic frame
    # dx = -delta[3]
    dy = delta[2]
    dz = delta[1]
    usx = Vs[3]
    usy = -Vs[2]
    usz = -Vs[1]
    #Note: The  velocities are negative of the transformation because as the structure moves in those directions, the wind moves in the opposite direction. 

    ### Remove the rotational velocities
    # rax = blade.rx[idx] 
    ray = blade.ry[idx] + dy
    raz = blade.rz[idx] + dz

    

    omega = env.RS(t)

    #Rotational velocities
    # urx = 0.0
    ury = omega*raz 
    urz = -omega*ray

    return (usx, usy-ury, usz-urz)
end

function update_mesh!(blade::Blade, mesh::Mesh, assembly::GXBeam.Assembly, gxstate, env::Environment, t, na)

    for j = 1:na
        mesh.delta[j] = interpolate_deflection(mesh.interpolationpoints[j], assembly, gxstate)

        mesh.def_theta[j] = interpolate_angle(mesh.interpolationpoints[j], assembly, gxstate)

        mesh.aerov[j] = convert_velocities(blade, env, assembly, gxstate, mesh.interpolationpoints, t, j)
    end
end



function transform_BC_G(rhr_x, rhr_y, rhr_z, azimuth, precone, tilt, yaw)

    rtx, rty, rtz = rotate_y(rhr_x, rhr_y, rhr_z, precone; T=false)
    rtx, rty, rtz = rotate_x(rtx, rty, rtz, azimuth; T=true)
    rtx, rty, rtz = rotate_y(rtx, rty, rtz, tilt; T=true)
    rg_x, rg_y, rg_z = rotate_z(rtx, rty, rtz, yaw; T=true)

    return rg_x, rg_y, rg_z
end

function transform_BC_HR(rbc_x, rbc_y, rbc_z, precone)

    rhr_x, rhr_y, rhr_z = rotate_y(rbc_x, rbc_y, rbc_z, precone; T=false)

    return rhr_x, rhr_y, rhr_z
end

function transform_HR_L(rhr_x, rhr_y, rhr_z, curve, sweep, precone)

    rt_x, rt_y, rt_z = rotate_x(rhr_x, rhr_y, rhr_z, sweep; T=true)
    rt_x, rt_y, rt_z = rotate_y(rt_x, rt_y, rt_z, curve; T=true)
    rl_x, rl_y, rl_z = rotate_y(rt_x, rt_y, rt_z, precone; T=false)

    return rl_x, rl_y, rl_z
end

function transform_G_L(rg_x, rg_y, rg_z, azimuth, curve, precone, sweep, tilt, yaw)

    rt_x, rt_y, rt_z = rotate_z(rg_x, rg_y, rg_z, yaw)
    rt_x, rt_y, rt_z = rotate_y(rt_x, rt_y, rt_z, tilt)
    rt_x, rt_y, rt_z = rotate_x(rt_x, rt_y, rt_z, azimuth)
    rt_x, rt_y, rt_z = rotate_x(rt_x, rt_y, rt_z, sweep; T=true)
    rt_x, rt_y, rt_z = rotate_y(rt_x, rt_y, rt_z, curve; T=true)
    rl_x, rl_y, rl_z = rotate_y(rt_x, rt_y, rt_z, precone; T=true)

    return rl_x, rl_y, rl_z
end