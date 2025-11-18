

export Rotor, Blade



"""
    Rotor(B, hubht, tilt, yaw, turbine, mach, re, rotation, tip)

A struct to hold information about the rotor. Modeled after CCBlade.Rotor. 

### Inputs
- B::Int - the number of blades. 
- hubht::Float - the height of the hub of the blade in the global reference frame. 
- tilt::Float - the tilt angle (radians). Positive about the negative X axis (tilts the turbine up). 
- yaw::Float - the yaw angle (radians). Positive about the negative Y axis (turns the turbine clockwise)
- turbine::Bool - Whether the rotor is a wind turbine or something else. 
- mach::CCBlade.MachCorrection - the mach number correction model (or nothing).
- re::CCBlade.ReCorrection - the Reynolds number correction model (or nothing). 
- rotation::CCBlade.RotationCorrection - the rotation correction model (or nothing).
- tip::CCBlade.TipCorrection - the tip correction model (or nothing). 
"""
struct Rotor{TI, TF, TB, T1 <: Union{Nothing, CCBlade.MachCorrection}, T2 <: Union{Nothing, CCBlade.ReCorrection}, T3 <: Union{Nothing, CCBlade.RotationCorrection}, T4 <: Union{Nothing, CCBlade.TipCorrection}}
    B::TI
    hubht::TF
    tilt::TF
    yaw::TF
    turbine::TB
    mach::T1
    re::T2
    rotation::T3
    tip::T4
end

"""
    Rotor(B, hubht, turbine; tilt=0.0, yaw=0.0, pitch=0.0, mach=nothing, re=nothing, rotation=nothing, tip=nothing)

A convinience constructor for the Rotor struct. 

### Inputs
- B::Int - number of blades
- hubht::Float - The height of the hub from the ground (meters).
- turbine::Bool - Whether the rotor is a turbine or not. 
- tilt::Float - the tilt angle (radians). Positive about the negative X axis (tilts the turbine up). 
- yaw::Float - the yaw angle (radians). Positive about the negative Y axis (turns the turbine clockwise)
- mach::CCBlade.MachCorrection - the mach number correction model (or nothing).
- re::CCBlade.ReCorrection - the Reynolds number correction model (or nothing). 
- rotation::CCBlade.RotationCorrection - the rotation correction model (or nothing).
- tip::CCBlade.TipCorrection - the tip correction model (or nothing). 
"""
function Rotor(B, hubht, turbine; tilt=0.0, yaw=0.0, mach=nothing, re=nothing, rotation=nothing, tip=nothing)
    return Rotor(B, hubht, tilt, yaw, turbine, mach, re, rotation, tip)
end


"""
    Blade(rhub, rtip, tilt, yaw, rx, ry, rz, twist, airfoils)

The blade struct... what else is there to say. 

### Inputs
- rhub::TF - the radial distance from the center of rotation to the edge of the hub (where the blade begins). (meters)
- rtip::TF - the blade length distance from the center of rotation to the tip of the blade. (i.e. rtip-rhub should be the total length of the blade). (meters)
- rx::AbstractVector{<:TF} - The x position of the blade in the rotor plane (not including tilt and yaw, but including precone, curve, and sweep). 
- ry::AbstractVector{<:TF} - The y position of the blade in the rotor plane (not including tilt and yaw, but including precone, curve, and sweep). 
- rz::AbstractVector{<:TF} - The z position of the blade in the rotor plane (not including tilt and yaw, but including precone, curve, and sweep). 
- rR::AbstractVector{<:TF} - The percent of the blade length that the aerodynamic node is located at. 
- thetax::AbstractVector{<:TF} - The angular position of the local coordinate frame from the blade root reference frame about the X axis. Also known as the sweep angle.  
- thetay::AbstractVector{<:TF} - The angular position of the local coordinate frame from the blade root reference frame about the Y axis. Also known as the curve angle.  
- twist::AbstractVector{<:TF} - the twist distribution (radians). 
- airfoils::AbstractVector{<:DynamicStallModels.Airfoil} - The vector of the airfoil structs (an airfoil for each node). 
"""
# struct Blade{TF, TVF, TVF2, TVDS}  
#     rhub::TF
#     rtip::TF
#     rx::TVF #Lead-lag direction (freestream) curve value
#     ry::TVF #Flapwise direction (Sweep value)
#     rz::TVF #Radial direction
#     r::TVF #Todo: Decide if I want this in here, or if I'll just use the norm of the vectorized version. 
#     rR::TVF
#     c::TVF2 #Chord length
#     thetax::TVF #Sweep angle
#     thetay::TVF #Curve angle
#     twist::TVF2
#     precone::TF
#     xcp::TVF 
#     airfoils::TVDS
# end
struct Blade{TF, TF2}  
    rhub::TF
    rtip::TF
    rx::AbstractVector{<:TF} #Lead-lag direction (freestream) curve value
    ry::AbstractVector{<:TF} #Flapwise direction (Sweep value)
    rz::AbstractVector{<:TF} #Radial direction
    r::AbstractVector{<:TF} #Todo: Decide if I want this in here, or if I'll just use the norm of the vectorized version. 
    rR::AbstractVector{<:TF}
    c::AbstractVector{<:TF2} #Chord length
    thetax::AbstractVector{<:TF} #Sweep angle
    thetay::AbstractVector{<:TF} #Curve angle
    twist::AbstractVector{<:TF2}
    precone::TF
    xcp::AbstractVector{<:TF} 
    airfoils::AbstractVector{<:DS.Airfoil}
end

"""
    Blade(rvec, twist, airfoils; rhub=rvec[1], rtip=rvec[end], precone=0.0, sweep=0.0, curve=0.0, rx=zero(rvec), ry=zero(rvec)) -> blade

A convinience constructor for the blade struct. 

**Arguments**
- `span::AbstractVector{<:Float}`: The radial distances of the aerodynamic nodes. 
- `twist::AbstractVector{<:Float}`
- airfoils::AbstractVector{<:DynamicStallModels.Airfoil}
- rhub::Float - the radial distance to the hub. 
- rtip::Float - the total length of the blade. 
- precone::Float - The precone angle of the entire blade.
- sweep::Union{Float, AbstractVector{Float}} - The local sweep angle
- curve::Union{Float, AbstractVector{Float}} - The local curve angle. 
- `rx::Float`: The curve distances, or the distances in the X direction of the aerodynamic nodes.
- `ry::Float`: The sweep distances, or the distances in the Y direction of the aerodynamic nodes.
"""
function Blade(span, chord, twist, xcp, airfoils::AbstractVector{<:DS.Airfoil}; rhub=span[1], rtip=span[end], precone=0.0, sweep=0.0, curve=0.0, rx=zero(span), ry=zero(span))
    n = length(airfoils)

    if length(span)!=length(twist)!=n
        error("Blade(): The number of airfoils and nodes (span) but be the same.")
    end

    if length(sweep)==1
        sweep = sweep.*ones(n)
    end

    if length(curve)==1
        curve = curve.*ones(n)
    end

    if length(sweep)!=length(curve)!=n
        error("Blade(): The length of the sweep and curve vectors must be as long as the radial node vector.")
    end

    pi2 = pi/2

    if precone>pi2
        @warn("Blade(): The precone angle you provided appears to be degrees, not radians.")
    elseif any(i->i>pi2, sweep)
        @warn("Blade(): The sweep angle(s) you provided appears to be degrees, not radians.")
    elseif any(i->i>pi2, curve)
        @warn("Blade(): The curve angle(s) you provided appears to be degrees, not radians.")
    end

    if any(i->i>pi2, twist)
        @warn("Blade(): The twist angle appears to be in degrees, not radians.")
    end

    # rx = zero(span)
    # ry = zero(span)
    # rz = zero(span)
    rvec = @. sqrt(rx^2 + ry^2 + span^2)
    rR = span./rtip

    # for i in eachindex(airfoils)
    #     rx[i] = span[i]*sin(curve[i])
    #     ry[i] = span[i]*sin(sweep[i])
    #     rz[i] = span[i]*cos(curve[i])*cos(sweep[i])
    #     rx[i], ry[i], rz[i] = rotate_vector(rx[i], ry[i], rz[i], 0, precone, 0)
    # end

    return Blade(rhub, rtip, rx, ry, span, rvec, rR, chord, sweep, curve,  twist, precone, xcp, airfoils)
end






# struct AeroStates{TF1, TF2, TA} #Todo: Parametric typing
#     azimuth::TF1 #Float
#     phi::TF2 #These are arrays of floats, all the same size     #Inflow angle
#     alpha::TF2   #Angle of attack
#     W::TF2       #Inflow velocity
#     cx::TF2
#     cy::TF2
#     cm::TF2
#     fx::TF2 #Todo: I should get rid of either the loads or the coefficients. I don't think that I need to store both... I think that I can have a function to fetch on based on the other. 
#     fy::TF2
#     mx::TF2
#     xds::TA #This is an array of floats, but different size than the rest. #Todo: I might be able to get rid of this. Do I really need the intermediate dynamic stall states? That would cut down my allocations by alot.  
# end

# function get_aerostate(aerostates::AeroStates, i)
#     nt = length(aerostates.azimuth) #Get the total number of time steps

#     #check and see if the requested aerostate is in the vectors. 
#     if i>nt
#         error("get_aerostates(): There aren't that many time steps.")
#     end

#     azimuth = aerostates.azimuth[i]
#     phi = aerostates.phi[i,:]
#     alpha = aerostates.alpha[i,:]
#     W = aerostates.W[i,:]
#     cx = aerostates.cx[i,:]
#     cy = aerostates.cy[i,:]
#     cm = aerostates.cm[i,:]
#     fx = aerostates.fx[i,:]
#     fy = aerostates.fy[i,:]
#     mx = aerostates.mx[i,:]
#     xds = aerostates.xds[i,:]

#     return azimuth, phi, alpha, W, cx, cy, cm, fx, fy, mx, xds
# end



# """
#     Mesh(interpolationpoints, delta, def_theta, aerov, xcc, xds_idxs, p_ds, 
#             system, prescribed_conditions, distributed_loads, paug, constants)

# **Fields**
# - interpolationpoints::Vector{InterpolationPoint}
# - delta::Vector{StaticArray(Float64, Float64, Float64)} - The linear structural deflection interpolated to the aerodynamic nodes.
# - def_theta::Vector{StaticArray(Float64, Float64, Float64)} - The angular structural deflection interpolated to the aerodynamic nodes.
# - aerov::Vector{StaticArray(Float64, Float64, Float64)} - The structural linear velocities interpolated to the aerodynamic nodes. 
# - xxc::Vector{Float64} - 
# - xds_idxs
# - p_ds
# - system
# - distributed_loads
# - paug
# - constants::NamedTuple 
# """
# struct Mesh{TIP, TSA, TAF, TAF2, TAI, TT, TS}
#     interpolationpoints::TIP
#     delta::TSA #The structural deflection at the aerodynamics nodes. 
#     def_theta::TSA
#     aerov::TSA
#     xcc::TAF
#     xds_idxs::TAI
#     p_ds::TAF
#     #GXBeam memory
#     assembly::TASS
#     system::TS
#     pcond::TPC #prescribed_condtions
#     dload::TDL #distributed_loads
#     lv::TLV #linear_velocity
#     av::TAV #angular_velocity
#     xpfunc::TXP
#     pfunc::TP
#     td::Bool
#     sd::Bool
# end

