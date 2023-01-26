#=
Supplemental code from Dr. Ning on coordinate system transformations. 


=#


using FLOWMath: linear, trapz
using LinearAlgebra: norm
using SpecialFunctions: gamma

"""
Coordinate systems are defined here:
https://wisdem.readthedocs.io/en/master/wisdem/commonse/csystem.html
NOTE: I wrote those docs a long time ago, it is now mostly deprecated, 
we should create our own along these lines in our own docs.
All angles are in radians.
"""

# -------- Coordinate Systems ------------

"""
    CartesianVectors(x, y, z)
Convenience struct to group cartesian directions.  Typically each is a vector, but can be scalars.
"""
struct CartesianVectors{TVF}
    x::TVF
    y::TVF
    z::TVF
end

# Overload: add two vectors together.
Base.:+(v1::CartesianVectors, v2::CartesianVectors) = CartesianVectors(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)

# Overload: index one compponent in these vectors.
Base.getindex(v::CartesianVectors, i) = [v.x[i], v.y[i], v.z[i]]

"""
(private function)
Rotate coordinate system about +z in the +theta direction.
Right hand c.s. with x cross y = z.
all angles in radians
"""
function rotate_about_z(x, y, z, theta)

    c = cos.(theta)
    s = sin.(theta)

    xnew = @. x*c + y*s
    ynew = @. -x*s + y*c

    return xnew, ynew, z
end

function yaw_to_hub(vec, tilt)
    zh, xh, yh = rotate_about_z(vec.z, vec.x, vec.y, tilt)
    return CartesianVectors(xh, yh, zh)
end

function hub_to_yaw(vec, tilt)
    zy, xy, yy = rotate_about_z(vec.z, vec.x, vec.y, -tilt)
    return CartesianVectors(xy, yy, zy)
end

function hub_to_azimuth(vec, azimuth)
    yz, zz, xz = rotate_about_z(vec.y, vec.z, vec.x, azimuth)
    return CartesianVectors(xz, yz, zz)
end

function azimuth_to_hub(vec, azimuth)
    yh, zh, xh = rotate_about_z(vec.y, vec.z, vec.x, -azimuth)
    return CartesianVectors(xh, yh, zh)
end

function azimuth_to_blade(vec, cone)
    zb, xb, yb = rotate_about_z(vec.z, vec.x, vec.y, -cone)
    return CartesianVectors(xb, yb, zb)
end

function blade_to_azimuth(vec, cone)
    za, xa, ya = rotate_about_z(vec.z, vec.x, vec.y, cone)
    return CartesianVectors(xa, ya, za)
end

function blade_to_airfoil(vec, twist)
    xa, ya, za = rotate_about_z(vec.x, vec.y, vec.z, -twist)
    return CartesianVectors(xa, ya, za)
end

function airfoil_to_blade(vec, twist)
    xb, yb, zb = rotate_about_z(vec.x, vec.y, vec.z, twist)
    return CartesianVectors(xb, yb, zb)
end

# ---------------------------------------

"""
    Curvature(x, y, z, cone, s)
(private function) Internal representation of a blade with curvature.
Represented by coordinates in azimuthal-coordinate system.
cone is total coning angle including curvature.
s is total path length accounting for curvature.
"""
struct Curvature{TVF}
    x_az::TVF
    y_az::TVF
    z_az::TVF
    cone::TVF
    s::TVF
end

"""
    bladecurvature(precurve, presweep, r, precone)
Create a curvature definition for blade. 
# Arguments
- `precurve::Vector{Float64}`: x position in blade coordinate system.
- `presweep::Vector{Float64}`: y position in blade coordinate system.
- `r::Vector{Float64}`: z position in blade coordinate system.
- `precone::Float64`: precone angle
# Returns
- `curvature::Curvature`: internally used curvature object
"""
function bladecurvature(precurve, presweep, r, precone)

    # coordinate in azimuthal coordinate system
    az = CartesianVectors(precurve, presweep, r) |> 
        (x -> blade_to_azimuth(x, precone))

    # total coning angle 
    n = length(r)
    cone = zeros(n)
    cone[1] = atan(-(az.x[2] - az.x[1]), az.z[2] - az.z[1])
    cone[2:n-1] = atan.(-(az.x[3:n] - az.x[1:n-2]), az.z[3:n] - az.z[1:n-2])
    cone[n] = atan(-(az.x[n] - az.x[n-1]), az.z[n] - az.z[n-1])

    # total path length of blade
    s = zeros(n)
    s[1] = r[1]
    for i = 2:n
        s[i] = s[i-1] + norm(az[i] .- az[i-1])
          #sqrt((precurve[i] - precurve[i-1])^2 + 
            # (presweep[i] - presweep[i-1])^2 + (r[i] - r[i-1])^2)
    end

    return Curvature(az.x, az.y, az.z, cone, s)
end

"""
    Angles(yaw, tilt, azimuth, precone, twist, pitch)
Define orientation of rotor geometry
"""
struct Angles{TF}
    yaw::TF
    tilt::TF
    azimuth::TF
    precone::TF
    twist::Vector{TF}
    pitch::TF
end

"""
    AeroLoads(r, Np, Tp, Rp)
Aerodynamic loads (force per unit length) in blade-aligned coordinate system (N', T', R') -> (x, y, z)
Note that +y is in opposite direction of typical wind tuyrbine convention for T' so be sure to negate if needed.
"""
struct AeroLoads{TVF}
    r::TVF
    Np::TVF
    Tp::TVF
    Rp::TVF
end

"""
    loads(r, rhoA, aero, angles, curvature, Omega, g=9.81)
Combine aerodynamic, gravitational, and inertial (centrifugal) loading (latter is commented out since GXBeam already includes this).  
Output is in airfoil coordinate system.
# Arguments
- `r::Vector{Float}`: radial locations where total loads should be specified at (discretization may be different from aero)
- `rhoA::Vector{Float}`: mass per unit length (rho * A) along beam, at corresponding r locations
- `aero::AeroLoads`: aerodynamic loading per unit length in blade-aligned coordinate system
- `angles::Angles`: relevant angles for the wind turbine
- `curvature::Curvature`: internal representation of blade curvature
- `Omega::Float`: rotation speed (rad/s)
- `g::Float`: acceleration of gravity (for self-weight)
# Returns
- `total_loads::CartesianVectors`: loads in airfoil coordinate system.
"""
function loads(r, rhoA, aero, angles, curvature, Omega, g=9.81)  # assumes rhoA spaced with r

    lzero = zeros(length(r))

    # --- aerodynamic loads.  already in blade-aligned ----
    Px_aero = linear(aero.r, aero.Np, r)
    Py_aero = linear(aero.r, aero.Tp, r)  # user needs to input negative sign if wind turbine in CCBlade
    Pz_aero = linear(aero.r, aero.Rp, r)
    aero_loads = CartesianVectors(Px_aero, Py_aero, Pz_aero)

    # ----- gravity loads -----
    # yaw c.s.
    weight_yaw = CartesianVectors(lzero, lzero, -rhoA * g)

    # rotate yaw to blade
    gravity_loads = weight_yaw |> 
        (x -> yaw_to_hub(x, angles.tilt)) |> 
        (x -> hub_to_azimuth(x, angles.azimuth)) |> 
        (x -> azimuth_to_blade(x, curvature.cone))

    # # ------- centrifugal loads ---------
    # # azimuthal c.s.
    # centrifugal_az = CartesianVectors(lzero, lzero, rhoA * Omega^2 .* curvature.z_az)
    # centrifugal_loads = centrifugal_az |>
    #     (x -> azimuth_to_blade(x, curvature.cone))

    # ------- combine --------
    total_loads = aero_loads + gravity_loads |>  # + centrifugal_loads
        (x -> blade_to_airfoil(x, angles.twist .+ angles.pitch))

    return total_loads
end

"""
    tip_deflection(dx, dy, dz, angles, curvature)
Compute tip deflection in yaw coordinate system for purpose of computing tower strike clearance
# Arguments
- `dx::Vector{Float}`: x deflections of blade in airfoil c.s.
- `dy::Vector{Float}`: y deflections of blade in airfoil c.s.
- `dz::Vector{Float}`: z deflections of blade in airfoil c.s.
- `angles::Angles`: orientation of wind turbine
- `curvature::Curvature`: internal curved blade representation
# Returns
- `tip_deflection::Float`: tip defleciton in +x direction of yaw c.s.
"""
function tip_deflection(dx, dy, dz, angles, curvature)

    theta = angles.twist[end] + angles.pitch
    azimuth = pi  # closest to tower

    dr = CartesianVectors(dx[end], dy[end], dz[end])  # airfoil c.s.
    delta = dr |>
        (x -> airfoil_to_blade(x, theta)) |>
        (x -> blade_to_azimuth(x, curvature.cone[end])) |>
        (x -> azimuth_to_hub(x, azimuth)) |>
        (x -> hub_to_yaw(x, angles.tilt))

    return delta.x
end


"""
    compute_AEP(CDF, P, loss_factor)
AEP integration
# Arguments
- `CDF::Vector{Float}`: cumulatative distribution function
- `P::Vector{Float}`: corresponding powers
- `loss_factor::Float`: multicative factor for array losses, availability, etc.
# Returns
- `AEP::Float`: annual energy production
"""
function compute_AEP(CDF, P, loss_factor)
    return trapz(CDF, P)/1e3 * loss_factor * 365*24
end

"""
    weibull_cdf(x, A, k)
standard weibull distribution (CDF)
# Arguments
- `x::Vector{Float}`: argument of CDF
- `A::Float`: scale parameter (A > 0)
- `k::Float`: shape parameter (k > 0)
# Returns
- 'cdf::Vector{Float}`: corresponding CDF
"""
weibull_cdf(x, A, k) = @. 1.0 - exp(-(x/A)^k)

"""
    weibull_cdf(x, A, k)
weibull distribution based on mean value
# Arguments
- `x::Vector{Float}`: argument of CDF
- `xbar::Float`: mean value of distribution
- `k::Float`: shape parameter (k > 0)
# Returns
- 'cdf::Vector{Float}`: corresponding CDF
"""
function weibull_cdf_with_mean(x, xbar, k)
    A = xbar / gamma(1.0 + 1.0/k)
    return weibull_cdf(x, A, k)
end

"""
    rayleigh_cdf(x, A, k)
rayleigh distribution
# Arguments
- `x::Vector{Float}`: argument of CDF
- `xbar::Float`: mean value of distribution
# Returns
- 'cdf::Vector{Float}`: corresponding CDF
"""
function rayleigh_cdf(x, xbar)
    A = 2/sqrt(pi)*xbar
    return weibull_cdf(x, A, 2)
end

"""
    drivetrain_losses(Paero, Prated, dtype)
computer losses from drivetrain (based on NREL method)
In other words, convert from aerodynamic to mechanical power
# Arguments
- `Paero::Float`: aerodynamic power
- `Prated::Float`: rated power
- `dtype::String`: drivetrain type
# Returns
- 'Pmech::Vector{Float}`: mechanical power (after losses)
"""
function drivetrain_losses(Paero, Prated, dtype)

    if dtype == "GEARED"
        constant = 0.01289
        linear = 0.08510
        quadratic = 0.0

    elseif dtype == "SINGLE_STAGE"
        constant = 0.01331
        linear = 0.03655
        quadratic = 0.06107

    elseif dtype == "MULTI_DRIVE"
        constant = 0.01547
        linear = 0.04463
        quadratic = 0.05790

    elseif dtype == "PM_DIRECT_DRIVE"
        constant = 0.01007
        linear = 0.02000
        quadratic = 0.06899

    else
        error("invalid drivetrain type specified")
    end


    Pbar = Paero / Prated

    # handle negative power case (with absolute value)
    Pbar = abs(Pbar)

    # truncate idealized power curve for purposes of efficiency calculation
    Pbar = min(Pbar, 1.0)

    # compute mechanical losses
    loss = Prated * (constant + linear*Pbar + quadratic*Pbar^2)

    # mechanical power
    Pmech = Paero - loss
    Pmech = max(Pmech, 0.0)

    return Pmech
end



# # quick tests (to remove)

# Fx = zeros(5)
# Fy = zeros(5)
# Fz = [0.1, 0.2, 0.4, 0.3, 0.2]

# F1 = CartesianVectors(Fx, Fy, Fz)

# Fx = zeros(5)
# Fy = zeros(5)
# Fz = [0.2, 0.2, 0.3, 0.2, 0.2]

# F2 = CartesianVectors(Fx, Fy, Fz)

# println(F1 + F2)
# println(F1[2])

# precurve = zeros(5)
# presweep = zeros(5)
# r = [0.1, 0.2, 0.4, 0.7, 1.0]
# precone = 5*pi/180
# curvature = bladecurvature(precurve, presweep, r, precone)
# println(curvature.z_az)
# println(curvature.cone*180/pi)
# println(curvature.s)

# r = [0.1, 0.2, 0.4, 0.7, 1.0]
# rhoA = 0.05*ones(5)
# ra = [0.1, 1.0]
# Np = [2, 8.0]
# Tp = [2.0, 2]
# Rp = [0.0, 0.0]
# aero = AeroLoads(ra, Np, Tp, Rp)
# angles = Angles(5.0*pi/180, 4.0*pi/180, float(pi), precone, zeros(5), 10*pi/180)
# Omega = 10*pi/30

# total_loads = loads(r, rhoA, aero, angles, curvature, Omega)
# println(total_loads.x)
# println(total_loads.y)
# println(total_loads.z)

# dx = ones(5)
# dy = 0.1*ones(5)
# dz = zeros(5)
# defl = tip_deflection(dx, dy, dz, angles, curvature)
# println(defl)

# rho = 1.0
# V = range(0.0, 25, length=30)
# R = 63.0
# Cp = 0.4
# P = @. Cp * 0.5 * rho * V^3 * pi * R^2
# Prated = 5e6
# Vbar = 8.0
# k = 2.0
# cdf = weibull_cdf_with_mean(V, Vbar, k)
# dtype = "GEARED"
# Pmech = drivetrain_losses.(P, Prated, dtype)
# Pmech = min.(Pmech, Prated)
# loss_factor = 0.7
# compute_AEP(cdf, Pmech, loss_factor)

# using PyPlot
# pygui(true)
# close("all")
# figure()
# plot(V, Pmech/Prated)
# plot(V, cdf)