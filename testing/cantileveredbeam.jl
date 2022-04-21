using GXBeam, Plots, StaticArrays

## Turbine information
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

## Analysis information
np = 21 #Number of points
ne = 20 #Number of elements





## Create points
z = range(rhub, rtip, length=np)

create_gxbeam_point(p) = SVector(p[1], p[2], p[3])

points = [create_gxbeam_point([0.0, 0.0, z[ip]]) for ip = 1:np]

## Starting and stoping point indices
start = 1:np-1
stop = 2:np

## Compliance matrix

## Mass Matrix
