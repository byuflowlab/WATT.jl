# Getting Started

Here we will go through a basic unsteady wind turbine blade analysis of a simplified version of the NREL 5MW wind turbine. 

## Preparing inputs
In order to simulate we need to set up the aerodynamic mesh, the structural mesh, and the environmental variables for simulation. 

### Aerodynamic mesh
Inside of the aerodynamic mesh, there are two key objects that we need. The `Blade` and the `Rotor`. We'll first go through building up the blade.

#### Blade object
The core of the blade object is a vector that defines a reference line spanning from the root to the tip. This line is defined relative to the center of rotation (starting from the hub radius). It can change in all 3 directions. The blade is described by a distribution of chord values, twist values (radians), the airfoil reference point (i.e. where the chord intersects the blade reference line), and dynamic airfoil polars along the reference line. 

Since the other distributions are fairly straightforward, we start with the airfoil polar distribution. The airfoil polar is a dynamic airfoil polar from [DynamicStallModels](https://github.com/byuflowlab/DynamicStallModels.jl.git). You can see more about generating your own custom dynamic airfoil polar in the documentation for that package, here we use the package [OpenFASTTools](https://github.com/byuflowlab/OpenFASTTools.jl.git) to read in OpenFAST airfoil files. 

```julia
using DynamicStallModels, OpenFASTTools

DS = DynamicStallModels
of = OpenFASTTools

#Read in the different airfoils used in the NREL 5MW
aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput("./Airfoils/Cylinder1.dat") 
aftypes[2] = of.read_airfoilinput("./Airfoils/Cylinder2.dat") 
aftypes[3] = of.read_airfoilinput("./Airfoils/DU40_A17.dat") 
aftypes[4] = of.read_airfoilinput("./Airfoils/DU35_A17.dat") 
aftypes[5] = of.read_airfoilinput("./Airfoils/DU30_A17.dat") 
aftypes[6] = of.read_airfoilinput("./Airfoils/DU25_A17.dat") 
aftypes[7] = of.read_airfoilinput("./Airfoils/DU21_A17.dat") 
aftypes[8] = of.read_airfoilinput("./Airfoils/NACA64_A17.dat") 

# indices correspond to which airfoil is used at which station
af_idx = [3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8]

# create airfoil array
afs = aftypes[af_idx]

n = length(af_idx)
airfoils = StructArray{DS.Airfoil}(undef, n)
xcp = Vector{Float64}(undef, n)
for i = 1:n
    airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
end
```

Here, we use a helper constructor to simplify generating the reference line and distribution. We simply need a vector of radial locations that is the same length of the distribution of other values. 

```julia

rvec = [1.501, 4.0625, 6.625, 9.1875, 11.75, 14.3125, 16.875, 19.4375, 22.0, 24.5625, 27.125, 29.6875, 32.25, 34.8125, 37.375, 39.9375, 42.5, 45.0625, 47.625, 50.1875, 52.75, 55.3125, 57.875, 60.4375, 62.9999] #meters

chordvec = [3.542, 3.6782526, 3.9713763, 4.2647549, 4.557, 4.6896847, 4.6057022, 4.4836912, 4.3556258, 4.2205945, 4.0700368, 3.9094362, 3.748, 3.59425, 3.4405, 3.28675, 3.133, 2.97925, 2.8255, 2.67175, 2.518, 2.3642514, 2.1822262, 1.7173923, 1.4190122] #meters

twistvec = [13.308, 13.308, 13.308, 13.308, 13.308, 12.257711, 11.137553, 10.311991, 9.5908574, 8.8603117, 8.1019334, 7.3238479, 6.544, 5.8027373, 5.0666879, 4.3341462, 3.6145353, 3.0214745, 2.5174169, 2.0216489, 1.526, 1.0287878, 0.54745436, 0.18292632, 0.10600483].*(pi/180) #radians

rhub = 1.5
rtip = 63.0
precone = 0.0 

blade = WATT.Blade(rvec, chordvec, twistvec, xcp, airfoils; rhub=rhub, rtip=rtip, precone)
```

#### Rotor object
The rotor object is significantly simpler:

```julia
B = 3 #Number of blades, integer
hubht = 90.0 #Height of the center of the hub, meters
tilt = 0.0 #Tilt of the rotor plane, radians
yaw = 0.0 #Yaw angle of the rotor plane, radians
turbine = true

rotor = WATT.Rotor(B, hubht, turbine; tilt, yaw)
```


### Structural mesh
Now that we have an aerodynamic mesh, we need a structural mesh. The structural mesh should follow the same reference line as the aerodynamic mesh. We directly use the `Assembly` object from [GXBeam](https://github.com/byuflowlab/GXBeam.jl.git). Creating the structural mesh requires knowledge of the cross-sectional properties, specifically mass and compliance properties. The package [GXBeamCS](https://github.com/byuflowlab/GXBeamCS.git) is helpful in obtaining those properties. To help us get up to speed faster with using WATT, we use `OpenFASTTools` again to read in a `BeamDyn` file. 

```julia
### Read in the dictionaries describing the OpenFAST files
ofpath = "./" 
edfile = of.read_edfile("sn5_edfile.dat", ofpath)
bdfile = of.read_bdfile("sn5_bdfile.dat", ofpath)
bdblade = of.read_bdblade("sn5_bdblade.dat", ofpath)

### Create the assembly
assembly = of.make_assembly(edfile, bdfile, bdblade)
```


### Environmental variables
Finally, we can create an environment that provides the wind speed as a function of time and the rotation rate of the wind turbine. Here, we let the incoming wind speed be a linear combination of sinusoids. 

```julia
windspeedfun(t) = 10.0 + 2*sin(2*pi*t/2) + 0.5*sin(2*pi*t/0.1)

ufun = (t) -> SVector(windspeedfun(t), 0.0, 0.0) #Velocity of a uniform wind field, m/s
omegafun = (t)-> SVector(0.0, 0.0, 0.0) #Swirling component of velocity field
udotfun = (t)-> SVector(0.0, 0.0, 0.0) #Acceleration of uniform wind field (needed for some stall models)
omegadotfun = (t) -> SVector(0.0, 0.0, 0.0) #Angular acceleration of the velocity field
Vinf = (t) -> windspeedfun(t) #Velocity magnitude of the windfield
Vinfdot = (t) -> 0.0
RS = (t) -> omega #Rotation rate of the wind turbine, radians
RSdot = (t) -> 0.0
env = Rotors.SimpleEnvironment(rho, mu, a, shearexp, ufun, omegafun, udotfun, omegadotfun, Vinf, RS, Vinfdot, RSdot)
```

## Simulating

Now we're ready to initialize the memory for the states and create the connectivity between the meshes. The connectivity is stored in the `mesh` named tuple. We'll need the amount of time so we can create large enough storage. 

```julia
tvec = 0:0.01:100

aerostates, gxhistory, mesh = WATT.initialize_sim(blade, assembly, tvec; verbose=true)
```

And we can simulate the aerostructural response of the wind turbine across time: 

```julia
WATT.run_sim!(rotor, blade, mesh, env, tvec, aerostates, gxhistory; verbose=true)
```

Now that the simulation is done, we can grab what values we want from the aerodynamic or structural state storage containers. 

### Aerodynamic states
```julia
# julia> keys(aerostates)
# (:azimuth, :phi, :alpha, :W, :Cx, :Cy, :Cm, :Fx, :Fy, :Mx, :xds)
using Plots, LaTeXStrings

tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)")
plot!(tiploads, tvec_r, aerostates.Fx[:,end], lab=L"$F_x$")
plot!(tiploads, tvec_r, aerostates.Fy[:,end], lab=L"$F_y$")
display(tiploads)
```


### Structural states
```julia
nt = length(tvec)
tipdef_x = [gxhistory[i].points[end].u[1] for i in 1:nt]
tipdef_y = [gxhistory[i].points[end].u[2] for i in 1:nt]
tipdef_z = [gxhistory[i].points[end].u[3] for i in 1:nt]


tiptheta_xof = zeros(nt)
tiptheta_yof = zeros(nt)
tiptheta_zof = zeros(nt)

for i = 1:nt
    theta = WATT.WMPtoangle(gxhistory[i].points[end].theta) #convert the Wiener-Milenkovic parameters to Euler angles. 
    tiptheta_x[i] = theta[1]
    tiptheta_y[i] = theta[2]
    tiptheta_z[i] = theta[3]
end

tipdefs = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)")
plot!(tvec_r, -tipdef_z, lab=L"$\delta x$") #Plot in OpenFAST reference frame
plot!(tvec_r, tipdef_y, lab=L"$\delta y$")
plot!(tvec_r, tipdef_x, lab=L"$\delta z$")
display(tipdefs)
```