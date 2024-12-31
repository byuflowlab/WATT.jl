using Revise
using OpenFASTTools, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, DynamicStallModels
using YAML
using Infiltrator
using StaticArrays
# using Plots

DS = DynamicStallModels
of = OpenFASTTools

localpath = @__DIR__
cd(localpath)

# Infiltrator.clear_disabled!()
# Infiltrator.toggle_async_check(false)

### Read in OpenFAST files
ofpath = "/Users/adamcardoza/.julia/dev/Rotors/testing/SNC" 
inputfile = of.read_inputfile("sn5_input.fst", ofpath)
inflowwind = of.read_inflowwind("sn5_inflowwind.dat", ofpath)
# addriver = of.read_addriver("sn5_ADdriver.dvr", ofpath)
adfile = of.read_adfile("sn5_ADfile.dat", ofpath)
adblade = of.read_adblade("sn5_ADblade.dat", ofpath)
edfile = of.read_edfile("sn5_edfile.dat", ofpath)
bdfile = of.read_bdfile("sn5_bdfile.dat", ofpath)
bdblade = of.read_bdblade("sn5_bdblade.dat", ofpath)

B = edfile["NumBl"] #Number of Blades

rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
rfrac = bdblade["rfrac"]
chordvec = adblade["BlChord"]
twistvec = adblade["BlTwist"]

### Simulation corrections? 
rvec[1] += 0.001 
rtip = 63.0

hubht = inflowwind["RefHt"]
pitch = edfile["BlPitch(1)"]
precone = edfile["PreCone(1)"]*pi/180 
yaw = edfile["NacYaw"]*(pi/180) 
tilt = edfile["ShftTilt"]*(pi/180)

rho = inputfile["AirDens"]
mu = inputfile["KinVisc"]*rho
vinf = inflowwind["HWindSpeed"]
a = inputfile["SpdSound"]
omega = edfile["RotSpeed"]*2*pi/60 #Convert to rads/s
shearexp = inflowwind["PLexp"]
tsr = omega*rtip/vinf

env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)



n = length(adblade["BlSpn"])
ne = Int(bdblade["station_total"])

if !@isdefined(resetflag)
    resetflag = false
end

if resetflag
    readflag = true
    defflag = true
    runflag = true
end

if !@isdefined(readflag)
    readflag = true
end




assembly = of.make_assembly(edfile, bdfile, bdblade)

### Prep the ASD rotor and operating conditions 
aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput(ofpath*"/Airfoils/Cylinder1.dat") 
aftypes[2] = of.read_airfoilinput(ofpath*"/Airfoils/Cylinder2.dat") 
aftypes[3] = of.read_airfoilinput(ofpath*"/Airfoils/DU40_A17.dat") 
aftypes[4] = of.read_airfoilinput(ofpath*"/Airfoils/DU35_A17.dat") 
aftypes[5] = of.read_airfoilinput(ofpath*"/Airfoils/DU30_A17.dat") 
aftypes[6] = of.read_airfoilinput(ofpath*"/Airfoils/DU25_A17.dat") 
aftypes[7] = of.read_airfoilinput(ofpath*"/Airfoils/DU21_A17.dat") 
aftypes[8] = of.read_airfoilinput(ofpath*"/Airfoils/NACA64_A17.dat") 

# indices correspond to which airfoil is used at which station
af_idx = Int.(adblade["BlAFID"])


# create airfoil array
afs = aftypes[af_idx]

n = length(rvec)
airfoils = Vector{DS.Airfoil}(undef, n)
for i = 1:n
    airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
end

rR = rvec./rtip
# blade = Rotors.Blade(rhub, rtip, rR, airfoils)
blade = Rotors.Blade(rvec, twistvec.*(pi/180), airfoils; rhub=rhub, rtip=rtip, precone)

turbine = true
rotor_r = Rotors.Rotor(Int(B), hubht, turbine; tilt, yaw)



if !@isdefined(defflag)
    defflag = true
end

defflag = true
if defflag
    println("initialize...")
    aerostates, gxstates, mesh = initialize(blade, assembly; verbose=true)
    defflag = false
end

if !@isdefined(runflag)
    runflag = true
end


pitch = 0.0
azimuth0 = 0.0

runflag = true
if runflag #Ran up to 801
    aerostates, gxstates, mesh = Rotors.fixedpoint!(aerostates, gxstates, azimuth0, rotor_r, env,  blade, mesh, pitch; verbose=true, iterations=1)

    runflag = false 
end

# tipdx = [gxstates.points[i].u[1] for i = 1:ne]
# tipdy = [gxstates.points[i].u[2] for i = 1:ne]
# tipdz = [gxstates.points[i].u[3] for i = 1:ne]


# using Plots, LaTeXStrings

# loadplt = plot(xaxis="Blade fraction", yaxis="Load (N)", legend=(0.9, 0.3))
# plot!(loadplt, rfrac, aerostates.Fx, lab=L"$F_x$")
# plot!(loadplt, rfrac, aerostates.Fy, lab=L"$F_y$")
# display(loadplt)


# defplt = plot(xaxis="Blade fraction", yaxis="Deflection (m)")
# plot!(defplt, rfrac, tipdx, lab=L"$\delta_x$")
# plot!(defplt, rfrac, tipdy, lab=L"$\delta_y$")
# plot!(defplt, rfrac, tipdz, lab=L"$\delta_z$")
# display(defplt)


nothing