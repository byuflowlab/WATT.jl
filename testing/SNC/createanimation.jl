#=
Compare the take_step function for gxbeam to the time_domain_analysis. 

=#

using Revise
using OpenFASTTools, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, DynamicStallModels
using StaticArrays
# using Plots

DS = DynamicStallModels
of = OpenFASTTools

localpath = @__DIR__
cd(localpath)


### Read in OpenFAST files
ofpath = "./" 
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




assembly = of.make_assembly(edfile, bdfile, bdblade)

### Prep the ASD rotor and operating conditions 
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
af_idx = Int.(adblade["BlAFID"])


# create airfoil array
afs = aftypes[af_idx]

n = length(rvec)
airfoils = Vector{DS.Airfoil}(undef, n)
for i = 1:n
    airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
end

rR = rvec./rtip
blade = Rotors.Blade(rvec, twistvec.*(pi/180), airfoils; rhub=rhub, rtip=rtip, precone)

turbine = true
rotor_r = Rotors.Rotor(Int(B), hubht, turbine; tilt, yaw)



if !@isdefined(defflag)
    defflag = true
end

# tvec = [0.000, 0.001]
tvec = 0:0.01:2
# tvec = 0:0.01:9.6

# if defflag
# println("initialize...")
aerostates, gxstates, gxhistory, mesh = initialize(blade, assembly, tvec; verbose=true)
    # defflag = false
# end

if !@isdefined(runflag)
    runflag = true
end


runflag = true
if runflag
    aerostates, gxstates, gxhistory = simulate!(rotor_r, env, tvec, aerostates, gxstates, gxhistory, blade, mesh; verbose=true)

    runflag = false 
end


tipdef_x = [gxhistory[i].points[end].u[1] for i in eachindex(tvec)]

# @show tipdef_x

include("randomairfoil.jl")

airfoil = airfoil.*5

section = zeros(3, size(airfoil, 1))
for ic = 1:size(airfoil, 1)
    section[1,ic] = airfoil[ic,1] - 0.5
    section[2,ic] = 0
    section[3,ic] = airfoil[ic,2]
end



mkpath("dynamic-simulation")
write_vtk("dynamic-simulation/dynamic-simulation", assembly, gxhistory, tvec, scaling=5;
    sections = section)

#Note: Correct solution. (Or rather the solution I was getting with time_domain_analysis) #tvec = 0:0.001:0.01
#[0.0, 3.9924964120256665e-5, 0.00015876402804123688, 0.0003542923417708522, 0.0006238931195408572, 0.0009641344237213083, 0.0013719055459919265, 0.001843945046587685, 0.002376257741316404, 0.0029647096623815136, 0.0036052055288560625]

nothing