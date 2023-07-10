using Revise
using Rotors, OpenFASTsr, Plots, Statistics, DelimitedFiles, DynamicStallModels, CCBlade

#=
Verify that the BEMT results from the Rotors code (not just CCblade) is the same for OpenFAST, CCBlade core, and possibly QBlade.

    Adam Cardoza 6/29/23

=#

DS = DynamicStallModels
of = OpenFASTsr

localpath = @__DIR__
cd(localpath)

ofpath = "./" 

addriver = of.read_addriver("sn5_ADdriver.inp", ofpath)
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

hubht = addriver["HubHt(1)"]
pitch = edfile["BlPitch(1)"]
precone = edfile["PreCone(1)"]*pi/180 
yaw = edfile["NacYaw"]*(pi/180) 
tilt = edfile["ShftTilt"]*(pi/180)

rho = addriver["FldDens"]
mu = addriver["KinVisc"]*rho 
vinf = addriver["HWndSpeed_mat"][1]  
a = addriver["SpdSound"]  
omega = addriver["RotSpd_mat"][1]*2*pi/60 #Converting to rads/s
shearexp = addriver["PLExp_mat"][1]   
tsr = omega*rtip/vinf

env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)

n = length(adblade["BlSpn"])

if !@isdefined(readflag)
    readflag = true
end

if readflag
    println("Reading OpenFAST files...")

    fullouts = readdlm("./sn5_ADdriver.1.out", skipstart=6)

    names = fullouts[1,:]
    data = Float64.(fullouts[3:end,:])
    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]


    nt = length(tvec)

    global Uxmat = zeros(nt, n)
    global Uymat = zeros(nt, n)
    
    global phimat = zeros(nt, n)
    global alphamat = zeros(nt, n)

    for i = 1:n
        if i<10
            number = "00$i"
        elseif i<100
            number = "0$i"
        else
            number = "$i"
        end


        nameux = "AB1N"*number*"Vx" 
        nameuy = "AB1N"*number*"Vy"

        namealpha = "AB1N"*number*"Alpha"
        namephi = "AB1N"*number*"Phi"


        Uxmat[:,i] = outs[nameux]
        Uymat[:,i] = outs[nameuy]

        alphamat[:,i] = outs[namealpha]
        phimat[:,i] = outs[namephi]
    end

    readflag = false
end


assembly = of.make_assembly(edfile, bdfile, bdblade)

### Prep the ASD rotor and operating conditions 
aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput("../simpleNREL/Airfoils/Cylinder1.dat") 
aftypes[2] = of.read_airfoilinput("../simpleNREL/Airfoils/Cylinder2.dat") 
aftypes[3] = of.read_airfoilinput("../simpleNREL/Airfoils/DU40_A17.dat") 
aftypes[4] = of.read_airfoilinput("../simpleNREL/Airfoils/DU35_A17.dat") 
aftypes[5] = of.read_airfoilinput("../simpleNREL/Airfoils/DU30_A17.dat") 
aftypes[6] = of.read_airfoilinput("../simpleNREL/Airfoils/DU25_A17.dat") 
aftypes[7] = of.read_airfoilinput("../simpleNREL/Airfoils/DU21_A17.dat") 
aftypes[8] = of.read_airfoilinput("../simpleNREL/Airfoils/NACA64_A17.dat") 

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
rotor = Rotors.Rotor(Int(B), hubht, turbine; tilt, yaw, tip=CCBlade.PrandtlTipHub())

azimuth = 0.0
delta = zeros(3)
def_theta = zeros(3)
aerov = zeros(3)
t = 0.0
xcc = zeros(11)

ccouts = Vector{CCBlade.Outputs{eltype(rvec)}}(undef, n)


for i = 1:n
    Vx, Vy = Rotors.get_aerostructural_velocities(rotor, blade, env, t, i, azimuth, delta, def_theta, aerov)
    ccouts[i] = Rotors.solve_BEM!(rotor, blade, env, i, Vx, Vy, pitch, xcc)
end

nt, na = size(phimat)

### Compare inflow angles
phiplt = plot(xaxis="Radius (m)", yaxis="Inflow Angle (deg)")
plot!(rvec, phimat[end,:], lab="OpenFAST")
# for i = 1:nt
#     plot!(rvec, phimat[end,:], lab="OpenFAST $i")
# end
plot!(rvec, ccouts.phi.*(180/pi), lab="Rotors")
display(phiplt)

#=
Curious, the inflow angle is fairly off. I wonder why that is... and it looks like it could be significantly off. Such that it would drastically affect my solution. Possible differences:
- Convergence criteria
- what polar is getting used
- Inflow conditions (radius, velocities)
- Corrections? 

All the solutions across time are off, so that tells me there is something about OpenFAST that is causing it not to converge. 

- Convergence criteria in OpenFAST didn't match convergence criteria in CCBlade. 
=#

nothing