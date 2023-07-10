using Revise
using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, DynamicStallModels
using YAML
using Infiltrator
using StaticArrays
# using Plots

DS = DynamicStallModels
of = OpenFASTsr

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
mu = inputfile["KinVisc"]
vinf = inflowwind["HWindSpeed"]
a = inputfile["SpdSound"]
omega = edfile["RotSpeed"]*2*pi/60 #Convert to rads/s
shearexp = inflowwind["PLexp"]
tsr = omega*rtip/vinf

env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)



n = length(adblade["BlSpn"])
ne = Int(bdblade["station_total"])

if !@isdefined(readflag)
    readflag = true
end

if readflag
    println("Reading OpenFAST files...")
    # fullouts = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=6)
    fullouts = readdlm("./sn5_input.out", skipstart=6)

    names = fullouts[1,:]

    # data = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=8)
    # data = readdlm("./simpleNREL/sn5_input.out", skipstart=8)
    data = Float64.(fullouts[3:end,:])

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]


    nt = length(tvec)


    global fxmat = zeros(nt, n)
    global fymat = zeros(nt, n)
    # global Uxmat = zeros(nt, n)
    # global Uymat = zeros(nt, n)
    global Mmat = zeros(nt, n)
    # global alphamat = zeros(nt, n)

    for i = 1:n
        if i<10
            number = "00$i"
        elseif i<100
            number = "0$i"
        else
            number = "$i"
        end
        namex = "AB1N"*number*"Fx"
        namey = "AB1N"*number*"Fy"
        nameux = "AB1N"*number*"Vx" 
        nameuy = "AB1N"*number*"Vy"
        namem = "AB1N"*number*"Mm"
        
        # namealpha = "AB1N"*number*"Alpha"

        fxmat[:,i] = outs[namex]
        fymat[:,i] = outs[namey]
        # Uxmat[:,i] = outs[nameux]
        # Uymat[:,i] = outs[nameuy]
        Mmat[:,i] = outs[namem]
        # alphamat[:,i] = outs[namealpha]
    end



    yamlfile = YAML.load_file("sn5_input.BD1.sum.yaml")
    



    nfea = 9
    rnodes = zeros(nfea)
    temp = yamlfile["Init_Nodes_E1"]

    for i = 1:nfea
        rnodes[i] = temp[i][3]
    end

    global dxmat = zeros(nt, nfea)
    global dymat = zeros(nt, nfea)
    global dzmat = zeros(nt, nfea)
    global dfx = zeros(nt, nfea)
    global dfy = zeros(nt, nfea)
    global dfz = zeros(nt, nfea)

    for i = 1:nfea
        if i<10
            number = "00$i"
        elseif i<100
            number = "0$i"
        else
            number = "$i"
        end

        namedx = "B1N"*number*"_TDxr"
        namedy = "B1N"*number*"_TDyr"
        namedz = "B1N"*number*"_TDzr"
        # namedfx = "B1N"*number*"_DFxR"
        # namedfy = "B1N"*number*"_DFyR"
        # namedfz = "B1N"*number*"_DFzR"

        dxmat[:,i] = outs[namedx]
        dymat[:,i] = outs[namedy]
        dzmat[:,i] = outs[namedz]
        # dfx[:,i] = outs[namedfx]
        # dfy[:,i] = outs[namedfy]
        # dfz[:,i] = outs[namedfz]
    end

    

    

    readflag = false
end

azimuth = outs["B1Azimuth"].*(pi/180)
tipdx = outs["B1TipTDxr"]
tipdy = outs["B1TipTDyr"]
tipdz = outs["B1TipTDzr"]


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
# blade = Rotors.Blade(rhub, rtip, rR, airfoils)
blade = Rotors.Blade(rvec, twistvec.*(pi/180), airfoils; rhub=rhub, rtip=rtip, precone)

turbine = true
rotor_r = Rotors.Rotor(Int(B), hubht, turbine; tilt, yaw)

if !@isdefined(defflag)
    defflag = true
end

if defflag
    aerostates, gxstates, mesh = initialize(blade, assembly, tvec; verbose=true)
    defflag = false
end

azimuth0 = 0.0
global system
system = initial_condition!(rotor_r, blade, assembly, env, aerostates, gxstates, mesh, tvec[1], azimuth0, pitch; verbose=true)

for i = 2:nt
    global system
    system = take_step!(aerostates, gxstates, mesh, rotor_r, blade, assembly, env, system, tvec, i, pitch; verbose=true, speakiter=100)
end

# aerostates, gxstates = simulate(rotor_r, blade, env, assembly, tvec; verbose=true)


#=
#Todo: #WorkLocation: 

1) system isn't getting passed correctly I don't think. Every other run errors because nlsolve can't handle what's going on. 
- OK.... now it's running every time. I'm still not convinced that system is getting passed properly. -> Well time to start comparing outputs/intermediate variables I guess. 
-> Why does Fx start negative????? -> I had double negatived the angle that the dynamic stall saw. 

2) Something is either getting updated too much, or isn't getting updated, or something. Because when the simulation does run, it's responding faster than it was before. 
=#

#Tip deflections
tipdef_x = [gxstates[i].points[end].u[1] for i in eachindex(tvec)]
tipdef_y = [gxstates[i].points[end].u[2] for i in eachindex(tvec)]
tipdef_z = [gxstates[i].points[end].u[3] for i in eachindex(tvec)]

tiptheta_x = zeros(nt)
tiptheta_y = zeros(nt)
tiptheta_z = zeros(nt)

tiptheta_xof = zeros(nt)
tiptheta_yof = zeros(nt)
tiptheta_zof = zeros(nt)

for i = 1:nt
    theta = Rotors.WMPtoangle(gxstates[i].points[end].theta)
    tiptheta_x[i] = theta[1]
    tiptheta_y[i] = theta[2]
    tiptheta_z[i] = theta[3]

    thetawmp = SVector(outs["B1TipRDxr"][i], outs["B1TipRDyr"][i], outs["B1TipRDzr"][i])
    theta = Rotors.WMPtoangle(thetawmp)
    tiptheta_xof[i] = theta[1]
    tiptheta_yof[i] = theta[2]
    tiptheta_zof[i] = theta[3]
end

# loads, cchistory, xds, gxstates, azimuth = Rotors.simulate(rotor_r, blade, env, assembly, tvec; verbose=true, speakiter=1000, g=inputfile["Gravity"], plotbool=false, plotiter=20)

using Plots, LaTeXStrings

tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)", legend=(0.9, 0.3))
plot!(tvec, fxmat[:,end], lab=L"$F_x$ - OF", seriescolor=:blue)
plot!(tvec, -fymat[:,end], lab=L"$F_y$ - OF", seriescolor=:red)
# plot!(tvec, Mmat[:,end], lab=L"$M_z$ - OF", seriescolor=:green)
plot!(tvec, aerostates.fx[:,end], lab=L"$F_x$ - R", linestyle=:dash)
plot!(tvec, aerostates.fy[:,end], lab=L"$F_y$ - R", linestyle=:dash)
# plot!(tvec, aerostates.M[:,end], lab=L"D_z", linestyle=:dash)
display(tiploads)

tipdefs = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)", legend=:topleft) #
plot!(tvec, tipdx, lab=L"$\delta x$ - OF", linestyle=:dash)
plot!(tvec, tipdy, lab=L"$\delta y$ - OF", linestyle=:dash)
plot!(tvec, tipdz, lab=L"$\delta z$ - OF", linestyle=:dash)
plot!(tvec, -tipdef_z, lab=L"\delta x")
plot!(tvec, tipdef_y, lab=L"\delta y")
plot!(tvec, tipdef_x, lab=L"\delta z")
display(tipdefs)


nothing