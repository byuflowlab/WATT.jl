using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, DynamicStallModels

DS = DynamicStallModels
of = OpenFASTsr

localpath = @__DIR__
cd(localpath)


### Read in OpenFAST files
ofpath = "./simpleNREL/" 
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
omega = edfile["RotSpeed"]*2*pi/60
shearexp = inflowwind["PLexp"]
tsr = omega*rtip/vinf

env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)



n = length(adblade["BlSpn"])
ne = Int(bdblade["station_total"])

if !@isdefined(readflag)
    readflag = true
end

if readflag
    # fullouts = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=6)
    fullouts = readdlm("./simpleNREL/sn5_input.out", skipstart=6)

    names = fullouts[1,:]

    # data = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=8)
    # data = readdlm("./simpleNREL/sn5_input.out", skipstart=8)
    data = Float64.(fullouts[3:end,:])

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]


    nt = length(tvec)


    global fxmat = zeros(nt, n)
    global fymat = zeros(nt, n)
    global Uxmat = zeros(nt, n)
    global Uymat = zeros(nt, n)
    global Mmat = zeros(nt, n)
    global alphamat = zeros(nt, n)

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
        
        namealpha = "AB1N"*number*"Alpha"

        fxmat[:,i] = outs[namex]
        fymat[:,i] = outs[namey]
        Uxmat[:,i] = outs[nameux]
        Uymat[:,i] = outs[nameuy]
        Mmat[:,i] = outs[namem]
        alphamat[:,i] = outs[namealpha]
    end

    global dxmat = zeros(nt, ne)
    global dymat = zeros(nt, ne)
    global dzmat = zeros(nt, ne)
    # global Dxmat = zeros(nt, ne)
    # global Dymat = zeros(nt, ne)
    # global Dzmat = zeros(nt, ne)
    global Vxmat = zeros(nt, ne)
    global Vymat = zeros(nt, ne)
    global Vzmat = zeros(nt, ne)
    global Wxmat = zeros(nt, ne)
    global Wymat = zeros(nt, ne)
    global Wzmat = zeros(nt, ne)

    for i = 1:ne
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
        # namebx = "B1N"*number*"_FxR"
        # nameby = "B1N"*number*"_FyR"
        # namebz = "B1N"*number*"_FzR"
        namevx = "B1N"*number*"_TVxr"
        namevy = "B1N"*number*"_TVyr"
        namevz = "B1N"*number*"_TVzr"
        namewx = "B1N"*number*"_RVxr"
        namewy = "B1N"*number*"_RVyr"
        namewz = "B1N"*number*"_RVzr"

        dxmat[:,i] = outs[namedx]
        dymat[:,i] = outs[namedy]
        dzmat[:,i] = outs[namedz]
        # Dxmat[:,i] = outs[namebx]
        # Dymat[:,i] = outs[nameby]
        # Dzmat[:,i] = outs[namebz]
        Vxmat[:,i] = outs[namevx]
        Vymat[:,i] = outs[namevy]
        Vzmat[:,i] = outs[namevz]
        Wxmat[:,i] = outs[namewx]
        Wymat[:,i] = outs[namewy]
        Wzmat[:,i] = outs[namewz]

    end
    readflag = false
end

azimuth = outs["B1Azimuth"].*(pi/180)


assembly = of.make_assembly(edfile, bdfile, bdblade)

### Prep the ASD rotor and operating conditions 
aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/Cylinder1.dat") 
aftypes[2] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/Cylinder2.dat") 
aftypes[3] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU40_A17.dat") 
aftypes[4] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU35_A17.dat") 
aftypes[5] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU30_A17.dat") 
aftypes[6] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU25_A17.dat") 
aftypes[7] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU21_A17.dat") 
aftypes[8] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/NACA64_A17.dat") 

# indices correspond to which airfoil is used at which station
af_idx = Int.(adblade["BlAFID"])


# create airfoil array
afs = aftypes[af_idx]

n = length(rvec)
airfoils = Vector{DS.Airfoil}(undef, n)
for i = 1:n
    airfoils[i] = make_dsairfoil(afs[i])
end

rR = rvec./rtip
blade = Rotors.Blade(rhub, rtip, rR, airfoils)


dsmodel = DS.BeddoesLeishman(DS.Indicial(), n, airfoils, 3)
dsmodelinit = Rotors.BeddoesLeishman()

if !@isdefined(runflag)
    runflag = true
end


if runflag
    tvec_r = range(0, 50, length=50001)
    loads, cchistory, xds, gxhistory, def_thetax = Rotors.simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade, env, assembly, tvec_r; verbose=true, dsmodel, dsmodelinit, speakiter=1000, g=inputfile["Gravity"], plotbool=false, plotiter=500, tipcorrection=nothing)
    runflag = false
end

    

#Tip deflections
tipdef_x = [gxhistory[i].points[end].u[1] for i in eachindex(tvec_r)]
tipdef_y = [gxhistory[i].points[end].u[2] for i in eachindex(tvec_r)]
tipdef_z = [gxhistory[i].points[end].u[3] for i in eachindex(tvec_r)]


Vx1 = [gxhistory[1].points[i].V[1] for i in eachindex(assembly.points)]
Vy1 = [gxhistory[1].points[i].V[2] for i in eachindex(assembly.points)]
Vz1 = [gxhistory[1].points[i].V[3] for i in eachindex(assembly.points)]


Uxhist = zeros(50001, n)
Uyhist = zeros(50001, n)
WR = zeros(50001, n)

for i in eachindex(tvec_r)
    Uxhist[i,:] = @. cchistory[i].W*sin(cchistory[i].phi)
    Uyhist[i,:] = @. cchistory[i].W*cos(cchistory[i].phi)
    WR[i,:] = @. sqrt(Uxhist[i,:].^2 + Uyhist[i,:].^2)
end

Wof = zeros(nt, n)
for i in eachindex(tvec)
    Wof[i,:] = @. sqrt(Uxmat[i,:].^2 + Uymat[i,:].^2)
end




using Plots, LaTeXStrings



# loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)")
# plot!(rvec, fxmat[end,:], lab=L"F_x")
# plot!(rvec, fymat[end,:], lab=L"F_y")
# # display(loadplt)



tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)", legend=:topleft)
plot!(tvec, fxmat[:,end], lab=L"$F_x$ - OF", seriescolor=:blue)
plot!(tvec, fymat[:,end], lab=L"$F_y$ - OF", seriescolor=:red)
# plot!(tvec, Mmat[:,end], lab=L"$M_z$ - OF", seriescolor=:green)
plot!(tvec_r, loads.Fx[:,end], lab=L"$F_x$ - R", linestyle=:dash)
plot!(tvec_r, loads.Fy[:,end], lab=L"$F_y$ - R", linestyle=:dash)
# plot!(tvec, loads.M[:,end], lab=L"D_z", linestyle=:dash)
display(tiploads)
# savefig("/Users/adamcardoza/Desktop/SimpleNRELTipLoads_10seconds_020823.png")


Uplt = plot(xaxis="Time (s)", yaxis="Tip Velocity (m/s)", legend=:outerright)
plot!(tvec, Uxmat[:,1], lab=L"$U_x$ - OF", seriescolor=:blue)
# plot!(tvec, Uymat[:,1], lab=L"$U_y$ - OF", seriescolor=:red)
plot!(tvec_r, Uxhist[:,1], lab=L"$U_x$ - Rotors", linestyle=:dash)
# plot!(tvec, Uyhist[:,1], lab=L"$U_y$ - Rotors", linestyle=:dash)
# display(Uplt)

#=
I think the windspeed values I'm plotting include the induced velocities. 

The wigglies in the velocities that I was getting from OpenFAST were due to the tower wiggling. 
=#


tipdefs2 = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)", legend=:outerright) #
plot!(tvec, dxmat[:,end], lab=L"\delta x - OF", linestyle=:dash)
plot!(tvec, dymat[:,end], lab=L"\delta y - OF", linestyle=:dash)
plot!(tvec, dzmat[:,end], lab=L"\delta z - OF", linestyle=:dash)
plot!(tvec_r, -tipdef_z, lab=L"\delta x - GX")
plot!(tvec_r, tipdef_y, lab=L"\delta y - GX")
plot!(tvec_r, tipdef_x, lab=L"\delta z - GX")
display(tipdefs2)
# savefig("/Users/adamcardoza/Desktop/SimpleNRELTipDeflections_10seconds_020823.png")

#=
2/9/23
The solutions start off quite similar, but they diverge over time. It appears that the OpenFAST solution is converging to a solution where the tip deflection is somewhere near 4.5 m, whereas Rotors quickly converges to a tip deflection just under 3 m. I don't want to be under predicting the loads. I don't know if I'm just neglecting an effect, or if potentially the stiffnesses aren't aligned (although I'm fairly confident that they are aligned). 

Additionally, my Y deflections are significantly larger than OpenFASTs'. I don't know what the source of that is. I don't know if they are the right sign either. 

I tried switching to a semi-implict approach (I had GXBeam deflect off of the i+1 aero state rather than the ith state.). I believe that currently I have OpenFAST running on an explicit approach. That didn't really make much of a difference, so I switched back. 

I'm trying switching the direction of the rotation (changing omega from negative to positive.). I don't know if that'll have any real effect on the final solution, but we'll see. -> That made the solution oscillate a lot, and eventually it looks like it'll go unstable. I need to plot more frequently to be sure that the thing goes unstable. -> Yeah, it definitely goes unstable.


Another thing that I could try is seeing if my loads application is off. Currently, I'm interpolating the aerodynamic loads, then integrating them across the GXBeam elements. Gasmi2013 or Sprague2014 (one of the two, I can't remember which) suggested that if you try to interpolate the loads, then you won't pass the entire load. 

Another thought that I had was what if there are forces that are getting applied to BeamDyn that I'm not applying to GXBeam. Right now I'm applying the aeroloads. Oh. I just realized I'm applying the moment to GXBeam, and not to OpenFAST. ... Well that didn't make a large difference. But let's revisit what loads are getting applied. Let's see. I'm applying the aero loads. I have gravity turned off currently for both. I have inertial loads applied because I'm applying a mass matrix. It's hard to say whether or not that is the issue because I'm not 100% confident that I'm interpolating correctly. I suppose I could interpolate for both so that I'm more confident that both are seeing the same mass matrix. That might make the solutions come closer, and might explain why one is converging to a different solution. 

Another approach is I could take a look at where the solutions start to converge and try and narrow down why they are diverging. I wonder if my DS states are the same. Because it looks like my loads are converging to a solution, but theirs kind of go off a little bit more before converging. Which... that could explain a major difference. ... Am I applying tip corrections? 


I should probably compare this against the static solution. 

=#




# airfoil = readdlm("./simpleturbine/Airfoils/NACA64_A17_coords.txt", skipstart=8)

# chordvec = adblade["BlChord"]

# sections = zeros(3, size(airfoil, 1), length(assembly.points))
# for ip = eachindex(assembly.points)
#     chord = chordvec[ip]
#     sections[1, :, ip] .= 0
#     sections[2, :, ip] .= chord .* (airfoil[:,1] .- 0.5)
#     sections[3, :, ip] .= chord .* airfoil[:,2]
# end

# mkpath("simpleturbine/viz/simpleturbine-simulation")
# write_vtk("simpleturbine/viz/simpleturbine-simulation/wind-turbine-blade-simulation", assembly, gxhistory[1:2500], tvec[1:2500]; sections = sections)

nothing