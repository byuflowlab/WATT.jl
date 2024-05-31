#=
Testing the simulate! function that only uses the aerodynamic models.

Adam Cardoza4/24/24
=#

using Revise
using OpenFASTTools, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, DynamicStallModels
# using Infiltrator
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

if !@isdefined(readflag)
    readflag = true
end

if readflag
    println("Reading OpenFAST files...")
    # fullouts = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=6)
    fullouts = readdlm("./sn5_ADdriver.out", skipstart=6)

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


    readflag = false
end



# assembly = of.make_assembly(edfile, bdfile, bdblade)

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

# tvec_r = tvec[1]:0.01:tvec[end]
# tvec_r = tvec[1:5]
tvec_r = tvec

defflag = true
if defflag
    println("initialize...")
    aerostates, mesh = Rotors.initialize(blade, tvec_r; verbose=true)
    defflag = false
end

if !@isdefined(runflag)
    runflag = true
end


runflag = true
if runflag 
    println("Running simulation...")
    Rotors.simulate!(aerostates, mesh, rotor_r, blade, env, tvec_r; verbose=true)
    

    runflag = false 
end





# Vx1 = [gxhistory[1].points[i].V[1] for i in eachindex(assembly.points)]
# Vy1 = [gxhistory[1].points[i].V[2] for i in eachindex(assembly.points)]
# Vz1 = [gxhistory[1].points[i].V[3] for i in eachindex(assembly.points)]


# Uxhist = zeros(50001, n)
# Uyhist = zeros(50001, n)
# WR = zeros(50001, n)

# for i in eachindex(tvec)
#     Uxhist[i,:] = @. cchistory[i].W*sin(cchistory[i].phi)
#     Uyhist[i,:] = @. cchistory[i].W*cos(cchistory[i].phi)
#     WR[i,:] = @. sqrt(Uxhist[i,:].^2 + Uyhist[i,:].^2)
# end

# Wof = zeros(nt, n)
# for i in eachindex(tvec)
#     Wof[i,:] = @. sqrt(Uxmat[i,:].^2 + Uymat[i,:].^2)
# end




using Plots, LaTeXStrings

println("Plotting... ")

# loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)", leg=:topleft)
# plot!(rvec, fxmat[end,:], lab=L"F_x")
# plot!(rvec, fymat[end,:], lab=L"F_y")
# display(loadplt)

# twistplt = plot(xaxis="Radius (m)", yaxis="Weiner-Milenkovic Parameter", leg=:topleft)
# plot!(rvec, rdx[end,:], lab=L"\theta_x")
# plot!(rvec, rdy[end,:], lab=L"\theta_y")
# plot!(rvec, rdz[end,:], lab=L"\theta_z")
# display(twistplt)

#The distributed load the BeamDyn seens. 
# loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)", leg=:topleft)
# plot!(rvec, dfx[end,:], lab=L"F_x")
# plot!(rvec, dfy[end,:], lab=L"F_y")
# plot!(rvec, dfy[end,:], lab=L"F_y")
# display(loadplt)

idxs = 1:length(tvec)


xtiploads = plot(yaxis=L"
$F_x$ Tip Load (N/m)")
plot!(xtiploads, tvec[idxs], fxmat[idxs,end], lab="OF", seriescolor=1)
plot!(xtiploads, tvec_r, aerostates.Fx[:,end], lab="R", linestyle=:dash, seriescolor=2, lw=2.5)

ytiploads = plot(tvec[idxs], -fymat[idxs,end], lab="OF", seriescolor=1, xaxis="Time (s)", yaxis=L"
$F_y$ Tip Load (N/m)")
plot!(ytiploads, tvec_r, aerostates.Fy[:,end], lab="R", linestyle=:dash, seriescolor=2, lw=2.5)

# plot!(tiploads, tvec, loads.M[:,end], lab=L"D_z", linestyle=:dash)
# plot!(tiploads, tvec, Mmat[:,end], lab=L"$M_z$ - OF", seriescolor=:green)

tiploadsplt = plot(xtiploads, ytiploads, layout=(2,1))
display(tiploadsplt)
# savefig(tiploadsplt, "/Users/adamcardoza/Desktop/SimpleNRELTipLoads_varyingairfoils_chords_twists_gravity_shear_1seconds_042624_aero_only.png")






# nodeidx = 300

# # alphaplt = plot(xaxis="Time (s)", yaxis="Angle of Attack (deg)")
# # plot!(tvec, aerostates.alpha[:,nodeidx].*(180/pi), lab="R")
# # plot!(tvec, outs["AB1N$nodeidx"*"Alpha"], lab="OF")
# # plot!(tvec[2:end], iouts["alpha"].*(180/pi), lab="AD")
# # display(alphaplt)


# Wof = @. sqrt(outs["AB1N$nodeidx"*"Vx"]^2 + outs["AB1N$nodeidx"*"Vy"]^2)


# # Uplt = plot(xaxis = "Time (s)", yaxis="Inflow Velocity (m/s)")
# # plot!(tvec, aerostates.W[:,nodeidx], lab="R")
# # plot!(tvec, Wof, lab="OF")
# # plot!(tvec[2:end], iouts["U"], lab="AD")
# # display(Uplt)

# # Vxr = cchistory[:,nodeidx].u./cchistory[:,nodeidx].G

# # Vxplt = plot(xaxis = "Time (s)", yaxis="Vx (m/s)", leg=:bottomright)
# # plot!(Vxmat[:,nodeidx], lab="R")
# # plot!(outs["AB1N$nodeidx"*"Vx"], lab="OF")
# # # plot!(bouts["Vx"], lab="AD")
# # # display(Vxplt)

# # Vyplt = plot(xaxis = "Time (s)", yaxis="Vy (m/s)", leg=:topright)
# # plot!(Vymat[:,nodeidx], lab="R")
# # plot!(outs["AB1N$nodeidx"*"Vy"], lab="OF")
# # # plot!(bouts["Vy"], lab="AD")
# # # display(Vyplt)

# # aplt = plot(xaxis="Time (s)", yaxis="Axial induction factor")
# # plot!(tvec, aerostates.a[:,nodeidx], lab="R")
# # plot!(tvec, outs["AB1N$nodeidx"*"AxInd"], lab="AD")
# # # display(aplt)

# # applt = plot(xaxis="Time (s)", yaxis="Tangential Induction factor")
# # plot!(tvec, aerostates.ap[:,nodeidx], lab="R")
# # plot!(tvec, outs["AB1N$nodeidx"*"TnInd"], lab="AD")
# # # display(applt)

# phiplt = plot(xaxis="Time (s)", yaxis="Inflow Angle (deg)", leg=:topleft)
# plot!(tvec, aerostates.phi[:,nodeidx].*(180/pi), lab="R")
# plot!(tvec, outs["AB1N$nodeidx"*"Phi"], lab="AD")
# display(phiplt)
# # savefig(phiplt, "simpleNREL_inflow_angle_1sec_7.12.23.pdf")

# # twist = @. (aerostates.phi[:,nodeidx] + -aerostates.alpha[:,nodeidx])*180/pi
# # twistof = @. outs["AB1N$nodeidx"*"Phi"] + -outs["AB1N$nodeidx"*"Alpha"]

# # twistplt = plot(xaxis="Time (s)", yaxis="Twist Angle (deg)")
# # plot!(tvec, twist, lab="R")
# # plot!(tvec, outs["AB1N$nodeidx"*"Theta"], lab="AD")
# # # plot!(tvec, twistof, lab="OF", linestyle=:dash)
# # display(twistplt)

# thetaplt = plot(xaxis = "Time (s)", yaxis=L"\delta_\theta (deg)", leg=:topleft)
# plot!(tvec, -tiptheta_z.*(180/pi), lab="R - x")
# plot!(tvec, tiptheta_y.*(180/pi), lab="R - y")
# plot!(tvec, tiptheta_x.*(180/pi), lab="R - z")
# # plot!(tvec, outs["B1TipRDxr"], lab="OF - x", linestyle=:dash)
# # plot!(tvec, outs["B1TipRDyr"], lab="OF - y", linestyle=:dash)
# # plot!(tvec, outs["B1TipRDzr"], lab="OF - z", linestyle=:dash)
# plot!(tvec, tiptheta_xof.*(180/pi), lab="OF - x", linestyle=:dash)
# plot!(tvec, tiptheta_yof.*(180/pi), lab="OF - y", linestyle=:dash)
# plot!(tvec, tiptheta_zof.*(180/pi), lab="OF - z", linestyle=:dash)
# display(thetaplt)
# # savefig(thetaplt, "simpleNREL_angular_def_1s_7.12.23.pdf")

# # aziplt = plot(tvec, azimuth, lab="AeroDyn", xaxis="Time (s)", yaxis="Azimuthal angle (deg)", leg=:top)
# # plot!(tvec, aerostates.azimuth.*(180/pi), lab="Rotors")
# # plot!(tvec, outs["Azimuth"], lab="ElastoDyn", linestyle=:dash)
# # display(aziplt)
# # savefig(aziplt, "/Users/adamcardoza/Desktop/azimuthplot_AD_ED.png") #The Aerodyn azimuth is a value for wake correction models. 

# # Fxamp_of = maximum(fxmat[2500:end,end]) - minimum(fxmat[2500:end,end]) #427.0
# # Fxamp_R = maximum(loads.Fx[2500:end, end])-minimum(loads.Fx[2500:end, end]) #448.517

# Mx_of = zeros(nt)
# Mx_r = zeros(nt)


# for i = 1:nt
#     Mx_of[i] = of.root_bending_moment(rvec, fxmat[i,:])
#     Mx_r[i] = of.root_bending_moment(rvec, aerostates.fx[i,:])
# end

# using Dierckx
# Mx_r_smooth_fit = Spline1D(tvec, Mx_r;w=ones(length(tvec)), k=3, bc="nearest", s=1.0e14)

# Mx_r_smooth = Mx_r_smooth_fit.(tvec)

# Mplt = plot(xaxis="Time (s)", yaxis=L"Root Bending Moment $(N\cdot m)$", leg=:topright)
# plot!(tvec, Mx_of, lab="OF")
# plot!(tvec, Mx_r, lab="R")
# # plot!(tvec, Mx_r_smooth, lab="smooth")
# display(Mplt)
# # savefig(Mplt, "/Users/adamcardoza/Desktop/SimpleNRELrootbendingmoment_varyingairfoils_chords_twists_gravity_shear_1seconds_071223.pdf")

# # Mxerr = @. 100*(Mx_r-Mx_of)/Mx_of

# # avgMerr = mean(Mxerr)

# # Mxamp_of = maximum(Mx_of[2500:end])-minimum(Mx_of[2500:end])
# # Mxamp_r = maximum(Mx_r[2500:end])-minimum(Mx_r[2500:end])



# # DEMx_of = of.damage_equivalent_load(Mx_of)
# # DEMx_r = of.damage_equivalent_load(Mx_r)
# # # DEMx_r = of.damage_equivalent_load(Mx_r_smooth)


# # DEMx_err = 100*(DEMx_r-DEMx_of)/DEMx_of

# #=
# I'm a little confused. I thought with the damage equivalent moment, if I smoothed the sucker out that it would become.... less. I guess not? Is it purely dependent on the the amplitude of the loads? 

# Also, with a difference in the max amplitudes of 100,000 I wouldn't expect the difference in the DEM to be 420,000.... 

# =#



# ### Comparing blade mass
# # Mass calculated by BeamDyn: 16844.336
# m_bd = 16844.336
# m_r = Rotors.get_blade_weight(assembly)

# masserr = 100*(m_r-m_bd)/m_bd

# @show masserr



### Computing gravity
# g = inputfile["Gravity"]
# gravvec = [SVector(-g*cos(aerostates.azimuth[i]), -g*sin(aerostates.azimuth[i]), 0.0) for i = 1:nt] #In the rotating hub reference frame. 

# gidx = nt

# rx = rtip*cos(aerostates.azimuth[gidx])
# ry = rtip*sin(aerostates.azimuth[gidx])

# graviplt = plot([0, -gravvec[gidx][2]], [0, gravvec[gidx][1]], leg=false, xaxis="-Y", yaxis="X")
# plot!([0, ry], [0, rx], linecolor=:black, lab="Blade position")
# display(graviplt)

# grot = Rotors.rotate_z(gravvec[gidx]..., aerostates.azimuth[gidx]; T=true)
#=
Wait... everything goes into the negative x direction... which means I have rx and ry wrong. 

Well... From what I can tell, I implemented gravity right. It's where I expect it to be... and at a value that I expect. So.... That's odd. I wonder if I'm spinning it the wrong way? Would that make a difference? I suppose it could push my deflections out of phase. 
=#


# anim = @animate for i in eachindex(tvec)
#     plot(xaxis="Radius (m)", yaxis="Load (N/m)", leg=:topleft)
#     plot!(rvec, cchistory[i, :].Np, lab="N")
#     plot!(rvec, cchistory[i,:].Tp, lab="T")
# end every 10
# gif(anim, "ccloads.gif", fps = 15)




# anim = @animate for i in eachindex(tvec)
#     plt1 = plot(leg=:topleft, yaxis="Distributed Load (N/m)", ylims=(-750, 8000))
#     plot!(rvec, fxmat[i,:], lab=L"$F_x$ - OF", linestyle=:dash)
#     plot!(rvec, -fymat[i,:], lab=L"$F_y$ - OF", linestyle=:dash)
#     plot!(rvec, aerostates.fx[i,:], lab=L"$F_x$ - R")
#     plot!(rvec, aerostates.fy[i,:], lab=L"$F_y$ - R")

#     plt2 = plot(xaxis="Radius (m)", yaxis="Deflection (m)", legend=:topleft, ylims=(-1.0, 5.1)) #
#     plot!(rnodes, dxmat[i,:], lab=L"x - OF", linestyle=:dash)
#     plot!(rnodes, dymat[i,:], lab=L"y - OF", linestyle=:dash)
#     plot!(rnodes, dzmat[i,:], lab=L"z - OF", linestyle=:dash)
#     plot!(rvec, -defz_gx[i,:], lab=L"x - GX")
#     plot!(rvec, defy_gx[i,:], lab=L"y - GX")
#     plot!(rvec, defx_gx[i,:], lab=L"z - GX")

#     plot(plt1, plt2, layout=(2,1))
# end #every 3
# gif(anim, "loads_defs.gif", fps = 15)

#=
When you look at things in the grand scheme of things, it really isn't that big of a difference. Now the overall loading might have some significant differences, but I don't know if I can expect to be any closer. Like... I'm using a different convergence tolerance I think (at least on the BEMT). Maybe implicitly I am. I guess I could check that. And check that the tip correction isn't getting applied. 

-[] Check tip and hub corrections (that they match)
-[] Check that the BEMT tolerances are the same. 
6/22/23
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