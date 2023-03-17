using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra

#=
Applying the aerodynamic loads from OpenFAST to GXBeam. 

Adding in the comparison using the decoupled results from BeamDyn

=#

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

rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]
rfrac = bdblade["rfrac"]

rho = inputfile["AirDens"]
mu = inputfile["KinVisc"]
vinf = inflowwind["HWindSpeed"]
a = inputfile["SpdSound"]
omega = edfile["RotSpeed"]*2*pi/60
shearexp = inflowwind["PLexp"]
tsr = omega*rtip/vinf


env = Rotors.environment(rho, mu, a, vinf, 0.0, shearexp)



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
    data = readdlm("./simpleNREL/sn5_input.out", skipstart=8)

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]


    nt = length(tvec)


    global fxmat = zeros(nt, n)
    global fymat = zeros(nt, n)
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
        namem = "AB1N"*number*"Mm"
        
        # namealpha = "AB1N"*number*"Alpha"

        fxmat[:,i] = outs[namex]
        fymat[:,i] = outs[namey]
        # Mmat[:,i] = outs[namem]
        # alphamat[:,i] = outs[namealpha]
    end

    # global dxmat = zeros(nt, ne)
    # global dymat = zeros(nt, ne)
    # global dzmat = zeros(nt, ne)
    # global Dxmat = zeros(nt, ne)
    # global Dymat = zeros(nt, ne)
    # global Dzmat = zeros(nt, ne)
    # global Vxmat = zeros(nt, ne)
    # global Vymat = zeros(nt, ne)
    # global Vzmat = zeros(nt, ne)
    # global Wxmat = zeros(nt, ne)
    # global Wymat = zeros(nt, ne)
    # global Wzmat = zeros(nt, ne)

    # for i = 1:ne
    #     if i<10
    #         number = "00$i"
    #     elseif i<100
    #         number = "0$i"
    #     else
    #         number = "$i"
    #     end

    #     namedx = "B1N"*number*"_TDxr"
    #     namedy = "B1N"*number*"_TDyr"
    #     namedz = "B1N"*number*"_TDzr"
    #     # namebx = "B1N"*number*"_FxR"
    #     # nameby = "B1N"*number*"_FyR"
    #     # namebz = "B1N"*number*"_FzR"
    #     namevx = "B1N"*number*"_TVxr"
    #     namevy = "B1N"*number*"_TVyr"
    #     namevz = "B1N"*number*"_TVzr"
    #     namewx = "B1N"*number*"_RVxr"
    #     namewy = "B1N"*number*"_RVyr"
    #     namewz = "B1N"*number*"_RVzr"

    #     dxmat[:,i] = outs[namedx]
    #     dymat[:,i] = outs[namedy]
    #     dzmat[:,i] = outs[namedz]
    #     # Dxmat[:,i] = outs[namebx]
    #     # Dymat[:,i] = outs[nameby]
    #     # Dzmat[:,i] = outs[namebz]
    #     Vxmat[:,i] = outs[namevx]
    #     Vymat[:,i] = outs[namevy]
    #     Vzmat[:,i] = outs[namevz]
    #     Wxmat[:,i] = outs[namewx]
    #     Wymat[:,i] = outs[namewy]
    #     Wzmat[:,i] = outs[namewz]

    # end

    bdread = readdlm("./simpleNREL/sn5_bddriver.out", skipstart=6)

    bdnames = bdread[1,:]

    # data = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=8)
    bddata = Float64.(bdread[3:end,:])

    bdouts = Dict(bdnames[i] => bddata[:,i] for i in eachindex(bdnames))

    tvecb = bdouts["Time"]


    ntb = length(tvecb)

    global dxmat = zeros(ntb, ne)
    global dymat = zeros(ntb, ne)
    global dzmat = zeros(ntb, ne)

    for i = 1:ne
        if i<10
            number = "00$i"
        elseif i<100
            number = "0$i"
        else
            number = "$i"
        end

        namedx = "N"*number*"_TDxr"
        namedy = "N"*number*"_TDyr"
        namedz = "N"*number*"_TDzr"

        dxmat[:,i] = bdouts[namedx]
        dymat[:,i] = bdouts[namedy]
        dzmat[:,i] = bdouts[namedz]

    end
    

    readflag = false
end




azimuth = outs["B1Azimuth"].*(pi/180)


assembly = of.make_assembly(edfile, bdfile, bdblade)

# E = 1.1e9 #6.83e10 #Young's Modulus
# nu = 0.3 #Poisson's Ratio
# m = 500 #2700.0 #kg/m^3 -> Decreasing the weight decreases the frequency of oscillation
# h = 1.5 # Thickness (meters)
# w = 1.5 # Width (meters)
# mu = 0.01 #Damping ratio

# G = E/(2*(1+nu))

# A = h*w
# Ay = A
# Az = A
# Iyy = w*h^3/12
# Izz = w^3*h/12
# J = Iyy + Izz

# compliance = Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)])

# newassembly = Rotors.update_assembly(assembly; compliance=compliance)

# Mmat = zero(fymat)
if !@isdefined(runflag)
    runflag = true
end



# if runflag
    # gxhistory = Rotors.simulate_gxbeam(rvec, rhub, rtip, tvec, azimuth, fxmat, fymat, Mmat, env, assembly; verbose=true, speakiter=100, structural_damping=true, linear=false, g=inputfile["Gravity"])
    
#     runflag = false
# end

nn,mm = size(fxmat) #Note: This is me just setting the loads to a constant distributed load. 
fxmat = zeros(nn, mm)
fymat = zeros(nn, mm)
Mmat = zeros(nn, mm)

# fxmat[1,:] .= 1000.0


gxstate = Rotors.steady_simulate_gxbeam(rvec, azimuth, fxmat, fymat, Mmat, env, assembly; verbose=true, speakiter=100, structural_damping=true, linear=false, g=inputfile["Gravity"])

    

#Tip deflections
# tipdef_x = [gxhistory[i].points[end].u[1] for i in eachindex(tvec)]
# tipdef_y = [gxhistory[i].points[end].u[2] for i in eachindex(tvec)]
# tipdef_z = [gxhistory[i].points[end].u[3] for i in eachindex(tvec)]


def_x = [gxstate.points[i].u[1] for i in eachindex(gxstate.points)]
def_y = [gxstate.points[i].u[2] for i in eachindex(gxstate.points)]
def_z = [gxstate.points[i].u[3] for i in eachindex(gxstate.points)]


# Vx1 = [gxhistory[1].points[i].V[1] for i in eachindex(assembly.points)]
# Vy1 = [gxhistory[1].points[i].V[2] for i in eachindex(assembly.points)]
# Vz1 = [gxhistory[1].points[i].V[3] for i in eachindex(assembly.points)]









using Plots, LaTeXStrings



# loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)")
# plot!(rvec, fxmat[end,:], lab=L"F_x")
# plot!(rvec, fymat[end,:], lab=L"F_y")
# # display(loadplt)



# tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)", legend=:topleft)
# plot!(tvec, fxmat[:,end], lab=L"F_x", seriescolor=:blue)
# plot!(tvec, fymat[:,end], lab=L"F_y", seriescolor=:red)
# plot!(tvec, Mmat[:,end], lab=L"M_z", seriescolor=:green)
# plot!(tvec, Dxmat[:,end], lab=L"D_x", linestyle=:dash)
# plot!(tvec, Dymat[:,end], lab=L"D_y", linestyle=:dash)
# plot!(tvec, Dzmat[:,end], lab=L"D_z", linestyle=:dash)
# display(tiploads)
# savefig("/Users/adamcardoza/Desktop/SimpleNRELTipLoads_10seconds.png")


# tipdefs = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)") #
# plot!(tvec, dxmat[:,end], lab=L"\delta x - OF", linestyle=:dash, seriescolor=:blue)
# plot!(tvec, dymat[:,end], lab=L"\delta y - OF", linestyle=:dash, seriescolor=:red)
# plot!(tvec, dzmat[:,end], lab=L"\delta z - OF", linestyle=:dash, seriescolor=:green)
# # display(tipdefs)
# # savefig("/Users/adamcardoza/Desktop/SimpleTurbineTipDeflections.png")


# tipdefs2 = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)", legend=:outerright) #
# plot!(tvecb, tipdx, lab=L"\delta x - OF")
# plot!(tvecb, tipdy, lab=L"\delta y - OF")
# plot!(tvecb, tipdz, lab=L"\delta z - OF")
# plot!(tvec, -tipdef_z, lab=L"\delta x - GX", linestyle=:dash)
# plot!(tvec, tipdef_y, lab=L"\delta y - GX", linestyle=:dash)
# plot!(tvec, tipdef_x, lab=L"\delta z - GX", linestyle=:dash)
# display(tipdefs2)
# savefig("/Users/adamcardoza/Desktop/SimpleNRELTipDeflections_10seconds.png")

#=
- 3/7/23 So at this point I've gotten the solution running without the RPM for statically or not statically initialized. Once I add the RPM in, then the only solution I can get running is a not initialized solution. My solution is nowhere near BeamDyn's... but without the rpm... it's pretty close. -> That makes me wonder if I'm applying the angular velocity incorrectly to BeamDyn or GXBeam. -> So I realized that if my initial condition isn't matching, then there ain't no way the rest of the simulation is going to match. So I'm going to try and get a steady state simulation to match. 
As I was just preparing the BeamDyn Driver input file for simulation I just realized that I didn't change the dt, so who knows what it was doing on the half steps during those time domain sims. No wonder the simulation was so different.... And took so long. 


- 3/6/23 So I've gotten pretty close to matching GXBeam and BeamDyn responding to loads from OpenFAST. (It took a while to learn the Fortran, and get things automated). I had things working earlier today, but now I've broken them again. I returned the things I had changed back to their original state... but alas, it be broken again. :| -> It's very possible that I wasn't writing all of the loadings to file for BeamDyn to read... so only a portion of the blade was seeing loads. I don't know how Fortran wasn't having a cow. 

- 3/4/23 A couple of weeks or so ago Dr. Ning wasn't super jazzed that I moved on from this point and wanted me to revisit the problem comparing a decoupled BeamDyn to the decoupled GXBeam. That'll help me see if it's the same problem. 


- 2/8/23 well, now it's matching much more closely. The problem was that my initial conditions were zero, instead of getting the initial velocity, but that was a problem inside of GXBeam. Additionally, I wasn't pre-integrating my loads... so it was getting wonky loads. Taylor seems to think that this is probably the best I'm going to get out of applying the OpenFAST aero loads to GXBeam, and that I should move onto "comparing apples to apples", i.e. applying my own aero loads to GXBeam. I think I'm going to run a longer simulation first and see if my deflections converge to the same spot. It is also interesting that my Y deflections are higher in magnitude and are positive. I mean... they could be positive because I'm rotating the dumb thing the wrong way... but we'll see. Also... they could be higher in magnitude because they lack the coupled aerodynamic damping (like what we see on the X deflection.)

Well, in a 50 second simulation, it doesn't converge... which is a bummer. I mean, I have no clue if it would get to a point where it would converge... or if it would just be off. I have no expectation that this should converge (well, now I don't after having talked with Taylor.) The Y deflections appear to be converging, which is good. 

Another interesting thing is that the OpenFAST loads don't seem to have converged, but the deflections have... which seems exceedingly odd. 

I also haven't looked at things like gravity yet. Wait... does GXBeam have gravity while OpenFAST doesn't? ... No, niether has gravity applied. Yeah... I think it might be good to try a simple combined case. 
=#



steadyplt = plot(xaxis="Radius (m)", yaxis="Deflection (m)", leg=:topleft)
plot!(rvec, -def_z, label="X GX")
plot!(rvec, def_y, label="Y GX")
plot!(rvec, def_x, label="Z GX")
plot!(rvec, dxmat[end,:], label="X OF", linestyle=:dash)
plot!(rvec, dymat[end,:], label="Y OF", linestyle=:dash)
plot!(rvec, dzmat[end,:], label="Z OF", linestyle=:dash)
# hline!([tipdx[end]], lab=L"\delta x - OF")
# hline!([tipdy[end]], lab=L"\delta y - OF")
# hline!([tipdz[end]], lab=L"\delta z - OF")
display(steadyplt)

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