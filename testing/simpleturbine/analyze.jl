using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra

#=
Analyze the results of the openFAST simulation. 

Adam Cardoza
=#

of = OpenFASTsr

localpath = @__DIR__
cd(localpath)


### Read in OpenFAST files
ofpath = "./" 
inputfile = of.read_inputfile("simple_input.fst", ofpath)
inflowwind = of.read_inflowwind("simple_inflowwind.dat", ofpath)
# addriver = of.read_addriver("simple_ADdriver.inp", ofpath)
adfile = of.read_adfile("simple_ADfile.dat", ofpath)
adblade = of.read_adblade("simple_ADblade.dat", ofpath)
edfile = of.read_edfile("simple_edfile.dat", ofpath)
bdfile = of.read_bdfile("simple_bdfile.dat", ofpath)
bdblade = of.read_bdblade("simple_bdblade.dat", ofpath)

rhub = edfile["HubRad"]
rvec = adblade["BlSpn"] .+ rhub
rtip = rvec[end]

rho = inputfile["AirDens"]
mu = inputfile["KinVisc"]
vinf = inflowwind["HWindSpeed"]
a = inputfile["SpdSound"]
omega = edfile["RotSpeed"]
shearexp = inflowwind["PLexp"]

env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)



n = length(adblade["BlSpn"])

if !@isdefined(readflag)
    readflag = true
end

readflag = true

if readflag
    # fullouts = readdlm("./simple_ADdriver.1.out", skipstart=6)
    fullouts = readdlm("./simple_input.out", skipstart=6)

    names = fullouts[1,:]

    # data = readdlm("./simple_ADdriver.1.out", skipstart=8)
    # data = readdlm("./simple_input.out", skipstart=8) #Todo: There's probably a faster way to do this.
    data = Float64.(fullouts[3:end,:]) 

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]


    nt = length(tvec)


    global fxmat = zeros(nt, n)
    global fymat = zeros(nt, n)
    # global Mmat = zeros(nt, n)

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
        # namem = "AB1N"*number*"Mm"
        

        fxmat[:,i] = outs[namex]
        fymat[:,i] = outs[namey]
        # Mmat[:,i] = outs[namem]
        
    end


    ne = Int(bdblade["station_total"])
    global dxmat = zeros(nt, ne) 
    global dymat = zeros(nt, ne)
    global dzmat = zeros(nt, ne)
    # global Dxmat = zeros(nt, ne)
    # global Dymat = zeros(nt, ne)
    # global Dzmat = zeros(nt, ne)
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

        dxmat[:,i] = outs[namedx] #Todo: I don't think these have been changing with my changes.... 
        dymat[:,i] = outs[namedy]
        dzmat[:,i] = outs[namedz]
        # Dxmat[:,i] = outs[namebx]
        # Dymat[:,i] = outs[nameby]
        # Dzmat[:,i] = outs[namebz]
    end
    # readflag = false
end

azimuth = outs["B1Azimuth"].*(pi/180)














using Plots, LaTeXStrings



loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)")
plot!(rvec, fxmat[end,:], lab=L"F_x")
plot!(rvec, fymat[end,:], lab=L"F_y")
# display(loadplt)



tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)")
plot!(tvec, fxmat[:,end], lab=L"F_x", seriescolor=:blue)
plot!(tvec, fymat[:,end], lab=L"F_y", seriescolor=:red)
# plot!(tvec, Mmat[:,end], lab=L"M_z", seriescolor=:green)
# plot!(tvec, Dxmat[:,end], lab=L"D_x", linestyle=:dash)
# plot!(tvec, Dymat[:,end], lab=L"D_y", linestyle=:dash)
# plot!(tvec, Dzmat[:,end], lab=L"D_z", linestyle=:dash)
display(tiploads)
# savefig("/Users/adamcardoza/Desktop/NREL5MWTipLoads.png")


tipdefs = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)") #
plot!(tvec, dxmat[:,end], lab=L"\delta x - OF", linestyle=:dash, seriescolor=:blue)
plot!(tvec, dymat[:,end], lab=L"\delta y - OF", linestyle=:dash, seriescolor=:red)
plot!(tvec, dzmat[:,end], lab=L"\delta z - OF", linestyle=:dash, seriescolor=:green)
display(tipdefs)
# savefig("/Users/adamcardoza/Desktop/NREL5MWTipDeflections.png")



nothing