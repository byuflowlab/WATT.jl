using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra

of = OpenFASTsr


### Read in OpenFAST files
ofpath = "./OpenFAST_NREL5MW/" 
inputfile = of.read_inputfile("NREL5MW_input.fst", ofpath)
inflowwind = of.read_inflowwind("NREL5MW_inflowwind.dat", ofpath)
# addriver = of.read_addriver("NREL5MW_ADdriver.inp", ofpath)
adfile = of.read_adfile("NREL5MW_ADfile.dat", ofpath)
adblade = of.read_adblade("NREL5MW_ADblade.dat", ofpath)
edfile = of.read_edfile("NREL5MW_edfile.dat", ofpath)
bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
bdblade = of.read_bdblade("NREL5MW_bdblade.dat", ofpath)

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

if readflag
    # fullouts = readdlm("./OpenFAST_NREL5MW/NREL5MW_ADdriver.1.out", skipstart=6)
    fullouts = readdlm("./OpenFAST_NREL5MW/NREL5MW_input.out", skipstart=6)

    names = fullouts[1,:]

    # data = readdlm("./OpenFAST_NREL5MW/NREL5MW_ADdriver.1.out", skipstart=8)
    data = readdlm("./OpenFAST_NREL5MW/NREL5MW_input.out", skipstart=8)

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]


    nt = length(tvec)


    global fxmat = zeros(nt, n)
    global fymat = zeros(nt, n)
    # global Mmat = zeros(nt, n)
    global dxmat = zeros(nt, n)
    global dymat = zeros(nt, n)
    global dzmat = zeros(nt, n)
    # global Dxmat = zeros(nt, n)
    # global Dymat = zeros(nt, n)
    # global Dzmat = zeros(nt, n)

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
        namedx = "B1N"*number*"_TDxr"
        namedy = "B1N"*number*"_TDyr"
        namedz = "B1N"*number*"_TDzr"
        # namebx = "B1N"*number*"_FxR"
        # nameby = "B1N"*number*"_FyR"
        # namebz = "B1N"*number*"_FzR"

        fxmat[:,i] = outs[namex]
        fymat[:,i] = outs[namey]
        # Mmat[:,i] = outs[namem]
        dxmat[:,i] = outs[namedx]
        dymat[:,i] = outs[namedy]
        dzmat[:,i] = outs[namedz]
        # Dxmat[:,i] = outs[namebx]
        # Dymat[:,i] = outs[nameby]
        # Dzmat[:,i] = outs[namebz]
    end
    readflag = false
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