using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra

of = OpenFASTsr


### Read in OpenFAST files
ofpath = "./simpleturbine/" 
inputfile = of.read_inputfile("simple_input.fst", ofpath)
inflowwind = of.read_inflowwind("simple_inflowwind.dat", ofpath)
addriver = of.read_addriver("simple_ADdriver.dvr", ofpath)
adfile = of.read_adfile("simple_ADfile.dat", ofpath)
adblade = of.read_adblade("simple_ADblade.dat", ofpath)
edfile = of.read_edfile("simple_edfile.dat", ofpath)
bdfile = of.read_bdfile("simple_bdfile.dat", ofpath)
bdblade = of.read_bdblade("simple_bdblade.dat", ofpath)

rhub = addriver["HubRad(1)"]
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

# fullouts = readdlm("./simpleturbine/simple_ADdriver.1.out", skipstart=6)
fullouts = readdlm("./simpleturbine/simple_input.out", skipstart=6)

names = fullouts[1,:]

# data = readdlm("./simpleturbine/simple_ADdriver.1.out", skipstart=8)
data = readdlm("./simpleturbine/simple_input.out", skipstart=8)

outs = Dict(names[i] => data[:,i] for i in eachindex(names))

tvec = outs["Time"]


nt = length(tvec)

fxmat = zeros(nt, n)
fymat = zeros(nt, n)
Mmat = zeros(nt, n)
dxmat = zeros(nt, n)
dymat = zeros(nt, n)
dzmat = zeros(nt, n)

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
    namedx = "B1N"*number*"_TDxr"
    namedy = "B1N"*number*"_TDyr"
    namedz = "B1N"*number*"_TDzr"

    fxmat[:,i] = outs[namex]
    fymat[:,i] = outs[namey]
    Mmat[:,i] = outs[namem]
    dxmat[:,i] = outs[namedx]
    dymat[:,i] = outs[namedy]
    dzmat[:,i] = outs[namedz]
end

azimuth = outs["B1Azimuth"]


assembly = of.make_assembly(edfile, bdfile, bdblade)

E = 1.1e9 #6.83e10 #Young's Modulus
nu = 0.3 #Poisson's Ratio
m = 500 #2700.0 #kg/m^3 -> Decreasing the weight decreases the frequency of oscillation
h = 1.5 # Thickness (meters)
w = 1.5 # Width (meters)
mu = 0.01 #Damping ratio

G = E/(2*(1+nu))

A = h*w
Ay = A
Az = A
Iyy = w*h^3/12
Izz = w^3*h/12
J = Iyy + Izz

compliance = Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)])

newassembly = Rotors.update_assembly(assembly; compliance=compliance)

gxhistory = Rotors.simulate_gxbeam(rvec, rhub, rtip, tvec, azimuth, fxmat, fymat, Mmat, env, newassembly; verbose=true, speakiter=100, structural_damping=true, linear=false, g=inputfile["Gravity"])

    

#Tip deflections
tipdef_x = [gxhistory[i].points[end].u[1] for i in eachindex(tvec)]
tipdef_y = [gxhistory[i].points[end].u[2] for i in eachindex(tvec)]
tipdef_z = [gxhistory[i].points[end].u[3] for i in eachindex(tvec)]













using Plots, LaTeXStrings



loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)")
plot!(rvec, fxmat[1,:], lab=L"F_x")
plot!(rvec, fymat[1,:], lab=L"F_y")
# display(loadplt)


tiploads = plot(xaxis="Time (s)", yaxis="Tip Load (N)")
plot!(tvec, fxmat[:,end], lab=L"F_x")
plot!(tvec, fymat[:,end], lab=L"F_y")
plot!(tvec, Mmat[:,end], lab=L"M")
# display(tiploads)

tipdefs = plot(xaxis="Time (s)", yaxis="Tip Deflection (m)", ylims=(-0.1, 0.1))
plot!(tvec, dxmat[:,end], lab=L"\delta x - OF", linestyle=:dash)
plot!(tvec, dymat[:,end], lab=L"\delta y - OF", linestyle=:dash)
plot!(tvec, dzmat[:,end], lab=L"\delta z - OF", linestyle=:dash)
plot!(tvec, -tipdef_x, lab=L"\delta x - GX")
plot!(tvec, tipdef_y, lab=L"\delta y - GX")
plot!(tvec, tipdef_z, lab=L"\delta z - GX")
display(tipdefs)

nothing