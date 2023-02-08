using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, StaticArrays

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
        Mmat[:,i] = outs[namem]
        # alphamat[:,i] = outs[namealpha]
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

t0 = tvec[1]
g=inputfile["Gravity"]
azimuth0 = azimuth[1]

distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(rvec)}}()
Rotors.update_forces!(distributed_loads, fxmat[1,:], fymat[1,:], Mmat[1,:], rvec, assembly)


Omega = SVector(0.0, 0.0, -env.RS(t0)) #Todo: This might be spinning the wrong way. 
gravity0 = SVector(g*sin(azimuth0), -g*cos(azimuth0), 0.0)
prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed
# V0 = [SVector(0.0, -rgxp[i]*env.RS(t0), 0.0) for i in eachindex(assembly.points)]


system, history0, converged = GXBeam.time_domain_analysis(assembly, tvec[1:1]; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, angular_velocity = Omega, gravity=gravity0) 
    

Vtot = getproperty.(history0[1].points, :V)
Vy = [Vtot[i][2] for i in eachindex(Vtot)]



nothing