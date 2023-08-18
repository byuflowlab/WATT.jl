using Revise
using OpenFASTsr, DelimitedFiles, GXBeam, Rotors, LinearAlgebra, DynamicStallModels
using FiniteDiff, ForwardDiff
# using Infiltrator
using StaticArrays
using Plots, LaTeXStrings
using Statistics

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
    
    fullouts = readdlm("./sn5_input.out", skipstart=6)

    names = fullouts[1,:]


    data = Float64.(fullouts[3:end,:])

    outs = Dict(names[i] => data[:,i] for i in eachindex(names))

    tvec = outs["Time"]
    nt = length(tvec)

    readflag = false
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
blade0 = Rotors.Blade(rvec, twistvec.*(pi/180), airfoils; rhub=rhub, rtip=rtip, precone)

# aerostates, gxstates, mesh = initialize(blade0, assembly, tvec; verbose=false)


turbine = true
rotor = Rotors.Rotor(Int(B), hubht, turbine; tilt, yaw)

azimuth0 = 0.0

nstep = 10

if nstep>nt
    println("There aren't that many time steps!.")
    nstep = nt
end

Mxr = zeros(nstep)

function objective(x)
    # println("I got called")

    blade = Rotors.Blade(rvec, x.*(pi/180), airfoils; rhub=rhub, rtip=rtip, precone)
    global aerostates, gxstates, mesh


    
    system = initial_condition!(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch)

    for i = 2:nstep
        system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch)
    end

    for i = 1:nstep
        Mxr[i] = of.root_bending_moment(rvec, aerostates.fx[i,:])
    end
    #Only calculate Mxr at nstep is an instantaneous derivative. 

    return mean(Mxr)
end

function objective2(x)
    # println("I got called")

    # @show typeof(x) #As expected. #Todo: The initialize phase isn't included in this.... bummer. 

    blade = Rotors.Blade(rvec, x.*(pi/180), airfoils; rhub=rhub, rtip=rtip, precone)

    # @show typeof(blade.twist)

    aerostates, gxstates, mesh = initialize(blade, assembly, tvec; verbose=false)

    
    system = initial_condition!(rotor, blade, assembly, env, aerostates, gxstates, mesh, tvec, azimuth0, pitch)

    for i = 2:nstep
        system = take_step!(aerostates, gxstates, mesh, rotor, blade, assembly, env, system, tvec, i, pitch)
    end

    for i = 1:nstep
        Mxr[i] = of.root_bending_moment(rvec, aerostates.fx[i,:])
    end
    #Only calculate Mxr at nstep is an instantaneous derivative. 

    return mean(Mxr)
end

# function grad_c(f, x, h=1e-100)
#     n = length(x)
#     dfdx = zero(x)

#     xcurr = complex.(x)

#     for i=1:n
#         xcurr[i] = xcurr[i] + h*im
#         dfdx[i] = image(f(xcurr))/h
#         xcurr[i] = complex(x[i])
#     end

#     return dfdx
# end

# deriv_cstep(f, x, h=1e-100) = image(f(complex(x, h)))/h

x0 = twistvec

# Mxravg = objective(x0)


# dMdx = zero(x0)


dMavgdx_fm = ForwardDiff.jacobian(objective2, x0)
# dMavgdx_fd = FiniteDiff.finite_difference_gradient(objective, x0)

# dMavgdx_c = grad_c(objective, x0) #Rotors doesn't appear to be compatible with complex numbers currently. 

println("dMdx_avg:")
for i = 1:length(x0)
    println(i, ": ", dMavgdx_fd[i], ",  ", dMavgdx_fm[i])
end
println("")

#Todo. How do I know that a derivative is accurate? -> I guess if finite differentiation and AD match.

nothing