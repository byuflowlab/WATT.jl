using Plots, OpenFASTsr, DelimitedFiles, Rotors, DynamicStallModels


#=
Test the loose coupling (like how OpenFAST couples the BEM and dynamic stall models).

Adam Cardoza 8/9/22
=#

path = dirname(@__FILE__)
cd(path)

of = OpenFASTsr
DS = DynamicStallModels

# include("../src/blades.jl")
# include("../src/environments.jl")
# include("../src/riso.jl")
# include("../src/solvers.jl")
# include("../src/utils.jl")
# include("../src/loosely.jl")



### Read in AeroDyn files
addriver = of.read_addriver("NREL5MW_ADdriver.inp", "./OpenFAST_NREL5MW")
adfile = of.read_adfile("NREL5MW_ADfile.dat","./OpenFAST_NREL5MW/")
adblade = of.read_adblade("NREL5MW_adblade.dat", "./OpenFAST_NREL5MW")


#### Define variables. 
rhub = addriver.hubrad
rtip = addriver.hubrad + adblade.span[end]

indices = 5:length(adblade.span)-1 #Skip the indices with cylinders cause that causes problems with the dynamic stall model. #Note: For some odd reason the final point was giving me crap... oh yeah, the tip correction drives things to zero. 
rvec = adblade.span[indices] .+ rhub #[11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.49].+rhub
chordvec = adblade.chord[indices] #[4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = (pi/180) .* adblade.twist[indices] #[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

B = 3
hubht = 90.0


## Airfoil Constants
# A = [0.3, 0.7] #Aerodyn defaults
A = [0.29, 0.33]
b = [0.13, 0.53]
Tp = 1.7
Tf = 3.0


## Turbine Control variables
pitch = addriver.pitch[1]
precone = addriver.precone*pi/180 #0.0 #2.5*pi/180 #TODO: !!!! I need to work in a way to include precone
yaw = addriver.yaw[1]*(pi/180) # 0.0*pi/180
tilt = addriver.shfttilt*(pi/180) #0.0 #5.0*pi/180
azimuth = 0.0

## Environmental variables
vinf = addriver.windspeed[1] #10.0
# tsr = 7.55
# rotorR = rtip*cos(precone)
rpm = addriver.rpm[1]
omega = rpm*(2*pi)/60 #vinf*tsr/rotorR

rho = adfile.airdens #1.225
mu = adfile.kinvisc #1.464e-5 #18.13e-6
a = adfile.spdsound #343.0
shearexp = addriver.shearexp[1] #0.0


### Prep the ASD rotor and operating conditions 
aftypes = Array{Array{Float64, 2}}(undef, 8)
aftypes[1] = readdlm("./OpenFAST_NREL5MW/Airfoils/Cylinder1.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder1.dat", radians=false)
aftypes[2] = readdlm("./OpenFAST_NREL5MW/Airfoils/Cylinder2.dat", skipstart=55) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder2.dat", radians=false)
aftypes[3] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU21_A17.dat", skipstart=54) # AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU40_A17.dat", radians=false)
aftypes[4] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU25_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU35_A17.dat", radians=false)
aftypes[5] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU30_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU30_A17.dat", radians=false)
aftypes[6] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU35_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU25_A17.dat", radians=false)
aftypes[7] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU40_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU21_A17.dat", radians=false)
aftypes[8] = readdlm("./OpenFAST_NREL5MW/Airfoils/NACA64_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = adblade.afid[indices] #[3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]


# create airfoil array
airfoils = aftypes[af_idx]

n = length(rvec)
afs = Array{DS.Airfoil}(undef, n)

for i = 1:n
    localpolar = hcat(airfoils[i][:,1].*(pi/180), airfoils[i][:,2:end])
    afs[i] = airfoil(localpolar; A=A)
end

rR = rvec./rtip
blade = Rotors.Blade(rhub, rtip, rR, afs) #Todo: Weird that this won't export. 



env = environment(rho, mu, a, vinf, omega, shearexp)


#### Define solution
tspan = (0.0, addriver.tmax[1])
dt = 0.01

tvec = tspan[1]:dt:tspan[2]



loads, cchistory, xds = simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade, env, tvec; verbose=true) #1.234718 seconds (5.51 M allocations: 262.152 MiB, 4.71% gc time, 83.64% compilation time) It's the createrisoode function that's taking time I think. There is probably a better way to do that. 


# Cl, Cd = parsesolution(xds, cchistory, tvec, rvec, chordvec, twistvec, pitch, blade, env)


#### Read in the OpenFAST solution (Must be run seperately)
filename = "./OpenFAST_NREL5MW/NREL5MW_aeroonly.1.out"

fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names))

tvec_of = outs["Time"]

nt_o = length(tvec_of)
na_o = length(adblade.span)


fxmat = zeros(nt_o, na_o)
fymat = zeros(nt_o, na_o)

for i = 1:na_o
    if i<10
        number = "00$i"
    elseif i<100
        number = "0$i"
    else
        number = "$i"
    end
    namex = "AB1N"*number*"Fx"
    namey = "AB1N"*number*"Fy"
    fxmat[:,i] = outs[namex]
    fymat[:,i] = outs[namey]
end


# @show loads.Fx[end, end]
# @show fxmat[end,end-1]
# println("")
# @show loads.Fy[end, end]
# @show fymat[end,end-1]

### Tip loading 
# -> Simulation isn't running OpenFAST inputs

tipplt = plot(xaxis="Time (s)", yaxis="Distributed Load (N/m)", ylims=(100, 7000))
# plot!(tvec, loads.N[:,end], lab="Normal - Unsteady")
# plot!(tvec, loads.T[:,end], lab="Tangent - Unsteady")
# hline!([cchistory[1].Np[end]], lab="Normal - Steady")
# hline!([cchistory[1].Tp[end]], lab="Tangent - Steady")
plot!(tvec, loads.Fx[:,end], lab="Fx - Unsteady")
plot!(tvec, loads.Fy[:,end], lab="Fy - Unsteady")
plot!(tvec_of, fxmat[:,end-1], markershape=:x, lab="Fx - OpenFAST")
plot!(tvec_of, fymat[:,end-1], markershape=:x, lab="Fy - OpenFAST")
display(tipplt)









nothing