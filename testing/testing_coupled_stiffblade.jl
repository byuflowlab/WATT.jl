using DifferentialEquations, FLOWMath, CCBlade, GXBeam, LinearAlgebra, Plots, StaticArrays, CurveFit, NLsolve, OpenFASTsr, DelimitedFiles

include("../src/blades.jl")
include("../src/environments.jl")
# include("../dev/bem.jl")
include("../src/riso.jl")
# include("../dev/bem-riso.jl")
include("../dev/gxbeam.jl")
include("../src/solvers.jl")
include("../src/loosely.jl")
include("../src/coupled.jl")

include("../src/extra.jl")
# include("../src/static.jl")

of = OpenFASTsr
ofpath = "./OpenFAST_NREL5MW" 


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
thickvec = chordvec.*0.08

B = 3
hubht = 90.0


## Airfoil Constants
# A = [0.3, 0.7] #Aerodyn defaults
A = [0.29, 0.33] #Hansen defaults
# A = [0.2, 0.7] # Random -> The final loading seems pretty insensitive to these. 
b = [0.13, 0.53] #Hansen defaults? 
# b = [0.01, 0.1] #Random -> The final loading also seems pretty insensitive to these. 
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
afs = Array{Airfoil}(undef, n)

for i = 1:n
    localpolar = hcat(airfoils[i][:,1].*(pi/180), airfoils[i][:,2:end])
    afs[i] = complexairfoil(localpolar; A=A)
end

rR = rvec./rtip
blade = Blade(rhub, rtip, rR, afs)



env = environment(rho, mu, a, vinf, omega, shearexp)



p_s, points, elements = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)
gxmodel = gxbeam(points, elements)


## Create system indices  
start = 1:gxmodel.ne
stop = 2:gxmodel.np

## Create assembly
assembly = create_gxbeam_assembly(gxmodel, p_s, start, stop) 




#### Define solution
tspan = (0.0, addriver.tmax[1])
dt = 0.01

tvec = tspan[1]:dt:tspan[2]



loads, cchistory, xds, gxhistory = simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, precone, tilt, yaw, blade, env, assembly, tvec; verbose=true, dsmodelinit=Steady(), solver=DiffEQ(Tsit5()))

#WorkLocation: 

#=
Todo.  I need to check that the motion of the deflections matches up with increases and decreases in the velocities that the blade sees. -> I'm not sure if it matches. -> I realized I needed to pull the precone, pitch, tilt, and shear differences out to see the motion differences. And it looks like it's correct. 

Todo: I need to compare the deflections of Rotors.jl to OpenFAST -> Need to test a more flexible blade.

TODO: I'd like to make a function that creates assemblies for the end user. 
=#

nt = length(tvec)


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




### Tip loading 
tipplt = plot(xaxis="Time (s)", yaxis="Distributed Load (N/m)", leg=:bottomright) #, ylims=(100, 7000)
plot!(tvec, loads.Fx[:,end], lab="Fx - Unsteady")
plot!(tvec, loads.Fy[:,end], lab="Fy - Unsteady")
plot!(tvec_of, fxmat[:,end-1], markershape=:x, lab="Fx - OpenFAST")
plot!(tvec_of, fymat[:,end-1], markershape=:x, lab="Fy - OpenFAST")
display(tipplt)

timeelapsed = tspan[2]-tspan[1]
framespersecond = round(Int, nt/timeelapsed)*10

function createdeflectionanimation()
    anim = @animate for i in 1:nt
        t = tvec[i]
        state = gxhistory[i]
        x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
    
        deflectionx = [state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
        deflectiony = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]
        deflectionz = [state.points[ipoint].u[3] for ipoint = 1:length(assembly.points)]

        plot(leg=:topleft, xaxis="Beam length (m)", yaxis="Deflection (m)", xlims=(-0.005, rtip*1.02), ylims=(-0.05, 0.09))
        plot!(x, deflectionx, lab="X")
        plot!(x, deflectiony, lab="Y")
        plot!(x, deflectionz, lab="Z")
        annotate!((5, 0.05, text("time = $t sec", :left)))
    end #every 10
    gif(anim, "bodydeflections_looselycoupled_stiff_081122.gif", fps = framespersecond)
end


nothing