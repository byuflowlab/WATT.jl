using Plots, StaticArrays, OpenFASTsr, DelimitedFiles, DynamicStallModels, Rotors

# println("Finished loading packages.")

DS = DynamicStallModels
of = OpenFASTsr


path = dirname(@__FILE__)
cd(path)



ofpath = "./OpenFAST_NREL5MW_modified" 

### Read in OpenFAST files
inputfile = of.read_inputfile("NREL5MW_input.fst", ofpath)
inflowwind = of.read_inflowwind("NREL5MW_inflowwind.dat", ofpath)
adfile = of.read_adfile("NREL5MW_ADfile.dat", ofpath)
adblade = of.read_adblade("NREL5MW_adblade.dat", ofpath)
edfile = of.read_edfile("NREL5MW_edfile.dat", ofpath)
# bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
bdblade = of.read_bdblade("NREL5MW_bdblade.dat", ofpath)


#### Define variables. 
rhub = edfile["HubRad"]
rtip = rhub + adblade["BlSpn"][end]

indices = 1:length(adblade["BlSpn"]) # 5:length(adblade.span)-1 #Skip the indices with cylinders cause that causes problems with the dynamic stall model. #Note: For some odd reason the final point was giving me crap... oh yeah, the tip correction drives things to zero....  
rvec = adblade["BlSpn"][indices] .+ rhub #[11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.49].+rhub
chordvec = adblade["BlChord"][indices] #[4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = (pi/180) .* adblade["BlTwist"][indices] #[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

rvec[1] += 0.001
rvec[end] -= 0.001

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
pitch = edfile["BlPitch(1)"]
precone = edfile["PreCone(1)"]*pi/180 #0.0 #2.5*pi/180 #TODO: !!!! I need to work in a way to include precone
yaw = edfile["NacYaw"]*(pi/180) # 0.0*pi/180
tilt = edfile["ShftTilt"]*(pi/180) #0.0 #5.0*pi/180
azimuth = 0.0

## Environmental variables
vinf = inflowwind["HWindSpeed"] #10.0
# tsr = 7.55
# rotorR = rtip*cos(precone)
rpm = edfile["RotSpeed"]
omega = rpm*(2*pi)/60 #vinf*tsr/rotorR

rho = inputfile["AirDens"] #1.225
mu = inputfile["KinVisc"] #1.464e-5 #18.13e-6
a = inputfile["SpdSound"] #343.0
shearexp = inflowwind["PLexp"] #0.0


### Prep the ASD rotor and operating conditions 
aftypes = Array{of.AirfoilInput}(undef, 8)
aftypes[1] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/Cylinder1.dat") #readdlm("./OpenFAST_NREL5MW_modified/Airfoils/Cylinder1.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder1.dat", radians=false)
aftypes[2] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/Cylinder2.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/Cylinder2.dat", skipstart=55) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder2.dat", radians=false)
aftypes[3] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU40_A17.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/DU21_A17.dat", skipstart=54) # AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU40_A17.dat", radians=false)
aftypes[4] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU35_A17.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/DU25_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU35_A17.dat", radians=false)
aftypes[5] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU30_A17.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/DU30_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU30_A17.dat", radians=false)
aftypes[6] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU25_A17.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/DU35_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU25_A17.dat", radians=false)
aftypes[7] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/DU21_A17.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/DU40_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU21_A17.dat", radians=false)
aftypes[8] = of.read_airfoilinput("./OpenFAST_NREL5MW_modified/Airfoils/NACA64_A17.dat") # readdlm("./OpenFAST_NREL5MW_modified/Airfoils/NACA64_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = Int.(adblade["BlAFID"][indices]) #[3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]


# create airfoil array
afs = aftypes[af_idx]

n = length(rvec)
airfoils = Vector{DS.Airfoil}(undef, n)
for i = 1:n
    airfoils[i] = make_dsairfoil(afs[i])
end

rR = rvec./rtip
blade = Rotors.Blade(rhub, rtip, rR, airfoils)



env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)

assembly = of.make_assembly(edfile, bdblade)
ne = length(assembly.elements)



dsmodel = DS.BeddoesLeishman(DS.Indicial(), n, airfoils, 3)
dsmodelinit = Rotors.BeddoesLeishman()


#### Define solution
tspan = (0.0, inputfile["TMax"]) #(0.0, 4.6)
dt = inputfile["DT"][1] #0.01

tvec = tspan[1]:dt:tspan[2]

# println("Got to the simulation.")

loads, cchistory, xds, gxhistory, def_thetax = Rotors.simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade, env, assembly, tvec; verbose=true, dsmodel, dsmodelinit, b=0.01, solver=Rotors.RK4(), speakiter=1000) #DiffEQ(Tsit5())



nt = length(tvec)


#### Read in the OpenFAST solution (Must be run seperately)
filename = "./OpenFAST_NREL5MW_modified/NREL5MW_aeroonly.1.out"
filename = "./OpenFAST_NREL5MW_modified/NREL5MW_input.out"

fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in eachindex(names))

tvec_of = outs["Time"]

nt_o = length(tvec_of)
na_o = length(adblade["BlSpn"])


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
tipplt = plot(xaxis="Time (s)", yaxis="Tip Load (N/m)", leg=:topright) #, ylims=(100, 7000)
plot!(tvec, loads.Fx[:,end], lab="Fx - Unsteady")
plot!(tvec, loads.Fy[:,end], lab="Fy - Unsteady")
plot!(tvec_of, fxmat[:,end], lab="Fx - OpenFAST")
plot!(tvec_of, fymat[:,end], lab="Fy - OpenFAST")
display(tipplt)

#=
- The current difference could be due to the fact that I'm not including angular deflection in the loop.

-> well, I added it in. I'm not sure that I added it in correctly. But that's just not trusting what I did in the past. Additionally, I've added some extra bugs from changing how I iterate through time. So I need to get those fixed. The solution was crazy unstable. 

-> I don't know if this makes me feel any better, but I went and tried to run the OpenFAST and it isn't able to handle it either... which is confusing. Because I'm pretty sure I ran this previously. So... that's grand. I'll have to see if the original files will run. 

-> Well.... the original models are turning off unsteady aerodynamics... so... that's great. I wonder if there is a different model that I can use. 

-> Well the trick to handling the instability was decreasing the timestep on the aero solver for OpenFAST, and that did the trick here as well. Now the main difference is the fact that OpenFAST accounts for moment loads and gravity, while mine doesn't. 

=#


timeelapsed = tspan[2]-tspan[1]
framespersecond = round(Int, nt/timeelapsed)*10

# function createdeflectionanimation()
#     anim = @animate for i in 1:nt
#         t = tvec[i]
#         state = gxhistory[i]
#         x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
    
#         deflectionx = [state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
#         deflectiony = [state.points[ipoint].u[2] for ipoint = 1:length(assembly.points)]
#         deflectionz = [state.points[ipoint].u[3] for ipoint = 1:length(assembly.points)]

#         defthetax = [state.points[ipoint].theta[1] for ipoint = 1:length(assembly.points)].*(180/pi)
#         defthetay = [state.points[ipoint].theta[2] for ipoint = 1:length(assembly.points)].*(180/pi)
#         defthetaz = [state.points[ipoint].theta[3] for ipoint = 1:length(assembly.points)].*(180/pi)

#         plot(leg=:topleft, xaxis="Beam length (m)", yaxis="Deflection (m)", xlims=(-0.005, rtip*1.02), ylims=(-12, 4), right_margin=20mm) #
#         plot!(x, deflectionx, lab="X")
#         plot!(x, deflectiony, lab="Y")
#         plot!(x, deflectionz, lab="Z")
#         annotate!((20, 2.5, text("time = $t sec", :left)))

#         # subplot = twinx()
#         # plot!(subplot, leg=:bottomleft, xlims=(-0.005, rtip*1.02), ylims=(-15, 35), yaxis="Angular Deflection (deg)")
#         # plot!(subplot, x, defthetax, lab="θx", linestyle=:dash)
#         # plot!(subplot, x, defthetay, lab="θy", linestyle=:dashdot)
#         # plot!(subplot, x, defthetaz, lab="θz", linestyle=:dot)
#     end every 10
#     gif(anim, "bodydeflections_looselycoupled_openfast_081222.gif", fps = framespersecond)
# end

# createdeflectionanimation()

nothing
