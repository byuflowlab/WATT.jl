using Plots, StaticArrays, OpenFASTsr, DelimitedFiles, DynamicStallModels, Rotors, LaTeXStrings

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
bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
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







env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)

assembly = of.make_assembly(edfile, bdfile, bdblade)
ne = length(assembly.elements)







#### Read in the OpenFAST solution (Must be run seperately)
filename = "./OpenFAST_NREL5MW_modified/NREL5MW_input.out"

fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in eachindex(names))

tvec = outs["Time"]

nt = length(tvec)
na = length(adblade["BlSpn"])


fxmat = zeros(nt, na)
fymat = zeros(nt, na)
Mmat = zeros(nt, na)

for i = 1:na
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
    fxmat[:,i] = outs[namex]
    fymat[:,i] = outs[namey]
    Mmat[:,i] = outs[namem]
    
end

azimuth = outs["B1Azimuth"]
delz = -outs["B1TipTDxr"] #Translate to my structural frame. 
dely = outs["B1TipTDyr"]
delx = outs["B1TipTDzr"]


#### Define solution

gxhistory = Rotors.simulate_gxbeam(rvec, rhub, rtip, tvec, azimuth, fxmat, fymat, Mmat, env, assembly; verbose=true, speakiter=100, structural_damping=true, linear=false, g=inputfile["Gravity"])

    

#Tip deflections
tipdef_x = [gxhistory[i].points[end].u[1] for i in eachindex(tvec)]
tipdef_y = [gxhistory[i].points[end].u[2] for i in eachindex(tvec)]
tipdef_z = [gxhistory[i].points[end].u[3] for i in eachindex(tvec)]


tipdefplt = plot(xaxis="Time (s)", yaxis="Deflection (m)", leg=:outerright)
plot!(tvec, tipdef_x, lab=L"$\delta_x$ GXBeam")
plot!(tvec, tipdef_y, lab=L"$\delta_y$ GXBeam")
plot!(tvec, tipdef_z, lab=L"$\delta_z$ GXBeam")
plot!(tvec, delx, lab=L"$\delta_x$ OpenFAST", linestyle=:dash)
plot!(tvec, dely, lab=L"$\delta_y$ OpenFAST", linestyle=:dash)
plot!(tvec, delz, lab=L"$\delta_z$ OpenFAST", linestyle=:dash)
display(tipdefplt)

#=
Well, at this point. I'm pretty sure that the stiffness matrix is correct. And I'm pretty sure I'm rotating the loads into the structural frame correctly... So what are the other options of what is wrong? I don't know. 
=#








nothing
