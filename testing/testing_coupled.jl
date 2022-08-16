using DifferentialEquations, Plots, StaticArrays, OpenFASTsr, DelimitedFiles, dynamicstallmodels, Rotors, Plots.PlotMeasures

DS = dynamicstallmodels
of = OpenFASTsr


path = dirname(@__FILE__)
cd(path)



ofpath = "./OpenFAST_NREL5MW" 

### Read in AeroDyn files
addriver = of.read_addriver("NREL5MW_ADdriver.inp", "./OpenFAST_NREL5MW")
adfile = of.read_adfile("NREL5MW_ADfile.dat","./OpenFAST_NREL5MW/")
adblade = of.read_adblade("NREL5MW_adblade.dat", "./OpenFAST_NREL5MW")
edfile = of.read_edfile("NREL5MW_edfile.dat", ofpath)
# bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
bdblade = of.read_bdblade("NREL5MW_bdblade.dat", ofpath)


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
afs = Array{DS.Airfoil}(undef, n)

for i = 1:n
    localpolar = hcat(airfoils[i][:,1].*(pi/180), airfoils[i][:,2:end])
    afs[i] = complexairfoil(localpolar; A=A)
end

rR = rvec./rtip
blade = Rotors.Blade(rhub, rtip, rR, afs)



env = Rotors.environment(rho, mu, a, vinf, omega, shearexp)

assembly = of.make_assembly(edfile, bdblade)
ne = length(assembly.elements)
#Todo: Well, something exploded. -> Potentially have a positive feedback loop

#=
Okay, yeah it has what looks like a positive feedback loop. 

What's interesting is that the tip appears to have 20+ degrees of deflection... which seems like a lot. 

It doesn't instantly shoot off, it shoots off after a period of time. So it could be due to a lack of numerical damping? 

Okay, I found two things:
     Todo: One, I have no rotational damping.
     
     Todo: Two, The thetas that GXBeam provides in the state are not actually the angle, but the components of a vector that the beam is rotating about. 

     I got rid of the deflection on twist in the meantime.

     - with b = 100, the damping forces aren't very large. They do seem fairly erratic, but I can't tell for sure. They might be fairly smooth. 

     - with b = 500, the solution blows up... so that's no bueno. That implies that the structural damping is actually kind of bad for numerical behaivor. 

     - reverting b back to 0.01 makes the simulation blow up again. 

     I gues that didn't fix the problem. So that tells us that the problem is in the model. I could push the step size smaller and see what happens. I could also try GXBeam with the new updates (we'd definitely want to test it before we put that into effect.)

    - Maybe I need to return the velocities. (Taylor confirmed that they are returned in the inertial frame. )

    - I wonder if part of it is the ode solver I'm using. 
    -> I figured out using callbacks to change the parameters.... Which is cool. So that's another method I can use for integration now. But it looks like the method is stalling out, at least it's probably the nonlinear solve that's bottomed out. So... that's a bummer. 

    Oh, it looks like the ODE solver bottomed out... not the nonlinear solver... so that's a bummer. 

    Yeah, I think it's a problem either with GXBeam going bonkers... or the velocities... or something. 
=#


#### Define solution
tspan = (0.0, addriver.tmax[1])
dt = 0.005 #0.01

tvec = tspan[1]:dt:tspan[2]



loads, cchistory, xds, gxhistory = Rotors.simulate(rvec, chordvec, twistvec, rhub, rtip, hubht, B, pitch, precone, tilt, yaw, blade, env, assembly, tvec; verbose=true, dsmodelinit=Rotors.Steady(), b=0.01, solver=Rotors.RK4()) #DiffEQ(Tsit5())

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

        defthetax = [state.points[ipoint].theta[1] for ipoint = 1:length(assembly.points)].*(180/pi)
        defthetay = [state.points[ipoint].theta[2] for ipoint = 1:length(assembly.points)].*(180/pi)
        defthetaz = [state.points[ipoint].theta[3] for ipoint = 1:length(assembly.points)].*(180/pi)

        plot(leg=:topleft, xaxis="Beam length (m)", yaxis="Deflection (m)", xlims=(-0.005, rtip*1.02), ylims=(-12, 4), right_margin=20mm) #
        plot!(x, deflectionx, lab="X")
        plot!(x, deflectiony, lab="Y")
        plot!(x, deflectionz, lab="Z")
        annotate!((20, 2.5, text("time = $t sec", :left)))

        # subplot = twinx()
        # plot!(subplot, leg=:bottomleft, xlims=(-0.005, rtip*1.02), ylims=(-15, 35), yaxis="Angular Deflection (deg)")
        # plot!(subplot, x, defthetax, lab="θx", linestyle=:dash)
        # plot!(subplot, x, defthetay, lab="θy", linestyle=:dashdot)
        # plot!(subplot, x, defthetaz, lab="θz", linestyle=:dot)
    end every 10
    gif(anim, "bodydeflections_looselycoupled_openfast_081222.gif", fps = framespersecond)
end

# createdeflectionanimation()

# nothing