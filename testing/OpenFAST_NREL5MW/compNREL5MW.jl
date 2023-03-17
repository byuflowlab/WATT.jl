#=
Compare CCBlade's solution against OpenFAST's solution of the NREL 5MW wind turbine. 

=#

using CCBlade, Plots, DelimitedFiles, FLOWMath

path = dirname(@__FILE__)
cd(path)

data = [0.0000000E+00  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.3308000E+01  3.5420000E+00        1
1.3667000E+00  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.3308000E+01  3.5420000E+00        1
4.1000000E+00  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.3308000E+01  3.8540000E+00        1
6.8333000E+00  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.3308000E+01  4.1670000E+00        2
1.0250000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.3308000E+01  4.5570000E+00        3
1.4350000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.1480000E+01  4.6520000E+00        4
1.8450000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.0162000E+01  4.4580000E+00        4
2.2550000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  9.0110000E+00  4.2490000E+00        5
2.6650000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  7.7950000E+00  4.0070000E+00        6
3.0750000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  6.5440000E+00  3.7480000E+00        6
3.4850000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  5.3610000E+00  3.5020000E+00        7
3.8950000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  4.1880000E+00  3.2560000E+00        7
4.3050000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  3.1250000E+00  3.0100000E+00        8
4.7150000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  2.3190000E+00  2.7640000E+00        8
5.1250000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.5260000E+00  2.5180000E+00        8
5.4666700E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  8.6300000E-01  2.3130000E+00        8
5.7400000E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  3.7000000E-01  2.0860000E+00        8
6.0133300E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.0600000E-01  1.4190000E+00        8
6.1499900E+01  0.0000000E+00  0.0000000E+00 0.0000000E+00  1.0600000E-01  1.4190000E+00        8]

rhub = 1.5
rvec = data[:,1] .+ rhub
cvec = data[:,6]
twistvec = data[:,5]
afid = Int.(data[:,7])
rtip = rvec[end] #63
hubht = 89.56
B = 3

precone = 0.0
tilt = 0.0
yaw = 0.0
shearexp = 0.0
azimuth = 0.0
pitch = 0.0

rho = 1.225
vinf = 10.0

rpm = 11.443998
omega = rpm*2*pi/60

afmats = Array{Array{Float64, 2}, 1}(undef, 8)
afmats[1] = readdlm("./OpenFAST_NREL5MW/Airfoils/Cylinder1.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder1.dat", radians=false)
afmats[2] = readdlm("./OpenFAST_NREL5MW/Airfoils/Cylinder2.dat", skipstart=55) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder2.dat", radians=false)
afmats[3] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU21_A17.dat", skipstart=54) # AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU40_A17.dat", radians=false)
afmats[4] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU25_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU35_A17.dat", radians=false)
afmats[5] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU30_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU30_A17.dat", radians=false)
afmats[6] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU35_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU25_A17.dat", radians=false)
afmats[7] = readdlm("./OpenFAST_NREL5MW/Airfoils/DU40_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU21_A17.dat", radians=false)
afmats[8] = readdlm("./OpenFAST_NREL5MW/Airfoils/NACA64_A17.dat", skipstart=54) #AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/NACA64_A17.dat", radians=false)

aftypes = [AlphaAF(afmats[i][:,1].*(pi/180), afmats[i][:,2], afmats[i][:,3]) for i = 1:8] #Okay, I'm pretty sure that 

airfoils = aftypes[afid]




rotor = Rotor(rhub, rtip, B, precone=precone, turbine=true, tip=nothing)

sections = Section.(rvec, cvec, twistvec.*(pi/180), airfoils)

op = windturbine_op.(vinf, omega, pitch, rvec, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

ccout = solve.(Ref(rotor), sections, op)







#### Read in the OpenFAST solution (Must be run seperately)
filename = "./OpenFAST_NREL5MW/NREL5MW_aeroonly.1.out"

fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names))

tvec_of = outs["Time"]

nt_o = length(tvec_of)
na_o = length(rvec) #length(adblade.span)


fxmat = zeros(nt_o, na_o)
fymat = zeros(nt_o, na_o)
phimat = zeros(nt_o, na_o)
alphamat = zeros(nt_o, na_o)

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
    namephi = "AB1N"*number*"Phi"
    namealpha = "AB1N"*number*"Alpha"

    fxmat[:,i] = outs[namex]
    fymat[:,i] = outs[namey]
    phimat[:,i] = outs[namephi]
    alphamat[:,i] = outs[namealpha]
end

phiplt = plot(xaxis="Radius (m)", yaxis="Inflow Angle (deg)")
vspan!([rvec[1], (rvec[4]+rvec[5])/2], fillalpha=0.3, lab="Cylinder Section")
plot!(rvec, ccout.phi.*(180/pi), lab="CCBlade", markershape=:cross)
plot!(rvec, phimat[end,:], lab="OpenFAST", markershape=:x)
display(phiplt)

#=
The first 4 points don't make sense, since they are all cylinders. ... There are still errors, which mean that the airfoil polars are off. Which can't be possible... Unless OpenFAST is doing some kind of treatment to the polars. 
=#
idx = 9
polar = afmats[afid[idx]]
plrplt = plot(xaxis="Angle of Attack (deg)", yaxis="Coefficient", title="Polar $idx")
plot!(polar[:,1], polar[:,2], lab="Lift")
plot!(polar[:,1], polar[:,3], lab="Drag")
vline!([ccout.alpha[idx]*180/pi], lab="CCBlade")
vline!([alphamat[1, idx]], lab="OpenFAST")
display(plrplt)

#=
Most of the polars aren't smooth. I wonder if Brent's method is what is causing the difference. Let's plot the residual and see what the space looks like. 
=#

phi = (-179:1:179).*(pi/180)
resids = zero(phi)

for i = 1:length(phi)
    resids[i], _ = CCBlade.residual(phi[i], rotor, sections[idx], op[idx])
end


#Todo. For some reason phi is not showing up. The inflow angle is found, and I plot that region, but the plot doesn't cross zero.. or show up at all. -> Silly changing the calculated angle to -179:179 fixed the problem. So some interesting plot behavior from Plots.jl
phiplt = plot(xaxis="Inflow Angle (deg)", yaxis="Residual", title="Section $idx", xlims=(0,1), ylims=(-2, 2), leg=:bottomright)
plot!(phi, resids, lab="residual()")
vline!([ccout.phi[idx]], lab="CCBlade")
vline!([phimat[1,idx]*(pi/180)], lab="OpenFAST")
display(phiplt)

@show ccout.alpha[idx], alphamat[1, idx]*(pi/180)

#=
Okay... So we know that the initial problem is something to do with the residual values. So either Brent's method is having an issue, or the residual functions are different.
Todo.
    - check that OpenFAST is using the single residual formulation.
    - check the method that OpenFAST is using to solve the residual. 
    - Try a set of smooth polars for both methods. 

In the future I may need to have options of how to solve the residual. Either by brute force, Newton's method, Brent's method, etc (a non-linear solver?). 

I hadn't converted my twist vector to radians. 

It looks. like CCBlade automatically sets the root and tip values to zero. There also appears to be a coouple of minor discrepancies. 
=#