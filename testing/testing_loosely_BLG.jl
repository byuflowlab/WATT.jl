using Plots, OpenFASTsr, DelimitedFiles, Rotors, DynamicStallModels, CCBlade, LaTeXStrings


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
addriver = of.read_addriver("NREL5MW_ADdriver.dvr", "./OpenFAST_NREL5MW_modified")
adfile = of.read_adfile("NREL5MW_ADfile.dat","./OpenFAST_NREL5MW_modified/")
adblade = of.read_adblade("NREL5MW_adblade.dat", "./OpenFAST_NREL5MW_modified")


#### Define variables. 
rhub = addriver["HubRad(1)"]
rtip = rhub + adblade["BlSpn"][end]

indices = 1:length(adblade["BlSpn"]) # 5:length(adblade.span)-1 #Skip the indices with cylinders cause that causes problems with the dynamic stall model. #Note: For some odd reason the final point was giving me crap... oh yeah, the tip correction drives things to zero....  
rvec = adblade["BlSpn"][indices] .+ rhub #[11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.49].+rhub
chordvec = adblade["BlChord"][indices] #[4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = (pi/180) .* adblade["BlTwist"][indices] #[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

rvec[1] += 0.001
rvec[end] -= 0.001

B = addriver["NumBlades(1)"]
hubht = addriver["HubHt(1)"]


## Airfoil Constants
 

## Turbine Control variables
pitch = addriver["Pitch_mat"][1]
precone = addriver["Precone(1)"]*pi/180 #0.0 #2.5*pi/180 #TODO: !!!! I need to work in a way to include precone
yaw = addriver["Yaw_mat"][1]*(pi/180) # 0.0*pi/180
tilt = addriver["ShftTilt(1)"]*(pi/180) #0.0 #5.0*pi/180
azimuth = 0.0

## Environmental variables
vinf = addriver["HWndSpeed_mat"][1] #10.0
# tsr = 7.55
# rotorR = rtip*cos(precone)
rpm = addriver["RotSpd_mat"][1]
omega = rpm*(2*pi)/60 #vinf*tsr/rotorR

rho = addriver["FldDens"] #1.225
mu = addriver["KinVisc"] #1.464e-5 #18.13e-6
a = addriver["SpdSound"] #343.0
shearexp = addriver["PLExp"][1] #0.0


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
    airfoils[i] = make_dsairfoil(afs[i], chordvec[i])
end

rR = rvec./rtip
# blade = Rotors.Blade(rhub, rtip, rR, airfoils) #TODO: Weird that this won't export. ???? What did I mean by this? 
blade = Blade(rvec, twistvec, airfoils; rhub, rtip, precone)





env = environment(rho, mu, a, vinf, omega, shearexp)

# dsmodel = DS.BeddoesLeishman(DS.Indicial(), n, airfoils, 3)
# dsmodelinit = Rotors.BeddoesLeishman()
turbine = true
rotor_r = Rotors.Rotor(B, hubht, turbine; tilt, yaw)





#### Define solution
tspan = (0.0, addriver["Tmax_mat"][1]) #(0.0, 0.02) #
dt = addriver["dT_mat"][1] #0.01

# tvec = collect(tspan[1]:dt:tspan[2])
tvec = collect(tspan[1]:dt:tspan[2]-dt)


# @show twistvec #Todo. Take Pcopy out. 
loads, coefs, cchistory, xds, azi_time = simulate(rotor_r, blade, env, tvec; verbose=true, azimuth0=0) #

loads2, coefs2, cchistory2, xds2, azi_time2 = simulate(rotor_r, blade, env, tvec; verbose=true, azimuth0=2*pi/3)

loads3, coefs3, cchistory3, xds3, azi_time3 = simulate(rotor_r, blade, env, tvec; verbose=true, azimuth0=4*pi/3)

na = length(rvec)
nt = length(tvec)

#### Read in the OpenFAST solution (Must be run seperately)
filename = "./OpenFAST_NREL5MW_modified/NREL5MW_ADdriver.1.out"

fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names))

tvec_of = outs["Time"]

nt_o = length(tvec_of)
na_o = length(adblade["BlSpn"])


fxmat = zeros(nt_o, na_o)
fymat = zeros(nt_o, na_o)
phimat = zeros(nt_o, na_o)
alphamat = zeros(nt_o, na_o)
Wmat = zeros(nt_o, na_o)
cnmat = zero(fxmat)
ccmat = zero(fymat)
Wdmat = zeros(nt_o, na_o)
Vrelmat = zeros(nt_o, na_o)

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
    nameVx = "AB1N"*number*"Vx"
    nameVy = "AB1N"*number*"Vy"
    namecn = "AB1N"*number*"Cn"
    namect = "AB1N"*number*"Ct"
    nameVdx = "AB1N"*number*"VDisx"
    nameVdy = "AB1N"*number*"VDisy"
    nameVrel = "AB1N"*number*"VRel"

    fxmat[:,i] = outs[namex]
    fymat[:,i] = outs[namey]
    phimat[:,i] = outs[namephi]
    alphamat[:,i] = outs[namealpha]
    Wmat[:,i] = @. sqrt(outs[nameVx]^2 + outs[nameVy]^2)
    cnmat[:,i] = outs[namecn]
    ccmat[:,i] = outs[namect]
    Wdmat[:,i] = @. sqrt(outs[nameVdx]^2 + outs[nameVdy]^2)
    Vrelmat[:,i] = outs[nameVrel]
end


# @show loads.Fx[end, end]
# @show fxmat[end,end-1]
# println("")
# @show loads.Fy[end, end]
# @show fymat[end,end-1]
fxsteady = [cchistory[i].Np[end]*cos(precone) for i in 1:length(tvec)]
fysteady = [cchistory[i].Tp[end]cos(precone) for i in 1:length(tvec)]



########## CCBlade simulation. 
rotor = CCBlade.Rotor(rhub, rtip, B, precone=precone, turbine=true, tip=nothing)

aft = Array{AlphaAF}(undef, 8)
aft[1] = AlphaAF(aftypes[1].aoa.*(pi/180), aftypes[1].cl, aftypes[1].cd)
aft[2] = AlphaAF(aftypes[2].aoa.*(pi/180), aftypes[2].cl, aftypes[2].cd)
aft[3] = AlphaAF(aftypes[3].aoa.*(pi/180), aftypes[3].cl, aftypes[3].cd)
aft[4] = AlphaAF(aftypes[4].aoa.*(pi/180), aftypes[4].cl, aftypes[4].cd)
aft[5] = AlphaAF(aftypes[5].aoa.*(pi/180), aftypes[5].cl, aftypes[5].cd)
aft[6] = AlphaAF(aftypes[6].aoa.*(pi/180), aftypes[6].cl, aftypes[6].cd)
aft[7] = AlphaAF(aftypes[7].aoa.*(pi/180), aftypes[7].cl, aftypes[7].cd)
aft[8] = AlphaAF(aftypes[8].aoa.*(pi/180), aftypes[8].cl, aftypes[8].cd)

# create airfoil array
# afscc = aft[adblade.afid] #aft[af_idx] 
afscc = aft[af_idx] 

# rcc = adblade.span .+ rhub
# ccc = adblade.chord
# tcc = adblade.twist .*(pi/180)

rcc = rvec
ccc = chordvec
tcc = twistvec

# define sections
sections = Section.(rcc, ccc, tcc, afscc)

op = windturbine_op.(vinf, omega, pitch, rcc, precone, yaw, tilt, azimuth, hubht, shearexp, rho)

ccout = CCBlade.solve.(Ref(rotor), sections, op)

# op = windturbine_op.(vinf, omega, pitch, rcc[5], precone, yaw, tilt, azimuth, hubht, shearexp, rho)
# ccout = solve(rotor, sections[5], op)

# T, Q = thrusttorque(rotor, sections, ccout)










### Tip loading 
# -> Simulation isn't running OpenFAST inputs

# tipplt = plot(xaxis="Time (s)", yaxis="Distributed Load (N/m)", leg=:outerright) #, ylims=(100, 7000)
# # plot!(tvec, loads.N[:,end], lab="Normal - Unsteady")
# # plot!(tvec, loads.T[:,end], lab="Tangent - Unsteady")
# # hline!([cchistory[1].Np[end]], lab="Normal - Steady")
# # hline!([cchistory[1].Tp[end]], lab="Tangent - Steady")
# plot!(tvec, loads.Fx[:,end], lab="Fx - Unsteady") #, markershape=:circle)
# plot!(tvec, loads.Fy[:,end], lab="Fy - Unsteady") #, markershape=:circle)
# plot!(tvec_of, fxmat[:,end-1], markershape=:x, lab="Fx - OpenFAST")
# plot!(tvec_of, fymat[:,end-1], markershape=:x, lab="Fy - OpenFAST")
# plot!(tvec, fxsteady, lab="Fx - CCBlade")
# plot!(tvec, fysteady, lab="Fy - CCBlade")
# display(tipplt)
# savefig("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/figures/coupling/mycoupling/loosely/aeroonly/beddoesleishmanaerodyn/tipload_shear.png")

xloadplt = plot() #yaxis="X Load (N/m)"
plot!(tvec, loads.Fx[:,end], lab="Fx - Unsteady")
plot!(tvec_of, fxmat[:,end-1], markershape=:x, lab="Fx - OpenFAST")
plot!(tvec, fxsteady, lab="Fx - CCBlade")

yloadplt = plot(xaxis="Time (s)", yaxis="Distributed Load (N/m)")
plot!(tvec, loads.Fy[:,end], lab="Fy - Unsteady")
plot!(tvec_of, fymat[:,end-1], markershape=:x, lab="Fy - OpenFAST")
plot!(tvec, fysteady, lab="Fy - CCBlade")

tipplt = plot(xloadplt, yloadplt, layout=(2,1))
# display(tipplt)

#Todo. Loads are much too high. I don't know if it's a polar input problem, a dimensionalization problem, or something else (like a states or CCBlade problem). -> Well, I've checked just CCBlade against the steady case, and they match fairly well (except for the very first and last points, but.... I think CCBlade treats them differently.) -> Means that the angle of attack that's getting passed to the DS model is either wrong, or the AOA isn't matching up across time. -> DSM wasn't handling multiple airfoils properly. 



index = 3
ofidx = index #+ 4 #Currently they are the same indices. 
nt = length(tvec)

phi = [cchistory[i, index].phi*180/pi for i = 1:nt]  

if ofidx<10
    num = "00$ofidx"
elseif ofidx<100
    num = "0$ofidx"
else
    num = "$ofidx"
end
name = "AB1N"*num*"Phi"
phiof = outs[name]

phiplt = plot(xaxis="Time (s)", yaxis="Inflow angle, ϕ (deg)", title="Station $index")
plot!(tvec, phi, lab="CCBlade", markershape=:x)
plot!(tvec_of, phiof, lab="OpenFAST")
# display(phiplt)






tix_local = 25

alpha1 = phi[tix_local]*pi/180 - twistvec[index]
alpha2 = phiof[tix_local]*pi/180 - twistvec[index]

cl1 = airfoils[index].cl(alpha1)
cl2 = airfoils[index].cl(alpha2)

# @show cl1, cl2, (cl1-cl2)/cl2 



# println("")
# # phicc = ccout.phi[ofidx]*180/pi
# phicc = ccout.phi*180/pi
# @show phi[1], phiof[1], phicc (phiof[1]-phi[1])/phiof[1]

# println("")

# alpha = cchistory[1].alpha[index]*180/pi
# alphaof = outs["AB1N"*num*"Alpha"][1]

# @show alpha, alphaof

# println("")

# theta = twistvec[index]*180/pi
# theta2 = cchistory[1].phi[index] - cchistory[1].alpha[index]*180/pi
# thetaof = outs["AB1N"*num*"Theta"][1]

# @show theta, theta2, thetaof #Twist values match decently well. 

# println("")

# u = cchistory[1].u[index]
# v = cchistory[1].v[index]
# uof = outs["AB1N"*num*"Vindx"][1]
# vof = outs["AB1N"*num*"Vindy"][1]

# @show u, v
# @show uof, vof

# println("")

# fx = cchistory[1].Np[index]*cos(precone)
# fy = cchistory[1].Tp[index]*cos(precone)
# fxof = outs["AB1N"*num*"Fx"][1]
# fyof = outs["AB1N"*num*"Fy"][1]

# @show fx, fy
# @show fxof, fyof

# println("")

# Tof = outs["RtAeroFxh"][1]



# Vx = op[ofidx].Vx
# Vy = op[ofidx].Vy
# Vxof = outs["AB1N"*num*"Vx"][1]
# Vyof = outs["AB1N"*num*"Vy"][1]

# # @show Vx, Vy #Todo. Aha! I found the culprit! Why are the Vy values so different? Psych... it's the same... I was using the wrong indices. 
# @show Vxof, Vyof

rof = adblade["BlSpn"] .+ rhub
alpharotmat = Array{Float64}(undef, nt, na)
Wrotmat = Array{Float64}(undef, nt, na)

# for i = 1:na
#     idx = 24*(i-1) + 24
#     alpharotmat[:,i] = pcopy[:,idx]

#     idx = 24*(i-1) + 23
#     Wrotmat[:,i] = pcopy[:,idx]
# end

for i = 1:nt
    alpharotmat[i, :] = cchistory[i, :].alpha
    Wrotmat[i,:] = cchistory[i, :].W
end


tix = 50

### Check the inflow angle. 
phitplt = plot(rvec, cchistory[tix, :].phi.*(180/pi), lab="CCBlade", markershape=:x, xaxis="Radius (m)", yaxis="Inflow angle (deg)", title="Time Step $tix")
plot!(rof, phimat[tix,:], lab="OpenFAST")
# display(phitplt)

#So interestingly, the blade Azimuth is recorded as changing in OpenFAST. So I can't compare against that. Well we aren't currently testing azimuthal changes, so it shouldn't be a problem. 



### Check the angle of attack. 
# alphatplt = plot(rvec, cchistory[tix].alpha.*(180/pi), lab="CCBlade", markershape=:x, xaxis="Radius (m)", yaxis="Angle of Attack (deg)", title="Time Step $tix")
# plot!(rof, alphamat[tix,:], lab="OpenFAST")
# plot!(rvec, alpharotmat[tix,:].*(180/pi), lab="Rotors")
# display(alphatplt)

W0 = @. sqrt(vinf^2 + (rvec*omega)^2)


velocitytplt = plot(rvec, cchistory[tix, :].W, lab="CCBlade", markershape=:x, xaxis="Radius (m)", yaxis="Inflow Velocity (m/s)", title="Time Step $tix", leg=:bottomright)
# plot!(rof, Wmat[tix,:], lab="OpenFAST", markershape=:circle, markeralpha=0.2) #I think this is the rotational and freestream, affected by windshear, with no induced velocity. 
# plot!(rof, Wdmat[tix,:], lab="OF dis", markershape=:utriangle, markeralpha=0.2) # I think this is the inflow affected by windshear, but no induced velocity. 
plot!(rof, Vrelmat[tix,:], lab="OpenFAST", markershape=:rect, markeralpha=0.2)
plot!(rvec, Wrotmat[tix,:], lab="Rotors", markershape=:cross)
# plot!(rvec, W0, lab="W0")
# display(velocitytplt)






state = 15 #1:22
index = 1 #1:14 #States 1-3 are experience separation. 

states = Array{Float64}(undef, nt, na)

for i = 1:na
    idx = 22*(i-1) + state
    states[:, i] = xds[:,idx]
end



sp1 = plot(tvec, states[:,index], lab="Rotors", title="State $state, Index $index", xaxis="Time (s)", markershape=:x)
# display(sp1)



fxblade = cchistory[tix, :].Np.*cos(precone)
fyblade = cchistory[tix, :].Tp.*cos(precone)


loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load", leg=:topleft, title="Time $tix")
plot!(rvec, loads.Fx[tix,:], lab="Fx - Unsteady") #, markershape=:circle)
plot!(rvec, loads.Fy[tix,:], lab="Fy - Unsteady") #, markershape=:circle)
plot!(rof, fxmat[tix,:], markershape=:x, lab="Fx - OpenFAST")
plot!(rof, fymat[tix,:], markershape=:x, lab="Fy - OpenFAST")
plot!(rvec, fxblade, lab="Fx - CCBlade", markershape=:cross)
plot!(rvec, fyblade, lab="Fy - CCBlade", markershape=:cross)
# display(loadplt)








errfun(x, xt) = 100*(x - xt)/xt
absavg(xvec) = sum(abs.(xvec))/length(xvec)


# tix = 30
rof = adblade["BlSpn"] .+ rhub
fxblade = cchistory[tix, :].Np.*cos(precone)
fyblade = cchistory[tix, :].Tp

fxerr = zero(tvec)
fyerr = zero(tvec)

for i in eachindex(tvec)
    fxerr[i] = absavg(errfun.(loads.Fx[i,:], fxmat[i,:]))
    fyerr[i] = absavg(errfun.(loads.Fy[i,:], -fymat[i,:]))
end


loadplt = plot(xaxis="Radius (m)", yaxis="Distributed Load", leg=:topleft, title="Time Step $tix")
plot!(rvec, loads.Fx[tix,:], lab="Fx - Unsteady") #, markershape=:circle)
plot!(rvec, loads.Fy[tix,:], lab="Fy - Unsteady") #, markershape=:circle)
plot!(rof, fxmat[tix,:], markershape=:x, lab="Fx - OpenFAST")
plot!(rof, fymat[tix,:], markershape=:x, lab="Fy - OpenFAST")
plot!(rvec, fxblade, lab="Fx - CCBlade", markershape=:cross)
plot!(rvec, fyblade, lab="Fy - CCBlade", markershape=:cross)
# display(loadplt)


# coefplt = plot(xaxis="Radius (m)", yaxis="Load Coefficient", leg=:topright, title="Time Step $tix")
# plot!(rvec, coefs.Cn[tix,:], lab=L"$C_n$ - Rotors") #, markershape=:circle)
# plot!(rvec, coefs.Ct[tix,:], lab=L"$C_t$ - Rotors") #, markershape=:circle) #Todo: Need a way to get Ct out. 
# plot!(rof, cnmat[tix,:], markershape=:x, lab=L"$C_n$ - OpenFAST")
# plot!(rof, ccmat[tix,:], markershape=:x, lab=L"$C_t$ - OpenFAST")
# display(coefplt)
# savefig(coefplt, "radialloads_timestep$tix.png")

loadingerrplt = plot(xaxis="Time (s)", yaxis="Relative Error (%)", title="Blade average Loading Error", ylims=(0,0.3))
plot!(tvec, fxerr, lab="Fx")
plot!(tvec, fyerr, lab="Fy")
display(loadingerrplt)

#Note: The error has increased from the coefficients here because there is a slight error on phi between CCBlade and OpenFAST (even though they use the same method, they have slightly different tolerances). This slight error is incorporated here through the inflow velocity (which is used to dimensionalize the loadings) in a exponential relationship which increases the error. 












####### Calculate the error fo the thrust and torque. 

T, Q = rotorloads(rhub, rtip, rvec, loads, loads2, loads3)

Tcc = zero(tvec)
Qcc = zero(tvec)

for i = 1:nt
    Tcc[i], Qcc[i] = thrusttorque(rotor, sections, cchistory[i, :])
end

Terr = zeros(nt-1) #max error is 0.121% 
Qerr = zeros(nt-1) #Max error is 0.515%. 



for i = 1:nt-1
    Terr[i] = errfun(T[i], outs["RtFldFxg"][i])
    Qerr[i] = errfun(Q[i], -outs["RtFldMxg"][i])
end



thrustplt = plot(xaxis="Time (s)", yaxis="Thrust (N)")
plot!(tvec, T, lab="Rotors")
plot!(tvec, Tcc, lab="CCBlade")
plot!(tvec_of, outs["RtFldFxg"], lab="OpenFAST")
display(thrustplt)  
# savefig(thrustplt, "thrust_loosely_BLADG_010623.png") 


torqueplt = plot(xaxis="Time (s)", yaxis=L"Torque $\left(N\cdot m\right)$")
plot!(tvec, -Q, lab="Rotors")
plot!(tvec, Qcc, lab="CCBlade")
plot!(tvec_of, outs["RtFldMxg"], lab="OpenFAST")
display(torqueplt) 
# savefig(torqueplt, "torque_loosely_BLADG_010623.png")






nothing