#=
Simplify the NREL 5MW to compare OpenFAST, GXBeam, and eventually Rotors.jl. 


Adam Cardoza 01/23/23

=#


using OpenFASTTools, Rotors, DynamicStallModels, FLOWMath

of = OpenFASTTools

path = dirname(@__FILE__)
cd(path)

### Read in OpenFAST files
ofpath = "../OpenFAST_NREL5MW/" 
inputfile = of.read_inputfile("NREL5MW_input.fst", ofpath)
inflowwind = of.read_inflowwind("NREL5MW_inflowwind.dat", ofpath)
# addriver = of.read_addriver("NREL5MW_ADdriver.inp", ofpath)
adfile = of.read_adfile("NREL5MW_ADfile.dat", ofpath)
adblade = of.read_adblade("NREL5MW_ADblade.dat", ofpath)
edfile = of.read_edfile("NREL5MW_edfile.dat", ofpath)
bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
bdblade = of.read_bdblade("NREL5MW_bdblade.dat", ofpath)
bddriver = of.read_bddriver("NREL5MW_bddriver.inp", ofpath)

### Simulation control
tmax = 3.0
dt = 0.001
# dt_out = 0.01 # "\"default\""
dt_out = dt

### Environmental variables
density = 1.225 #density (kg/m^3)
kinvisc = 1.5e-5 #Kinematic Viscosity (m^2/s)
a = 335 #Speed of Sound (m/s)
gravity = 9.81 #Gravity (m/s^2)

Uinf = 10 #Freestream velocity (m/s)
patm = 101325 #Atmospheric pressure (Pa)
rpm = 11.44 #angular velocity (rotations per minute)
omega = rpm*2*pi/60 #angular velocity (radians/second)

refht = 90 #Reference height for the power law exponent shear rule
plexp = 0.01


### Turbine description

pitch = 0.0 
tilt = 0.0
precone = 0.0

B = 3 #Number of blades
hubht = 90 #Hub height (m)
overhang = -1 

rhub = 1.5 #Hub radius
rtip = 63.0 #Tip radius

# rvec = adblade["BlSpn"] .+ rhub

# n = 200
# rvec = collect(range(rhub, rtip, length=n))
# rfrac = (rvec .- rhub)/(rtip-rhub)

nelem = 1
np_elem = 25

# chordvec = ones(n)
# twistvec = zeros(n)
# # twistvec = -ones(n).*3.0  
# AFID = ones(n).*8

E = 1.1e10 #6.83e10 #Young's Modulus
nu = 0.3 #Poisson's Ratio
m = 10.0 #2700.0 #kg/m^3 -> Decreasing the weight decreases the frequency of oscillation
h = 0.25 # Thickness (meters)
w = 0.25 # Width (meters)
mu = 0.001 #Damping ratio

A = w*h #Area

Ix = w*(h^3)/12  #Area moment of inertia
Iy = h*(w^3)/12
J = Ix + Iy

G=E/(2*(1+nu))

# K = of.stiffness_matrix(E, nu, h*w, Ix, Iy) #Stiffness matrix

#[3, 1, 2, 6, 4, 5] #BD to GX indexing
# gxtobdidx = [2, 3, 1, 5, 6, 4] #GX to BD indexing
# K1 = [(E*A) 0 0 0 0 0
#     0 (G*A) 0 0 0 0
#     0 0 (G*A) 0 0 0
#     0 0 0 (G*J) 0 0
#     0 0 0 0 (E*Iy) 0
#     0 0 0 0 0 (E*Ix)] #GXBeam stiffness matrix
# K = K[gxtobdidx, gxtobdidx] #Convert to BD stiffness matrix
#Note: both K and K1 produce the same matrix. :| Which means that I did it right the first go aroound. Which means that BeamDyn just can't solve my stiffness matrix. 

# l = rvec[2]-rvec[1] # Element length

# dm = m*w*h #Distributed mass of the section #TODO: Dr. Ning suggested setting this to zero to simplify terms.  
# Mix = dm*Ix #Mass moment of inertia
# Miy = dm*Iy
# M = of.mass_matrix(dm, Mix, Miy)




inflowfile = 0

let
    ### Set OpenFAST file values
    ##### ADdriver
    # addriver["TMax"] = tmax
    # addriver["DT"] = dt
    # addriver["FldDens"] = density
    # addriver["KinVisc"] = kinvisc
    # addriver["SpdSound"] = a
    # addriver["CompInflow"] = inflowfile
    # addriver["HWindSpeed"] = Uinf
    # addriver["RefHt"] = refht
    # addriver["PLExp"] = plexp
    # addriver["NumBlades(1)"] = B
    # addriver["HubRad(1)"] = rhub
    # addriver["HubHt(1)"] = hubht
    # addriver["Overhang(1)"] = overhang
    # addriver["ShftTilt(1)"] = tilt
    # addriver["PreCone(1)"] = precone
    # addriver["HubHt"] = hubht
    # addriver["NumCases"] = 1
    # addriver["HWndSpeed_mat"] = [Uinf]
    # addriver["PLExp_mat"] = [plexp]
    # addriver["RotSpd_mat"] = hubht
    # addriver["Pitch_mat"] = [pitch]
    # addriver["dT_mat"] = dt
    # addriver["Tmax_mat"] = tmax

    
    bdfile["member_total"] = nelem
    members = zeros(nelem, 2)
    for i = 1:nelem
        members[i,1] = i
        members[i,2] = np_elem
    end
    bdfile["KeyPairs"] = members

    n = Int(sum(members[:,2]) - nelem + 1)
    rvec = collect(range(rhub, rtip, length=n))
    rfrac = (rvec .- rhub)/(rtip-rhub)
    rnew = rvec.-rhub
    # chordvec = ones(n)
    # twistvec = zeros(n)
    # twistvec = -ones(n).*3.0  
    # AFID = ones(n).*8
    chordfit = Akima(adblade["BlSpn"], adblade["BlChord"])
    twistfit = Akima(adblade["BlSpn"], adblade["BlTwist"])

    chordvec = chordfit.(rnew)
    twistvec = twistfit.(rnew)

    AFID = of.integerfit(adblade["BlSpn"], adblade["BlAFID"], rnew)
    for i = 1:n
        if AFID[i]<3
            AFID[i]=3
        end
    end



    ##### AD Primary file
    adfile["IndToler"] = 2e-12
    adfile["MaxIter"] = 10000
    adfile["WakeMod"] = 1
    adfile["AFAeroMod"] = 2
    # adfile["AFNames"] = ["./Airfoils/NACA64_A17.dat"]
    # adfile["NumAFfiles"] = length(adfile["AFNames"])
    adfile["ADBlFile(1)"] = "sn5_ADblade.dat"
    adfile["ADBlFile(2)"] = "sn5_ADblade.dat"
    adfile["ADBlFile(3)"] = "sn5_ADblade.dat"
    adfile["OutList"] = ["B1Azimuth"]
    adfile["NodeOutList"] = ["Alpha", "Vx", "Vy", "Fx", "Fy", "Mm", "AxInd", "TnInd", "Phi", "Theta"]
    # append!(adfile["NodeOutList"], ["Alpha"])



    #### AD Blade file
    adblade["NumBlNds"] = n
    # ltemp = rtip - rhub
    # rtemp = ltemp.*adblade["BlSpn"]./61.5
    # adblade["BlSpn"] = rtemp
    adblade["BlSpn"] = rvec .- rhub
    adblade["BlSpn"][end] -= 0.0001
    adblade["BlCrvAC"] = zeros(n)
    adblade["BlSwpAC"] = zeros(n)
    adblade["BlCrvAng"] = zeros(n)

    # adblade["BlTwist"] = zero(adblade["BlTwist"])
    # adblade["BlChord"] = ones(length(adblade["BlChord"]))
    # adblade["BlAFID"] = ones(length(adblade["BlAFID"])).*8

    adblade["BlTwist"] = twistvec
    adblade["BlChord"] = chordvec
    adblade["BlAFID"] = AFID



    ### ED Primary
    edfile["TeetDOF"] = false ## Rotor-Teeter DOF (flag)
    edfile["DrTrDOF"] = false ## Drivetrain rotation flexibility DOF (flag)
    edfile["GenDOF"]    = false ## Generator DOF (flag)
    edfile["YawDOF"]    = false ## Yaw DOF (flag)
    edfile["TwFADOF1"]  = false ## First fore-aft tower bending-mode DOF (flag)
    edfile["TwFADOF2"]  = false ## Second fore-aft tower bending-mode DOF (flag)
    edfile["TwSSDOF1"]  = false ## First side-to-side tower bending-mode DOF (flag)
    edfile["TwSSDOF2"]  = false ## Second side-to-side tower bending-mode DOF (flag)
    edfile["PtfmSgDOF"] = false ## Platform horizontal surge translation DOF (flag)
    edfile["PtfmSwDOF"] = false ## Platform horizontal sway translation DOF (flag)
    edfile["PtfmHvDOF"] = false ## Platform vertical heave translation DOF (flag)
    edfile["PtfmRDOF"]  = false ## Platform roll tilt rotation DOF (flag)
    edfile["PtfmPDOF"]  = false ## Platform pitch tilt rotation DOF (flag)
    edfile["PtfmYDOF"]  = false ## Platform yaw rotation DOF (flag)

    # #Initial Conditions
    # edfile["BlPitch(1)"] = pitch
    # edfile["BlPitch(2)"] = pitch
    # edfile["BlPitch(3)"] = pitch
    edfile["RotSpeed"] = rpm 

    # #Turbine Config
    # edfile["NumBl"] = B
    edfile["TipRad"] = rtip
    edfile["HubRad"] = rhub
    edfile["PreCone(1)"] = precone
    edfile["PreCone(2)"] = precone
    edfile["PreCone(3)"] = precone
    # edfile["OverHang"] = overhang
    edfile["ShftTilt"] = tilt

    #Blade
    edfile["BldNodes"] = n
    edfile["BldFile(1)"] = "sn5_EDblade.dat"
    edfile["BldFile(2)"] = "sn5_EDblade.dat"
    edfile["BldFile(3)"] = "sn5_EDblade.dat"

    edfile["TwrFile"] = "sn5_EDtower.dat"
    edfile["OutList"] = ["Azimuth"]


    ### BD Driver
    bddriver["t_final"] = tmax
    bddriver["dt"] = dt


    ### BD Primary #Todo: Add stuff to automatically update the BDDriver as well. 
    bdfile["QuasiStaticInit"] = true
    bdfile["NRMax"] = 10000
    bdfile["quadrature"] = 2
    bdfile["DTBeam"] = "DEFAULT"
    bdfile["order_elem"] = 8
    bdfile["stop_tol"] = 1e-8
    bdfile["rhoinf"] = 1.0

    # nelem = 1
    # bdfile["member_total"] = nelem
    # np_elem = Int(n/nelem)
    # members = zeros(nelem, 2)
    # for i = 1:nelem
    #     members[i,1] = i
    #     if i==1
    #         members[i,2] = np_elem
    #     else
    #         members[i,2] = np_elem + 1
    #     end
    # end
    # bdfile["KeyPairs"] = members

    bdfile["kp_total"] = n
    # bdfile["KeyPairs"] = [1 n]
    bdfile["kp_xr"] = zeros(n)
    bdfile["kp_yr"] = zeros(n)
    
    # bdfile["kp_zr"] = ltemp.*bdfile["kp_zr"]./61.5
    # bdfile["initial_twist"] = zero(bdfile["initial_twist"])

    bdfile["kp_zr"] = rvec .- rhub
    bdfile["initial_twist"] = twistvec

    

    bdfile["BldFile"] = ["sn5_BDblade.dat"]
    bdfile["BldNd_BlOutNd"] = 99

    bdfile["OutNd"] = Int[1]
    bdfile["OutList"] = [" ", "TipTDxr", "TipTDyr", "TipTDzr", "TipRDxr", "TipRDyr", "TipRDzr"]
    # bdfile["NodeOutList"] = ["TDxr"]
    append!(bdfile["NodeOutList"], ["RDxr", "RDyr", "RDzr", "TVxr", "TVyr", "TVzr", "RVxr", "RVyr", "RVzr", "DFxr", "DFyr", "DFzr" ])

    # ### BD Blade file
    Kmat = zeros(6,6,Int(bdblade["station_total"]))
    Mmat = zeros(6,6,Int(bdblade["station_total"]))
    
    for i = 1:Int(bdblade["station_total"])
        Kmat[:,:,i] = bdblade["K$i"]
        Mmat[:,:,i] = bdblade["M$i"]
    end

    Kfit = Rotors.interpolate_matrix_symmetric(bdblade["rfrac"], rfrac, Kmat; fit=Linear)
    Mfit = Rotors.interpolate_matrix_symmetric(bdblade["rfrac"], rfrac, Mmat; fit=Linear)

    for i = 1:n
        bdblade["K$i"] = Kfit[:,:,i]
        bdblade["M$i"] = Mfit[:,:,i]
    end

    
    bdblade["station_total"] = n
    bdblade["rfrac"] = rfrac
    bdblade["mu"] .= mu


    bddriver["InputFile"] = "sn5_BDfile.dat"
    bddriver["DynamicSolve"] = true
    bddriver["t_final"] = tmax
    bddriver["dt"] = dt

    ### Inflowwind
    inflowwind["WindType"] = 1
    inflowwind["NWindVel"] = 0
    inflowwind["HWindSpeed"] = Uinf
    inflowwind["RefHt"] = hubht
    inflowwind["PLexp"] = plexp




    ### Input file
    inputfile["CompElast"] = 2
    inputfile["CompInflow"] = 1
    inputfile["CompAero"] = 2

    inputfile["EDFile"] = "sn5_EDfile.dat"
    inputfile["AeroFile"] = "sn5_ADfile.dat"
    inputfile["InflowFile"] = "sn5_inflowwind.dat"
    inputfile["BDBldFile(1)"] = "sn5_BDfile.dat"
    inputfile["BDBldFile(2)"] = "sn5_BDfile.dat"
    inputfile["BDBldFile(3)"] = "sn5_BDfile.dat"


    inputfile["TMax"] = tmax
    inputfile["DT"] = dt
    # inputfile["DT_Out"] = dt
    inputfile["DT_Out"] = dt_out
    inputfile["NumCrctn"] = 0

    inputfile["Gravity"] = gravity
    inputfile["AirDens"] = density
    inputfile["KinVisc"] = kinvisc
    inputfile["SpdSound"] = a
    inputfile["Patm"] = patm

    



    ### Write OpenFAST files
    # of.write_addriver(addriver, "sn5_ADdriver.dvr"; outputpath="./")
    of.write_adfile(adfile, "sn5_ADfile.dat"; outputpath="./")
    of.write_adblade(adblade, "sn5_ADblade.dat"; outputpath="./")
    of.write_edfile(edfile, "sn5_EDfile.dat"; outputpath="./")
    of.write_bdfile(bdfile, "sn5_BDfile.dat"; outputpath="./")
    of.write_bdblade(bdblade, "sn5_BDblade.dat"; outputpath="./")
    of.write_bddriver(bddriver, "sn5_bddriver.inp"; outputpath="./")
    of.write_inflowwind(inflowwind, "sn5_inflowwind.dat"; outputpath="./")
    of.write_inputfile(inputfile, "sn5_input.fst"; outputpath="./")
end





nothing