#=
A simple turbine to compare OpenFAST, GXBeam, and eventually Rotors.jl. 


Adam Cardoza 01/23/23

=#


using OpenFASTsr

of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

### Simulation control
tmax = 2.0
dt = 0.0002
dt_out = "\"default\""

### Environmental variables
density = 1.225 #density (kg/m^3)
kinvisc = 1.5 #Kinematic Viscosity (m^2/s)
a = 335 #Speed of Sound (m/s)
gravity = 9.81 #Gravity (m/s^2)

Uinf = 10 #Freestream velocity (m/s)
patm = 101325 #Atmospheric pressure (Pa)
rpm = 50 #angular velocity (rotations per minute)
omega = rpm*2*pi/60 #angular velocity (radians/second)

refht = 20 #Reference height for the power law exponent shear rule
plexp = 0.0


### Turbine description

pitch = 0.0 
tilt = 0.0
precone = 0.0

B = 3 #Number of blades
hubht = 20 #Hub height (m)
overhang = -1 

rhub = 1.0 #Hub radius
rtip = 10.0 #Tip radius

rvec = collect(rhub:.1:rtip)

n = length(rvec)

chordvec = 0.25*ones(n)
twistvec = zeros(n)
AFID = ones(n)

E = 1.1e9 #6.83e10 #Young's Modulus
nu = 0.3 #Poisson's Ratio
m = 500 #2700.0 #kg/m^3 -> Decreasing the weight decreases the frequency of oscillation
h = 1.5 # Thickness (meters)
w = 1.5 # Width (meters)
mu = 0.01 #Damping ratio



Ix = w*(h^3)/12  #Area moment of inertia
Iy = h*(w^3)/12
K = of.stiffness_matrix(E, nu, h*w, Ix, Iy) #Stiffness matrix

l = rvec[2]-rvec[1] # Element length

dm = m*l #Distributed mass of the section 
Mix = dm*Ix #Mass moment of inertia
Miy = dm*Iy
M = of.mass_matrix(dm, Mix, Miy)


### Read in OpenFAST files
ofpath = "./simpleturbine/" 
inputfile = of.read_inputfile("simple_input.fst", ofpath)
inflowwind = of.read_inflowwind("simple_inflowwind.dat", ofpath)
addriver = of.read_addriver("simple_ADdriver.dvr", ofpath)
adfile = of.read_adfile("simple_ADfile.dat", ofpath)
adblade = of.read_adblade("simple_ADblade.dat", ofpath)
edfile = of.read_edfile("simple_edfile.dat", ofpath)
bdfile = of.read_bdfile("simple_bdfile.dat", ofpath)
bdblade = of.read_bdblade("simple_bdblade.dat", ofpath)

inflowfile = 0

let
    ### Set OpenFAST file values
    ##### ADdriver
    addriver["TMax"] = tmax
    addriver["DT"] = dt
    addriver["FldDens"] = density
    addriver["KinVisc"] = kinvisc
    addriver["SpdSound"] = a
    addriver["CompInflow"] = inflowfile
    addriver["HWindSpeed"] = Uinf
    addriver["RefHt"] = refht
    addriver["PLExp"] = plexp
    addriver["NumBlades(1)"] = B
    addriver["HubRad(1)"] = rhub
    addriver["HubHt(1)"] = hubht
    addriver["Overhang(1)"] = overhang
    addriver["ShftTilt(1)"] = tilt
    addriver["PreCone(1)"] = precone
    addriver["HubHt"] = hubht
    addriver["NumCases"] = 1
    addriver["HWndSpeed_mat"] = [Uinf]
    addriver["PLExp_mat"] = [plexp]
    addriver["RotSpd_mat"] = hubht
    addriver["Pitch_mat"] = [pitch]
    addriver["dT_mat"] = dt
    addriver["Tmax_mat"] = tmax



    ##### AD Primary file
    adfile["AFNames"] = ["./Airfoils/NACA64_A17.dat"]
    adfile["NumAFfiles"] = length(adfile["AFNames"])
    adfile["ADBlFile(1)"] = "simple_ADblade.dat"
    adfile["ADBlFile(2)"] = "simple_ADblade.dat"
    adfile["ADBlFile(3)"] = "simple_ADblade.dat"



    #### AD Blade file
    adblade["NumBlNds"] = n
    adblade["BlSpn"] = rvec .- rhub
    adblade["BlCrvAC"] = zeros(n)
    adblade["BlSwpAC"] = zeros(n)
    adblade["BlCrvAng"] = zeros(n)
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

    #Initial Conditions
    edfile["BlPitch(1)"] = pitch
    edfile["BlPitch(2)"] = pitch
    edfile["BlPitch(3)"] = pitch
    edfile["RotSpeed"] = rpm 

    #Turbine Config
    edfile["NumBl"] = B
    edfile["TipRad"] = rtip
    edfile["HubRad"] = rhub
    edfile["PreCone(1)"] = precone
    edfile["PreCone(2)"] = precone
    edfile["PreCone(3)"] = precone
    edfile["OverHang"] = overhang
    edfile["ShftTilt"] = tilt

    #Blade
    edfile["BldNodes"] = n
    edfile["BldFile(1)"] = "simple_BDblade.dat"
    edfile["BldFile(2)"] = "simple_BDblade.dat"
    edfile["BldFile(3)"] = "simple_BDblade.dat"

    edfile["TwrFile"] = "simple_EDtower.dat"




    ### BD Primary
    bdfile["member_total"] = 1
    bdfile["kp_total"] = n
    bdfile["KeyPairs"] = [1 n]
    bdfile["kp_xr"] = zeros(n)
    bdfile["kp_yr"] = zeros(n)
    bdfile["kp_zr"] = rvec .- rhub
    bdfile["initial_twist"] = twistvec

    bdfile["BldFile"] = ["simple_BDblade.dat"]


    ### BD Blade file
    bdblade["station_total"] = n
    bdblade["rfrac"] = (rvec .- rhub)/(rtip-rhub)
    bdblade["mu"] .= mu
    
    for i = 1:n
        bdblade["K$i"] = K
        bdblade["M$i"] = M
    end


    ### Inflowwind
    inflowwind["WindType"] = 1
    inflowwind["NWindVel"] = 0
    inflowwind["HWindSpeed"] = Uinf
    inflowwind["RefHt"] = hubht
    inflowwind["PLExp"] = plexp



    ### Input file
    inputfile["CompElast"] = 2
    inputfile["CompInflow"] = 1
    inputfile["CompAero"] = 2

    inputfile["TMax"] = tmax
    inputfile["DT"] = dt

    inputfile["Gravity"] = gravity
    inputfile["AirDens"] = density
    inputfile["KinVisc"] = kinvisc
    inputfile["SpdSound"] = a
    inputfile["Patm"] = patm

    inputfile["DT_Out"] = dt_out



    ### Write OpenFAST files
    of.write_addriver(addriver, "simple_ADdriver.dvr"; outputpath="./simpleturbine")
    of.write_adfile(adfile, "simple_ADfile.dat"; outputpath="./simpleturbine")
    of.write_adblade(adblade, "simple_ADblade.dat"; outputpath="./simpleturbine")
    of.write_edfile(edfile, "simple_EDfile.dat"; outputpath="./simpleturbine")
    of.write_bdfile(bdfile, "simple_BDfile.dat"; outputpath="./simpleturbine")
    of.write_bdblade(bdblade, "simple_BDblade.dat"; outputpath="./simpleturbine")
    of.write_inflowwind(inflowwind, "simple_inflowwind.dat"; outputpath="./simpleturbine")
    of.write_inputfile(inputfile, "simple_input.fst"; outputpath="./simpleturbine")
end





nothing