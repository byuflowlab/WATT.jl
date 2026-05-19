#=
Shared NREL 5MW test fixture.

Builds a single `blade`, `rotor`, `env_rated`, and `assembly` from the OpenFAST
input files in `../data/openfast`, matching the setup pattern in
test_bem.jl / test_environments.jl / test_gxbeam.jl. New Phase 2 test files
`include` this file to avoid duplicating the ~60-line setup boilerplate.

Pre-existing test files (test_bem.jl, test_environments.jl, ...) keep their
inline setup — they double as transform regressions and shouldn't be coupled
to this fixture.

Exposes (at the include-site scope):
    aftypes, afs, airfoils, xcp, chordvec, twistvec, rvec, rhub, rtip, n,
    blade, rotor, env_rated, omega_rated, assembly, ne

Adam Cardoza
=#

using WATT, FLOWMath, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, StructArrays

const _of = OpenFASTTools
const _DS = DynamicStallModels

let ofpath = joinpath(@__DIR__, "..", "..", "data", "openfast")
    global aftypes, afs, airfoils, xcp, chordvec, twistvec, rvec, rhub, rtip, n
    global blade, rotor, env_rated, omega_rated, assembly, ne

    adblade = _of.read_adblade("sn5_adblade.dat", ofpath)
    edfile  = _of.read_edfile("sn5_EDfile.dat", ofpath)
    bdfile  = _of.read_bdfile("sn5_BDfile.dat", ofpath)
    bdblade = _of.read_bdblade("sn5_BDblade.dat", ofpath)

    aftypes = Array{_of.AirfoilInput}(undef, 8)
    aftypes[1] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder1.dat"))
    aftypes[2] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder2.dat"))
    aftypes[3] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU40_A17.dat"))
    aftypes[4] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU35_A17.dat"))
    aftypes[5] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU30_A17.dat"))
    aftypes[6] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU25_A17.dat"))
    aftypes[7] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU21_A17.dat"))
    aftypes[8] = _of.read_airfoilinput(joinpath(ofpath, "airfoils", "NACA64_A17.dat"))

    af_idx = Int.(adblade["BlAFID"])
    afs    = aftypes[af_idx]

    chordvec = adblade["BlChord"]
    twistvec = adblade["BlTwist"]            # degrees (from OpenFAST)
    rhub = edfile["HubRad"]
    rvec = adblade["BlSpn"] .+ rhub
    rtip = rvec[end]
    n    = length(rvec)

    airfoils = StructArray{_DS.Airfoil}(undef, n)
    xcp = Vector{Float64}(undef, n)
    for i = 1:n
        airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
    end

    # Twist must be in radians for WATT.Blade / CCBlade.Section.
    blade = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; rhub, rtip)

    # NREL 5MW rated: 10 m/s, TSR ≈ 7.55.
    B        = 3
    hubht    = 80.0
    turbine  = true
    vinf     = 10.0
    tsr      = 7.55
    omega_rated = vinf*tsr/rtip
    rho      = 1.225
    mu       = 1.464e-5
    a        = 343.0
    shearexp = 0.0

    rotor     = WATT.Rotor(B, hubht, turbine)
    env_rated = WATT.environment(rho, mu, a, vinf, omega_rated, shearexp)

    assembly = _of.make_assembly(edfile, bdfile, bdblade)
    ne       = length(assembly.elements)
end

nothing
