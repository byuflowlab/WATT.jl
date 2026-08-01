#=
Shared setup for the GPU dynamic-stall prototype scripts.

Both `gpu_dsmodel_reference.jl` (CPU golden-trace generator) and
`gpu_dsmodel_validate.jl` (GPU-vs-reference check) include this file so they
build the *identical* NREL 5MW blade and the *identical* batch of prescribed
(U, aoa) time histories. The reference script marches the CPU Beddoes-Leishman
v3 (AeroDyn + Gonzalez) model over that batch and serializes the result; the
validate script rebuilds the same blade, loads the stored inputs, and runs the
GPU kernels against them.

The dynamic-stall model is driven by *prescribed* section-wise (U, aoa)
histories rather than a coupled BEM solve. This decouples the DS validation
from BEMT and lets us dial in scenarios (attached -> deep stall) that
deliberately exercise the model's data-dependent branches (sigma1/sigma3, LE/TE
separation, vortex shedding) for the warp-divergence study.

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels
const _of = OpenFASTTools
const _DS = DynamicStallModels

const OFPATH = joinpath(@__DIR__, "..", "data", "openfast")
const REF_FILE = joinpath(@__DIR__, "data", "gpu_dsmodel_reference.jld2")

"""
    build_nrel5mw_blade() -> (blade, rotor)

Build the NREL 5MW `Blade` (per-section DS airfoils, chord, twist, xcp) and a
turbine `Rotor` from the OpenFAST inputs in `data/openfast/`. Mirrors the setup
in `aero_only_loads.jl` / `gpu_bemt_benchmark.jl` so the DS prototype runs on
the exact same geometry.
"""
function build_nrel5mw_blade()
    adblade = _of.read_adblade("sn5_adblade.dat", OFPATH)
    edfile  = _of.read_edfile("sn5_EDfile.dat", OFPATH)

    aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                    "DU30_A17.dat", "DU25_A17.dat", "DU21_A17.dat", "NACA64_A17.dat")
    aftypes = [_of.read_airfoilinput(joinpath(OFPATH, "airfoils", name)) for name in aftype_names]
    af_idx  = Int.(adblade["BlAFID"])
    afs     = aftypes[af_idx]

    chordvec = adblade["BlChord"]
    twistvec = adblade["BlTwist"] .* (pi / 180)
    rhub = edfile["HubRad"]
    rvec = adblade["BlSpn"] .+ rhub
    rtip = rvec[end]
    n    = length(rvec)

    airfoils = Vector{_DS.Airfoil}(undef, n)
    xcp = Vector{Float64}(undef, n)
    for i = 1:n
        airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
    end

    blade = WATT.Blade(rvec, chordvec, twistvec, xcp, airfoils; rhub, rtip)
    rotor = WATT.Rotor(3, 80.0, true)
    return blade, rotor
end

# Prescribed operating scenarios spanning attached flow through deep stall.
# Each scenario is one "sim" in the batch; angles in degrees, freq in Hz,
# Uscale multiplies the reference relative-wind magnitude.
const DS_SCENARIOS = (
    (name = "attached_low",  aoa_mean =  3.0, aoa_amp =  1.0, freq = 0.5, Uscale = 0.6),
    (name = "attached_mid",  aoa_mean =  6.0, aoa_amp =  2.0, freq = 1.0, Uscale = 0.8),
    (name = "light_stall",   aoa_mean = 10.0, aoa_amp =  4.0, freq = 1.5, Uscale = 1.0),
    (name = "moderate",      aoa_mean = 13.0, aoa_amp =  6.0, freq = 2.0, Uscale = 0.9),
    (name = "deep_stall",    aoa_mean = 16.0, aoa_amp = 10.0, freq = 2.5, Uscale = 1.0),
    (name = "deep_stall_hf", aoa_mean = 17.0, aoa_amp = 12.0, freq = 3.0, Uscale = 0.7),
)

"""
    make_ds_batch(blade; scenarios=DS_SCENARIOS, nt=250, dt=0.01, Uref=60.0)
        -> (U, aoa, tvec, dt)

Build the prescribed input histories for the DS batch.

`U` and `aoa` are `(n_sections, n_sims, nt)` arrays. Within a sim, `aoa` follows
`aoa_mean + aoa_amp*sin(2π·freq·t)` (constant across span); `U` is constant in
time but varies across span as `Uref·Uscale·(0.35 + 0.65·(r-rhub)/(rtip-rhub))`
so the batch spans a realistic root->tip range of Mach and DS time constants.
Only `U` and `aoa` are needed on device (Udot/alphadot are unused by the ADG
model).
"""
function make_ds_batch(blade; scenarios = DS_SCENARIOS, nt::Int = 250, dt::Float64 = 0.01, Uref::Float64 = 60.0)
    r    = blade.r
    rhub = blade.rhub
    rtip = blade.rtip
    n_sections = length(r)
    n_sims = length(scenarios)

    tvec = collect(0:dt:(dt * (nt - 1)))
    U   = Array{Float64}(undef, n_sections, n_sims, nt)
    aoa = Array{Float64}(undef, n_sections, n_sims, nt)

    for (s, sc) in enumerate(scenarios)
        span_frac = @. 0.35 + 0.65 * (r - rhub) / (rtip - rhub)
        Uj = @. Uref * sc.Uscale * span_frac                 # (n_sections,)
        amean = sc.aoa_mean * pi / 180
        aamp  = sc.aoa_amp  * pi / 180
        for i in 1:nt
            a_i = amean + aamp * sin(2pi * sc.freq * tvec[i])
            for j in 1:n_sections
                U[j, s, i]   = Uj[j]
                aoa[j, s, i] = a_i
            end
        end
    end
    return U, aoa, tvec, dt
end
