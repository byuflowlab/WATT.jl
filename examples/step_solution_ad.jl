#=
Frozen-start windowed sensitivity through the coupled single-step primitive.

Demonstrates the downstream surrogate workflow: march the coupled solution rest→s in Float64
(the "warm-up"), snapshot the native state at s, then restart from that frozen snapshot and march a
short window s→s+k with ForwardDiff to get ∂(deflection)/∂(structural parameters) — with the
snapshot held at ZERO partials (a constant initial condition, not itself a function of the params).

The plumbing:
  • `initialize_from_state` injects the snapshot AssemblyState into a caller-supplied mesh WITHOUT a
    rest-IC solve (it reuses GXBeam's `initialize_system!(...; initial_state=…)` path).
  • `run_from_state!` marches `step_solution!` over the window.
  • The structural parameters `p` (per-element compliance + mass) are threaded through GXBeam via a
    `pfunc` (rebuilds the assembly from `p`) and a `prepp!` callback (writes aero loads into `p`) —
    the same WATT-native AD mechanism `run_sim!` uses. Passing a Float64 snapshot with a Dual `p`
    makes GXBeam promote the injected state to Dual-with-zero-partials automatically.

Companion to `examples/coupled_run_sim.jl` (the un-differentiated reference solve).

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, JLD2
using ForwardDiff, LinearAlgebra

const of = OpenFASTTools
const DS = DynamicStallModels

datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")

# --- Blade fixture + dynamic-stall airfoils (identical setup to coupled_run_sim.jl). ---
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)

aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat",  "DU25_A17.dat",  "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]
af_idx  = of.integerfit(raf, afidx, rvec)
afs     = aftypes[af_idx]

dsairfoils = Vector{DS.Airfoil}(undef, nr)
xcp        = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    if cylinder_mask[i]
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; dsmodel=DS.NoModel(), polar=polars[i])
    else
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
    end
end

blade    = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)
rotor    = WATT.Rotor(3, 90.0, true)
assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list, midpoints=xp, mass=mass_list)

rho, mu, a, shearexp = 1.225, 1.837e-5, 335.0, 0.2
omega    = 11.4 * 7.55 / (Rtip * cos(precone))
env      = environment(joinpath(ofpath, "TurbSim.dat"), rho, mu, a, omega, shearexp)

nelem  = length(assembly.elements)
n_aero = nr

# ------------------------------------------------------------------
# Structural parameter vector p = [vec(C_1)…vec(C_nelem); vec(M_1)…vec(M_nelem); Fx…; Fy…].
# pfunc rebuilds the assembly from p and maps aero loads (per node) to per-element loads;
# prepp! writes the current aero loads into the load slots of p. (Same pattern as test/test_ad.jl.)
# ------------------------------------------------------------------
const NCE = 36
n_comp_(ne)    = NCE * ne
n_mass_(ne)    = NCE * ne
n_loads_(na)   = 2 * na
n_p_(ne, na)   = n_comp_(ne) + n_mass_(ne) + n_loads_(na)

function pack_p(assembly, n_aero)
    ne = length(assembly.elements)
    p  = Vector{Float64}(undef, n_p_(ne, n_aero))
    for i in 1:ne
        p[(i-1)*NCE+1 : i*NCE] = vec(Matrix(assembly.elements[i].compliance))
    end
    off = n_comp_(ne)
    for i in 1:ne
        p[off + (i-1)*NCE+1 : off + i*NCE] = vec(Matrix(assembly.elements[i].mass))
    end
    p[off + n_mass_(ne) + 1 : end] .= 0.0
    return p
end

function make_pfunc(assembly, n_aero)
    ne     = length(assembly.elements)
    pts    = assembly.points
    st     = 1:ne
    sp     = 2:length(pts)
    xp     = [assembly.elements[i].x for i in 1:ne]
    nc     = n_comp_(ne)
    nm     = n_mass_(ne)
    function pfunc(p, t)
        comp = [SMatrix{6,6}(reshape(p[(i-1)*NCE+1 : i*NCE], 6, 6)) for i in 1:ne]
        ms   = [SMatrix{6,6}(reshape(p[nc + (i-1)*NCE+1 : nc + i*NCE], 6, 6)) for i in 1:ne]
        Fx   = @view p[nc + nm + 1          : nc + nm + n_aero]
        Fy   = @view p[nc + nm + n_aero + 1 : nc + nm + 2*n_aero]
        asm  = GXBeam.Assembly(pts, st, sp; compliance=comp, mass=ms, midpoints=xp)
        TF   = eltype(p)
        z    = SVector{3,TF}(0.0, 0.0, 0.0)
        dl = Dict(i => begin
            L = asm.elements[i].L
            GXBeam.DistributedLoads{TF}(
                SVector{3,TF}(0.0, Fy[i]*L/2, -Fx[i]*L/2),
                SVector{3,TF}(0.0, Fy[i]*L/2, -Fx[i]*L/2),
                z, z, z, z, z, z)
        end for i in 1:ne)
        return (; assembly=asm, distributed_loads=dl)
    end
    return pfunc
end

function make_prepp(nelem, n_aero)
    nc = n_comp_(nelem); nm = n_mass_(nelem)
    function prepp!(pp, Fx, Fy, mx)
        pp[nc + nm + 1          : nc + nm + n_aero]   .= Fx
        pp[nc + nm + n_aero + 1 : nc + nm + 2*n_aero] .= Fy
    end
    return prepp!
end

p0     = pack_p(assembly, n_aero)
pfunc  = make_pfunc(assembly, n_aero)
prepp! = make_prepp(nelem, n_aero)
a_skel = pfunc(p0, 0.0).assembly

# ------------------------------------------------------------------
# 1. Warm up rest→s in Float64 and freeze the native snapshot at s.
# ------------------------------------------------------------------
tvec = collect(0.0:0.05:1.5)
s    = 20
aw, gw, mw = initialize_sim(blade, a_skel, tvec; verbose=false, pfunc=pfunc, p=p0)
run_sim!(rotor, blade, mw, env, tvec, aw, gw; verbose=false, prepp=prepp!, p=p0)

state_s   = gw[s]
xds_s     = copy(aw.xds[s, :])
azimuth_s = aw.azimuth[s]
println("Warmed up to step s=$s (t=$(tvec[s]) s). Frozen tip deflection: ",
        round(state_s.elements[end].u[3]; digits=5), " m")

# ------------------------------------------------------------------
# 2. Frozen-start windowed map, differentiated w.r.t. one mid-span element's
#    compliance block (36 entries) for a fast, legible demo. The remaining
#    entries of p are held at their p0 values.
# ------------------------------------------------------------------
tvec_win = tvec[s:s+4]          # 4-step window
t_s      = tvec[s]
ielem    = nelem ÷ 2
c_idx    = (ielem-1)*NCE+1 : ielem*NCE
c0       = p0[c_idx]

function windowed_tip(cblock)
    TF = eltype(cblock)
    p  = convert(Vector{TF}, p0)
    p[c_idx] .= cblock
    bld = WATT.Blade(blade.r, convert(Vector{TF}, blade.c), convert(Vector{TF}, blade.twist),
                     blade.xcp, blade.airfoils; rhub=blade.rhub, rtip=blade.rtip, precone=blade.precone)
    _, _, meshw = initialize_sim(bld, a_skel, tvec_win; verbose=false, pfunc=pfunc, p=p)
    init = initialize_from_state(state_s, xds_s, azimuth_s, meshw, bld, env, tvec_win, t_s; p=p)
    out  = Vector{TF}(undef, length(tvec_win) - 1)
    run_from_state!(init..., meshw, rotor, bld, env, tvec_win;
                    prepp=prepp!, p=p, out=(st, j) -> (out[j-1] = st.elements[end].u[3]))
    return out                  # windowed tip out-of-plane deflection, one entry per step
end

# ∂(tip deflection at each window step) / ∂(element compliance) : (k × 36)
J = ForwardDiff.jacobian(windowed_tip, c0)
println("\n∂u_tip/∂(compliance of element $ielem): Jacobian size = ", size(J))
println("  ‖J‖ per window step: ", round.(vec(sqrt.(sum(abs2, J; dims=2))); sigdigits=4))

# ------------------------------------------------------------------
# 3. Directional finite-difference check. Per-entry FD is ill-conditioned here — compliance
#    entries span many orders of magnitude (≈1e-9) — so we probe a single direction scaled by the
#    entry magnitudes and take a small FRACTIONAL central step (mirrors test/test_ad.jl).
# ------------------------------------------------------------------
using Random
v    = randn(MersenneTwister(7), length(c0)) .* abs.(c0)   # scale to each entry's magnitude
frac = 1e-4                                                 # small relative step along v
d_fd = (windowed_tip(c0 .+ frac .* v) .- windowed_tip(c0 .- frac .* v)) ./ (2 * frac)
d_ad = J * v                                               # directional derivative from AD
println("\nDirectional AD-vs-FD along a |c0|-scaled probe (window step k):")
for k in eachindex(d_fd)
    rel = abs(d_ad[k] - d_fd[k]) / max(abs(d_fd[k]), eps())
    println("  step $k:  AD = ", round(d_ad[k]; sigdigits=6),
            "   FD = ", round(d_fd[k]; sigdigits=6),
            "   rel.err = ", round(rel; sigdigits=3))
end
