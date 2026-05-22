#=
Diagnostic: plot the CCBlade BEM residual R(φ) for a cylindrical airfoil
station, the way WATT's `solve_BEMT` actually evaluates it.

Constructs a fake `DS.Airfoil` with a flat polar (Cl ≡ 0, Cd ≡ cd_const at
α = -π, 0, π), wires it into the same `(xv, pv)` packing used by
[src/bemt.jl:73](../src/bemt.jl#L73), and sweeps φ ∈ (ε, π-ε) through
`CCBlade.residual_and_outputs`. The intent is to check empirically whether
the residual actually has a clean root for a cylinder, or whether WATT is
just getting lucky on the NREL 5MW geometry.

Three curves are plotted: no tip correction, Prandtl tip+hub correction,
and the tip+hub correction with the cylinder station moved closer to the
hub (the worst-case F→0 scenario).

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels, CCBlade, FLOWMath
using StructArrays, Plots

const of = OpenFASTTools
const DS = DynamicStallModels

ofpath = joinpath(@__DIR__, "..", "data", "openfast")

# --- Build a flat-polar DS.Airfoil from Cylinder1.dat ---------------------
afinput = of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder1.dat"))
af, _   = of.make_dsairfoil(afinput)

cd_const = afinput.cd[1]   # ≈ 0.50 for Cylinder1
polar_   = [-pi  0.0  cd_const  0.0;
             0.0  0.0  cd_const  0.0;
             pi   0.0  cd_const  0.0]
af_flat  = DS.update_airfoil(af; dsmodel=DS.NoModel(), polar=polar_)

# --- Operating-point: representative NREL 5MW cylinder station ------------
# r ≈ 2.5 m, chord ≈ 3.5 m, rhub = 1.5, rtip = 63.
# Vx = 11.4 m/s, Vy = r*Ω where Ω = Vrated*tsr/(rtip*cos(precone)) ≈ 1.37 rad/s.
r_base     = 2.5
chord      = 3.5
twist      = 13.0 * pi/180
rhub       = 1.5
rtip       = 63.0
precone    = 2.5 * pi/180
Vrated     = 11.4
tsr        = 7.55
Ω          = Vrated*tsr/(rtip*cos(precone))
Vx         = Vrated
rho        = 1.225
mu         = 1.837e-5
asound     = 335.0
pitch      = 0.0
B          = 3
turbine    = true

# --- Sweep R(φ) under three correction settings ---------------------------
phis = collect(range(1e-3, pi - 1e-3, length=2000))

function residual_curve(af, r; tip_corr=nothing)
    xv = [r, chord, twist, rhub, rtip, Vx, r*Ω, rho, pitch, mu, asound]
    pv = (af, B, turbine, nothing, nothing, nothing, tip_corr)
    return [CCBlade.residual_and_outputs(phi, xv, pv)[1] for phi in phis]
end

R_no_corr   = residual_curve(af_flat, r_base; tip_corr=nothing)
R_with_pth  = residual_curve(af_flat, r_base; tip_corr=CCBlade.PrandtlTipHub())
R_near_hub  = residual_curve(af_flat, rhub + 0.2; tip_corr=CCBlade.PrandtlTipHub())

# Brent root (what WATT.solve_BEMT would converge to) for each case
function brent_root(af, r; tip_corr=nothing)
    xv = [r, chord, twist, rhub, rtip, Vx, r*Ω, rho, pitch, mu, asound]
    pv = (af, B, turbine, nothing, nothing, nothing, tip_corr)
    residual(phi) = CCBlade.residual_and_outputs(phi, xv, pv)[1]
    phistar, _ = FLOWMath.brent(residual, 1e-3, pi/2 - 1e-3)
    return phistar
end

phi_brent_no   = brent_root(af_flat, r_base; tip_corr=nothing)
phi_brent_pth  = brent_root(af_flat, r_base; tip_corr=CCBlade.PrandtlTipHub())
phi_brent_near = brent_root(af_flat, rhub + 0.2; tip_corr=CCBlade.PrandtlTipHub())

# CCBlade cylinder-shortcut solution: φ = atan(Vx, Vy), no induction
phi_short_base = atan(Vx, r_base * Ω)
phi_short_near = atan(Vx, (rhub + 0.2) * Ω)

# --- Plot -----------------------------------------------------------------
plt = plot(phis, R_no_corr,                label="no tip/hub corr",       lw=2,
           xlabel="φ (rad)", ylabel="R(φ)", title="BEM residual for a cylinder station",
           legend=:topright)
plot!(plt, phis, R_with_pth,                label="PrandtlTipHub, r=2.5", lw=2)
# plot!(plt, phis, R_near_hub,                label="PrandtlTipHub, r=rhub+0.2", lw=2, linestyle=:dash)
hline!(plt, [0.0], label=false, color=:black, alpha=0.4)

vline!(plt, [phi_brent_no],   label="Brent root (no corr)",            color=1, linestyle=:dot, lw=2)
vline!(plt, [phi_brent_pth],  label="Brent root (PTH, r=2.5)",         color=2, linestyle=:dot, lw=2)
# vline!(plt, [phi_brent_near], label="Brent root (PTH, r=rhub+0.2)",    color=3, linestyle=:dot, lw=2)
vline!(plt, [phi_short_base], label="CCBlade shortcut (r=2.5)",        color=:black, linestyle=:dashdot, lw=2)
# vline!(plt, [phi_short_near], label="CCBlade shortcut (r=rhub+0.2)",   color=:gray,  linestyle=:dashdot, lw=2)

display(plt)

# --- Report sign-changes (rough root locations) ---------------------------
function root_intervals(R)
    out = Tuple{Float64, Float64}[]
    for i in 1:length(R)-1
        if sign(R[i]) != sign(R[i+1]) && isfinite(R[i]) && isfinite(R[i+1])
            push!(out, (phis[i], phis[i+1]))
        end
    end
    return out
end

println("Sign changes (no corr):       ", root_intervals(R_no_corr))
println("Sign changes (PTH, r=2.5):    ", root_intervals(R_with_pth))
println("Sign changes (PTH, near hub): ", root_intervals(R_near_hub))

println()
println("Brent roots:")
println("  no corr           φ = $phi_brent_no")
println("  PTH, r=2.5        φ = $phi_brent_pth")
# println("  PTH, r=rhub+0.2   φ = $phi_brent_near")
# println("CCBlade shortcut:")
println("  CCS, r=2.5        φ = $phi_short_base")
# println("  r=rhub+0.2        φ = $phi_short_near")
