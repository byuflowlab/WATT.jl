#=
Unsteady aerostructural example using the NREL 5MW dissertation blade,
driven by a trained Lux NeuralKoopman surrogate in place of GXBeam.

Mirrors `aerostructural_nrel5mw5seg.jl` but swaps the structural solve:
`initialize_sim` / `run_sim!`  →  `initialize_sim_surrogate` / `run_sim_surrogate!`.

The surrogate parameters are stored in `data/global_model_selstates_...jld2`
as plain NamedTuples of Float32 matrices (Lux's serialization of a
NeuralKoopman model). The forward pass is reimplemented here so the
example needs no Lux dependency.

**Caveat**: this example sets the design vector `x_norm = zeros(nx)` (the
training-set mean in normalized space), because we don't have the
mapping from this specific NREL 5MW blade geometry to the surrogate's
57-D design vector. Treat the resulting trajectory as an integration
demonstration, not a validated structural prediction.

Adam Cardoza
=#

using WATT, OpenFASTTools, DynamicStallModels, GXBeam
using StaticArrays, StructArrays, JLD2
using LinearAlgebra
using FLOWMath
using Printf
using Plots

const of = OpenFASTTools
const DS = DynamicStallModels

datadir = joinpath(@__DIR__, "..", "data")
ofpath  = joinpath(datadir, "openfast")

# ---------------------------------------------------------------------------
# Surrogate forward primitives (Lux-free reimplementation)
# ---------------------------------------------------------------------------

mish(x) = x * tanh(log1p(exp(x)))
sigmoid(x) = one(x) / (one(x) + exp(-x))

# Recursive cast of Lux NamedTuple parameter trees to Float64.
f64(x::AbstractArray{<:AbstractFloat}) = Float64.(x)
f64(x::AbstractArray) = x
f64(x::NamedTuple) = NamedTuple{keys(x)}(map(f64, values(x)))
f64(x::Tuple) = map(f64, x)
f64(x::Real) = Float64(x)
f64(x) = x

# Hermite cubic regrid: y[s_src] → y[s_dst], scalar series.
function hermite_regrid(y::AbstractVector{T}, s_src::AbstractVector,
                        s_dst::AbstractVector) where {T}
    n = length(s_src)
    ds = s_src[2:end] .- s_src[1:end-1]
    m_pair = (y[2:end] .- y[1:end-1]) ./ ds
    m = vcat(m_pair[1:1],
             (m_pair[2:end] .+ m_pair[1:end-1]) ./ T(2),
             m_pair[end:end])
    out = Vector{T}(undef, length(s_dst))
    for (k, sn) in enumerate(s_dst)
        I = clamp(searchsortedfirst(view(s_src, 2:n), sn), 1, n - 1)
        dx = s_src[I+1] - s_src[I]
        t  = (sn - s_src[I]) / dx
        h00 = (1 + 2t) * (1 - t)^2
        h10 = t * (1 - t)^2
        h01 = t^2 * (3 - 2t)
        h11 = t^2 * (t - 1)
        out[k] = h00 * y[I] + h10 * m[I] * dx + h01 * y[I+1] + h11 * m[I+1] * dx
    end
    return out
end

# Apply a FiLM-conditioned MLP (hidden Denses + final out Dense, FiLM on every layer).
# `gammas`, `betas` are tuples of length `length(hidden) + 1`.
function film_mlp_forward(x, ps_mlp, gammas, betas, act)
    h = x
    n_hidden = length(gammas) - 1
    for i in 1:n_hidden
        l = ps_mlp.hidden[Symbol("layer_$i")]
        h = l.weight * h .+ l.bias
        h = gammas[i] .* h .+ betas[i]
        h = act.(h)
    end
    h = ps_mlp.out.weight * h .+ ps_mlp.out.bias
    h = gammas[end] .* h .+ betas[end]
    return h
end

# Apply a FiLM-hypernet (build_mlp: hidden Denses with activation, output Dense linear)
# and split the raw output into per-target-layer (gamma, beta) tuples.
function hyper_forward(x_b, ps_hyper, target_widths::Vector{Int})
    layer_keys = collect(keys(ps_hyper.net))
    n = length(layer_keys)
    h = x_b
    for i in 1:n-1
        l = ps_hyper.net[Symbol("layer_$i")]
        h = l.weight * h .+ l.bias
        h = mish.(h)
    end
    last = ps_hyper.net[Symbol("layer_$n")]
    raw = last.weight * h .+ last.bias

    total = sum(target_widths)
    raw_g = raw[1:total, :]
    raw_b = raw[total+1:end, :]

    gammas = ntuple(length(target_widths)) do i
        offset = sum(target_widths[1:i-1])
        w = target_widths[i]
        1.0 .+ raw_g[offset+1:offset+w, :]
    end
    betas = ntuple(length(target_widths)) do i
        offset = sum(target_widths[1:i-1])
        w = target_widths[i]
        raw_b[offset+1:offset+w, :]
    end
    return gammas, betas
end

function get_eigs(ps_kq_hyper, ps_kq_global, x_b, ncp::Int, n_real::Int,
                  nlatent::Int, svd_buffer::Real)
    gammas, betas = hyper_forward(x_b, ps_kq_hyper, [ncp, ncp, n_real, nlatent])
    g_r, g_th, g_kr, g_q = gammas
    b_r, b_th, b_kr, b_q = betas

    c = 1.0 - svd_buffer
    r_raw_b      = g_r  .* ps_kq_global.r_raw      .+ b_r
    theta_b      = g_th .* ps_kq_global.theta      .+ b_th
    K_real_raw_b = g_kr .* ps_kq_global.K_real_raw .+ b_kr
    Q_raw_b      = g_q  .* ps_kq_global.Q_raw      .+ b_q

    r      = sigmoid.(r_raw_b) .* c
    mu     = r .* cos.(theta_b)
    omega  = r .* sin.(theta_b)
    K_real = tanh.(K_real_raw_b) .* c

    r_pair = reshape(r, 1, size(r, 1), size(r, 2))
    eig_mag_complex = reshape(vcat(r_pair, r_pair), 2 * size(r, 1), size(r, 2))
    eig_mag_real    = abs.(K_real)
    eig_mag         = vcat(eig_mag_complex, eig_mag_real)

    q_max = sqrt.(max.(c^2 .- eig_mag .^ 2, 0.0))
    Q     = tanh.(Q_raw_b) .* q_max

    return mu, omega, K_real, Q
end

function apply_jordan_A(z, mu, omega, K_real, M::Int, ncp::Int)
    B = size(z, 2)
    if ncp > 0
        z1 = z[1:2:M, :]
        z2 = z[2:2:M, :]
        z1n = mu    .* z1 .- omega .* z2
        z2n = omega .* z1 .+ mu    .* z2
        stacked = vcat(reshape(z1n, 1, ncp, B), reshape(z2n, 1, ncp, B))
        Az_c = reshape(stacked, M, B)
    else
        Az_c = similar(z, 0, B)
    end
    if size(K_real, 1) > 0
        Az_r = K_real .* z[M+1:end, :]
        return vcat(Az_c, Az_r)
    else
        return Az_c
    end
end

# ---------------------------------------------------------------------------
# Surrogate normalization constants (must match training script)
# ---------------------------------------------------------------------------

# 18 selected state components: (u_xyz, theta_xyz, V_xyz, Omega_xyz, F_xyz, M_xyz)
const NORM_SCALES_VEC = Float64[
    1e+1, 1e+0, 1e+0,        # u
    1e+1, 1e+1, 1e+1,        # theta
    1e-1, 1e-1, 1e+0,        # V
    1e+0, 1e+0, 1e+0,        # Omega
    1e-6, 1e-6, 1e-6,        # F
    1e-5, 1e-7, 1e-7,        # M
]
const FNORM     = 1e4
const N_REGRID  = 6
const N_SEL     = 18    # selected state components

# ---------------------------------------------------------------------------
# Conditioned surrogate struct + AbstractStructuralSurrogate interface
# ---------------------------------------------------------------------------

struct ConditionedKoopman{TG, TB, TPE, TPD, TPF} <: WATT.AbstractStructuralSurrogate
    # FiLM-baked (gamma, beta) tuples for encoder / decoder / force_nn
    enc_g::TG;   enc_b::TB
    dec_g::TG;   dec_b::TB
    force_g::TG; force_b::TB
    # Baked Jordan blocks (each is a (·, 1) column matrix)
    mu::Matrix{Float64}
    omega::Matrix{Float64}
    K_real::Matrix{Float64}
    Q::Matrix{Float64}
    # MLP parameters (cast to Float64)
    ps_encoder::TPE
    ps_decoder::TPD
    ps_force::TPF
    # Sizes
    n_complex_pairs::Int
    n_complex_dims::Int     # = 2 * ncp
    nlatent::Int
    # Inference-time grids
    s_elem_infer::Vector{Float64}    # nelem_WATT element midpoint arc-length (normalized)
    s_regrid::Vector{Float64}        # range(0, 1, length=N_REGRID)
    s_state_infer::Vector{Float64}   # range(0, 1, length=np) — matches training regrid_state convention
end

function build_conditioned_koopman(jld_path::AbstractString, x_norm::AbstractVector,
                                    assembly::GXBeam.Assembly)
    ckpt = JLD2.load(jld_path)
    # Cast the entire parameter tree to Float64 up front — the model was
    # trained in Float32 but inference is happily Float64 (no AD precision
    # issues at the WATT boundary, and downstream blade/assembly are Float64).
    ps = f64(ckpt["ps"])

    nlatent         = ckpt["nlatent"]
    ncp             = ckpt["n_complex_pairs"]
    nstates         = ckpt["nstates"]
    nforce_regrid   = ckpt["nforce_regrid"]
    svd_buffer      = Float64(ckpt["svd_buffer"])
    @assert nstates       == N_SEL * N_REGRID
    @assert nforce_regrid == 6     * N_REGRID

    M      = 2 * ncp
    n_real = nlatent - M

    # FiLM target widths (must match training script's MLP topology).
    ae_widths    = ckpt["ae_widths"]
    force_widths = ckpt["force_widths"]
    enc_widths        = vcat(collect(ae_widths),         nlatent)
    dec_widths        = vcat(reverse(collect(ae_widths)), nstates)
    force_widths_full = vcat(collect(force_widths),      nlatent)

    x_b = reshape(Float64.(x_norm), :, 1)

    enc_g,   enc_b   = hyper_forward(x_b, ps.enc_hyper,   enc_widths)
    dec_g,   dec_b   = hyper_forward(x_b, ps.dec_hyper,   dec_widths)
    force_g, force_b = hyper_forward(x_b, ps.force_hyper, force_widths_full)

    mu, omega, K_real, Q = get_eigs(ps.kq_hyper, ps.kq_global, x_b,
                                     ncp, n_real, nlatent, svd_buffer)

    # Inference-time grids (Float64)
    nelem_infer = length(assembly.elements)
    elem_r = Float64[norm(assembly.elements[j].x) for j in 1:nelem_infer]
    s_elem_infer = (elem_r .- elem_r[1]) ./ (elem_r[end] - elem_r[1])
    s_regrid     = collect(range(0.0, 1.0, length=N_REGRID))

    np = length(assembly.points)
    s_state_infer = collect(range(0.0, 1.0, length=np))

    return ConditionedKoopman(
        enc_g, enc_b, dec_g, dec_b, force_g, force_b,
        mu, omega, K_real, Q,
        ps.encoder, ps.decoder, ps.force_nn,
        ncp, M, nlatent,
        s_elem_infer, s_regrid, s_state_infer,
    )
end

function WATT.encode_initial(surr::ConditionedKoopman, u0_struct::WATT.SurrogateAssemblyState)
    np = length(u0_struct.points)
    u_raw = Matrix{Float64}(undef, np, N_SEL)
    for j in 1:np
        p = u0_struct.points[j]
        u_raw[j, 1:3]   = p.u
        u_raw[j, 4:6]   = p.theta
        u_raw[j, 7:9]   = p.V
        u_raw[j, 10:12] = p.Omega
        u_raw[j, 13:15] = p.F
        u_raw[j, 16:18] = p.M
    end
    u_regrid = Matrix{Float64}(undef, N_REGRID, N_SEL)
    for k in 1:N_SEL
        # Hermite-regrid each component from np structural points → 6 surrogate
        # nodes, then SCALE (multiply) by NORM_SCALES_VEC[k] to enter the
        # network's scaled-units convention. Matches training data prep:
        #   u_scaled = u_phys .* scale_vec
        u_regrid[:, k] = hermite_regrid(u_raw[:, k], surr.s_state_infer, surr.s_regrid) .*
                         NORM_SCALES_VEC[k]
    end
    u_flat = reshape(u_regrid, N_REGRID * N_SEL, 1)
    return film_mlp_forward(u_flat, surr.ps_encoder, surr.enc_g, surr.enc_b, mish)
end

function WATT.step_latent(surr::ConditionedKoopman, z, f_per_element::AbstractMatrix)
    f_regrid = Matrix{Float64}(undef, N_REGRID, 6)
    for k in 1:6
        # Force scaling: divide by FNORM (=1e4) to match training (`f_raw ./ fnorm`),
        # then Hermite-regrid each load component from nelem element midpoints
        # → 6 surrogate nodes using s_elem_infer (the WATT-side element grid).
        f_regrid[:, k] = hermite_regrid(Float64.(f_per_element[:, k]) ./ FNORM,
                                         surr.s_elem_infer, surr.s_regrid)
    end
    f_flat = reshape(f_regrid, 6 * N_REGRID, 1)
    fi_enc = film_mlp_forward(f_flat, surr.ps_force, surr.force_g, surr.force_b, mish)
    return apply_jordan_A(z, surr.mu, surr.omega, surr.K_real,
                          surr.n_complex_dims, surr.n_complex_pairs) .+ surr.Q .* fi_enc
end

function WATT.decode(surr::ConditionedKoopman, z)
    # Promote the output to Float64 so it matches the host blade/assembly's
    # element type (run_sim_surrogate! infers history `TF` from `blade.c[1]`).
    u_flat   = film_mlp_forward(z, surr.ps_decoder, surr.dec_g, surr.dec_b, mish)
    u_scaled = reshape(u_flat, N_REGRID, N_SEL)   # (6 surrogate nodes × 18 components) in SCALED units
    np = length(surr.s_state_infer)
    u_pts = Matrix{Float64}(undef, np, N_SEL)
    for k in 1:N_SEL
        # DESCALE: divide by NORM_SCALES_VEC[k] to leave scaled-units land and
        # return to physical units. Then Hermite-regrid each component from
        # 6 surrogate nodes → np structural points.
        # (Scale and regrid commute since NORM_SCALES_VEC[k] is a scalar.)
        u_pts[:, k] = hermite_regrid(u_scaled[:, k] ./ NORM_SCALES_VEC[k],
                                      surr.s_regrid, surr.s_state_infer)
    end
    points = Vector{WATT.SurrogatePointState{Float64}}(undef, np)
    for j in 1:np
        points[j] = WATT.SurrogatePointState{Float64}(
            SVector{3,Float64}(u_pts[j, 1],  u_pts[j, 2],  u_pts[j, 3]),   # u
            SVector{3,Float64}(u_pts[j, 4],  u_pts[j, 5],  u_pts[j, 6]),   # theta
            SVector{3,Float64}(u_pts[j, 7],  u_pts[j, 8],  u_pts[j, 9]),   # V
            SVector{3,Float64}(u_pts[j, 10], u_pts[j, 11], u_pts[j, 12]),  # Omega
            SVector{3,Float64}(u_pts[j, 13], u_pts[j, 14], u_pts[j, 15]),  # F
            SVector{3,Float64}(u_pts[j, 16], u_pts[j, 17], u_pts[j, 18]),  # M
        )
    end
    return WATT.SurrogateAssemblyState{Float64}(points)
end

# ---------------------------------------------------------------------------
# Blade / Rotor / Assembly  (identical to aerostructural_nrel5mw5seg.jl)
# ---------------------------------------------------------------------------
@load joinpath(datadir, "nrel5mw_5seg.jld2") rvec cvec twistvec le_loc compliance_list mass_list points xp start stop Rhub Rtip precone raf afidx polars cylinder_mask
nr = length(rvec)

aftype_names = ("Cylinder1.dat", "Cylinder2.dat", "DU40_A17.dat", "DU35_A17.dat",
                "DU30_A17.dat",  "DU25_A17.dat",  "DU21_A17.dat", "NACA64_A17.dat")
aftypes = [of.read_airfoilinput(joinpath(ofpath, "airfoils", name)) for name in aftype_names]

af_idx = of.integerfit(raf, afidx, rvec)
afs    = aftypes[af_idx]

dsairfoils = StructArray{DS.Airfoil}(undef, nr)
xcp        = Vector{Float64}(undef, nr)
for i = 1:nr
    dsairfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
    if cylinder_mask[i]
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; dsmodel=DS.NoModel(), polar=polars[i])
    else
        dsairfoils[i] = DS.update_airfoil(dsairfoils[i]; polar=polars[i])
    end
end

blade = WATT.Blade(rvec, cvec, twistvec, xcp, dsairfoils; rhub=Rhub, rtip=Rtip, precone)

B     = 3
hubHt = 90.0
rotor = WATT.Rotor(B, hubHt, true; tilt=0.0, yaw=0.0)

assembly = GXBeam.Assembly(points, start, stop; compliance=compliance_list,
                           midpoints=xp, mass=mass_list)

rho      = 1.225
mu_air   = 1.837e-5
a        = 335.0
shearexp = 0.2
Vrated   = 11.4
tsr      = 7.55
rotorR   = Rtip * cos(precone)
omega    = Vrated * tsr / rotorR

turbfile = joinpath(ofpath, "TurbSim.dat")
env      = environment(turbfile, rho, mu_air, a, omega, shearexp)

# ---------------------------------------------------------------------------
# Build the conditioned surrogate
# ---------------------------------------------------------------------------
surr_path = joinpath(datadir,
    "global_model_selstates_0_1_2_6_7_8_12_13_14_18_19_20_24_25_26_27_28_29_x_jordan_film_globalKQ_svdstable_multidata_lux_20260522_121111.jld2")

# x_norm = zeros(nx) → training-set mean in normalized space. Replace with the
# normalized design vector that actually corresponds to this NREL 5MW blade if
# you have the design-to-x mapping from the surrogate's training pipeline.
ckpt_meta = JLD2.load(surr_path, "nx")
x_norm    = zeros(Float64, ckpt_meta)

surr = build_conditioned_koopman(surr_path, x_norm, assembly)

println("Conditioned surrogate built: nlatent=$(surr.nlatent)  ncp=$(surr.n_complex_pairs)  ",
        "nelem_infer=$(length(surr.s_elem_infer))  np=$(length(surr.s_state_infer))")

# ---------------------------------------------------------------------------
# Simulate via the surrogate path
# ---------------------------------------------------------------------------
tvec = collect(0:0.05:10.0)   # dt matches training (0.05 s)

### Baseline GXBeam simulation FIRST — its IC state is needed to seed the
### surrogate, since the training data's first frame is the GXBeam
### structural equilibrium (not a zero-state).
println("\n=== Running GXBeam baseline ===")
aerostates_gx, gxhistory, mesh_gx = WATT.initialize_sim(blade, assembly, tvec; verbose=true)
WATT.run_sim!(rotor, blade, mesh_gx, env, tvec, aerostates_gx, gxhistory; verbose=true)

# Convert gxhistory[1] → SurrogateAssemblyState, matching the training-data
# extraction convention in unsteady_analysis (Cardoza2026):
#   u/theta/V/Omega come straight from points
#   F/M at the root point = -reaction (negated)
#   F/M at interior points (j≥2) = Akima(element midpoint Fi/Mi)
function gx_state_to_surrogate(gx::GXBeam.AssemblyState, assembly::GXBeam.Assembly)
    np = length(assembly.points)
    nelem = length(assembly.elements)
    span_elem = Float64[norm(e.x) for e in assembly.elements]
    span_node = Float64[norm(p) for p in assembly.points]

    Fi = zeros(Float64, nelem, 3)
    Mi = zeros(Float64, nelem, 3)
    for j in 1:nelem
        Fi[j, :] = gx.elements[j].Fi
        Mi[j, :] = gx.elements[j].Mi
    end

    F_pts = zeros(Float64, np, 3)
    M_pts = zeros(Float64, np, 3)
    F_pts[1, :] = -gx.points[1].F
    M_pts[1, :] = -gx.points[1].M
    for k in 1:3
        Fi_spl = FLOWMath.Akima(span_elem, Fi[:, k])
        Mi_spl = FLOWMath.Akima(span_elem, Mi[:, k])
        for j in 2:np
            F_pts[j, k] = Fi_spl(span_node[j])
            M_pts[j, k] = Mi_spl(span_node[j])
        end
    end

    points = Vector{WATT.SurrogatePointState{Float64}}(undef, np)
    for j in 1:np
        p = gx.points[j]
        points[j] = WATT.SurrogatePointState{Float64}(
            SVector{3,Float64}(p.u...),
            SVector{3,Float64}(p.theta...),
            SVector{3,Float64}(p.V...),
            SVector{3,Float64}(p.Omega...),
            SVector{3,Float64}(F_pts[j, 1], F_pts[j, 2], F_pts[j, 3]),
            SVector{3,Float64}(M_pts[j, 1], M_pts[j, 2], M_pts[j, 3]),
        )
    end
    return WATT.SurrogateAssemblyState{Float64}(points)
end

u0_struct_gx = gx_state_to_surrogate(gxhistory[1], assembly)



### Now run the surrogate, seeded with the GXBeam IC.
println("\n=== Running surrogate ===")
aerostates, surr_history, mesh =
    WATT.initialize_sim_surrogate(blade, assembly, tvec; verbose=true)

WATT.run_sim_surrogate!(rotor, blade, mesh, env, tvec, aerostates, surr_history, surr; u0_struct=u0_struct_gx, verbose=true)
# WATT.run_sim_surrogate!(rotor, blade, mesh, env, tvec, aerostates, surr_history, surr; verbose=true)

# ---------------------------------------------------------------------------
# Extract comparable quantities
# ---------------------------------------------------------------------------
# Surrogate
tip_def_z_surr = [s.points[end].u[3] for s in surr_history]
tip_F_y_surr   = [s.points[end].F[2] for s in surr_history]
root_M_x_surr  = [s.points[1].M[1]   for s in surr_history]

# GXBeam baseline.
# Training-data sign convention (see Cardoza2026 unsteady_analysis):
#   u[:, 1, F/M] = -gxhistory.points[1].F/M    (root reaction, sign-flipped)
#   u[:, j, F/M] = Akima-interp of element midpoint Fi/Mi at structural nodes (j ≥ 2)
# Use elements[end].Fi/Mi at the tip and -points[1].M at the root to match the
# surrogate's per-point output convention.
tip_def_z_gx = [gx.points[end].u[3]   for gx in gxhistory]
tip_F_y_gx   = [gx.elements[end].Fi[2] for gx in gxhistory]
root_M_x_gx  = [-gx.points[1].M[1]    for gx in gxhistory]


############ Compare initial states

z0 = encode_initial(surr, u0_struct_gx)
u0 = decode(surr, z0)


println("Initial state comparison (surrogate decode of surrogate encode of GXBeam IC):")
println("  - u0_struct_gx (GXBeam IC): $(u0_struct_gx.points[end].u)")
println("  - u0 (surrogate decode):    $(u0.points[end].u)")
println("  - surrogate initial state:  $(surr_history[1].points[end].u)")
println("  - GXBeam initial state:     $(gxhistory[1].points[end].u)")












# ---------------------------------------------------------------------------
# Plot — surrogate vs GXBeam overlaid
# ---------------------------------------------------------------------------
tipplt = plot(tvec, tip_def_z_gx, label="GXBeam",
              xlabel="Time (s)", ylabel="Tip deflection z (m)",
              title="Out-of-plane tip deflection", lw=2, color=:black)
plot!(tipplt, tvec, tip_def_z_surr, label="Surrogate", lw=2, color=:crimson, linestyle=:dash)

loadplt = plot(blade.r, aerostates_gx.Fx[end, :], label="Fx (GXBeam)", lw=2, color=:black)
plot!(loadplt, blade.r, aerostates_gx.Fy[end, :], label="Fy (GXBeam)", lw=2, color=:gray)
plot!(loadplt, blade.r, aerostates.Fx[end, :], label="Fx (surrogate)", lw=2,
      color=:crimson, linestyle=:dash)
plot!(loadplt, blade.r, aerostates.Fy[end, :], label="Fy (surrogate)", lw=2,
      color=:orange,  linestyle=:dash)
plot!(loadplt, title="Final-step spanwise aero loads",
      xlabel="Blade span (m)", ylabel="Force (N/m)")

mplt = plot(tvec, root_M_x_gx, label="GXBeam",
            xlabel="Time (s)", ylabel="Root M_x (N·m)",
            title="Root internal moment (x)", lw=2, color=:black)
plot!(mplt, tvec, root_M_x_surr, label="Surrogate", lw=2, color=:crimson, linestyle=:dash)

fplt = plot(tvec, tip_F_y_gx, label="GXBeam",
            xlabel="Time (s)", ylabel="Tip F_y (N)",
            title="Tip internal force (y)", lw=2, color=:black)
plot!(fplt, tvec, tip_F_y_surr, label="Surrogate", lw=2, color=:crimson, linestyle=:dash)

plt = plot(tipplt, loadplt, mplt, fplt, layout=(4, 1), size=(800, 1100))
display(plt)
