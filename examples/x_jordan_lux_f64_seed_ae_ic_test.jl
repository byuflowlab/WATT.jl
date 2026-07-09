# AE round-trip diagnostic for the trained seed-only Float64 model.
#
# Sister of `x_jordan_lux_f64_seed.jl`. Reuses the same data preprocessing,
# model struct definitions, and Lux forward pass, but DOES NOT TRAIN. It:
#
#   1. Loads the trained checkpoint (`MODEL_JLD2`).
#   2. Loads case 1 of the training HDF5 and runs the same preprocessing as
#      training (select 18 components → scale → regrid 6 nodes).
#   3. Computes the autoencoder round-trip on the IC frame:
#         z0     = encoder(u_train[:, 1, 1:1], enc_g, enc_b)
#         u_rec  = decoder(z0,                 dec_g, dec_b)
#      and prints max|u_rec - u_b_train|.
#   4. Saves u_b_train, x_b, enc_g/b, dec_g/b, z0, u_rec, and the kq/force
#      conditioning to `OUT_JLD2`, so the WATT reimplementation can load these
#      ground-truth tensors and compare byte-for-byte at every stage of its
#      own forward pass.
#
# To run:  julia --project=examples examples/x_jordan_lux_f64_seed_ae_ic_test.jl

using Lux, LuxCUDA, Optimisers
using Lux: mish, sigmoid_fast, tanh_fast
using HDF5, JLD2
using Random, Statistics, LinearAlgebra
using Printf

# ---------------------------------------------------------------------------
# User Configuration
# ---------------------------------------------------------------------------
const ALL_STATE_TYPES = ("u", "udot", "theta", "thetadot",
                         "V", "Vdot", "Omega", "Omegadot", "F", "M")

const STATE_COMPONENT_INDICES =
    (0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20, 24, 25, 26, 27, 28, 29)

const NORM_SCALES = Dict(
    "u"        => Float64[1e+1, 1e+0, 1e+0],
    "udot"     => Float64[1e+1, 1e+0, 1e+0],
    "theta"    => Float64[1e+1, 1e+1, 1e+1],
    "thetadot" => Float64[1e+0, 1e+0, 1e+0],
    "V"        => Float64[1e-1, 1e-1, 1.0],
    "Vdot"     => Float64[1e-1, 1e-1, 1e-1],
    "Omega"    => Float64[1e+0, 1e+0, 1e+0],
    "Omegadot" => Float64[1e+0, 1e+0, 1e+0],
    "F"        => Float64[1e-6, 1e-6, 1e-6],
    "M"        => Float64[1e-5, 1e-7, 1e-7],
)

const LATENT_MULT    = 6
const N_REGRID_NODES = 6

const AE_WIDTHS    = (256, 256, 256, 256)
const AE_ACT       = mish
const FORCE_WIDTHS = (256, 256, 256, 256)
const FORCE_ACT    = mish
const HYPER_WIDTHS = (128, 128)
const HYPER_ACT    = mish
const SVD_BUFFER   = 0.1

# const DATA_FILES = [
#     "/Users/adamcardoza/repos/Cardoza2026_FastBladeOpt/exploring/unsteadywind/allstates/serial_xencdec/uw_boundary_random_cases_993.h5",
# ]
const DATA_FILES = [
    "/Users/adamcardoza/repos/Cardoza2026_FastBladeOpt/exploring/unsteadywind/uw_boundary_random_cases_17_16.h5",
]
const MAX_CASES = 1

# const MODEL_JLD2 = joinpath(@__DIR__, "..", "data",
#     "model_f64_seed_20260526_115153.jld2")
const MODEL_JLD2 = joinpath(@__DIR__, "..", "data",
    "model_new.jld2")
# const OUT_JLD2   = joinpath(@__DIR__, "..", "data",
#     "ae_ic_diagnostic.jld2")
const OUT_JLD2   = joinpath(@__DIR__, "..", "data",
    "ae_ic_diagnostic_new.jld2")

const dev  = cpu_device()
const cdev = cpu_device()

println("=== AE IC round-trip diagnostic ===")
println("checkpoint: $MODEL_JLD2")
println("output:     $OUT_JLD2")

# ---------------------------------------------------------------------------
# Hermite spline helpers  (copied verbatim from training)
# ---------------------------------------------------------------------------
function h_poly(t::AbstractVector)
    A = Float64[
        1 0 -3  2;
        0 1 -2  1;
        0 0  3 -2;
        0 0 -1  1
    ]
    tt = (ones(eltype(t), length(t)), t, t .^ 2, t .^ 3)
    return ntuple(i -> A[i, 1] .* tt[1] .+ A[i, 2] .* tt[2] .+
                       A[i, 3] .* tt[3] .+ A[i, 4] .* tt[4], 4)
end

function regrid_with_positions(data::AbstractArray{T,3}, s_old::AbstractVector,
                                n_nodes_out::Int) where {T}
    n_cases, n_tsteps, n_nodes = size(data)
    s_old_f = convert(Vector{T}, collect(s_old))
    s_new = collect(range(zero(T), one(T), length=n_nodes_out))

    N = n_cases * n_tsteps
    y_flat = reshape(data, N, n_nodes)

    ds = s_old_f[2:end] .- s_old_f[1:end-1]
    m_pair = (y_flat[:, 2:end] .- y_flat[:, 1:end-1]) ./ reshape(ds, 1, :)
    m = hcat(m_pair[:, 1:1],
             (m_pair[:, 2:end] .+ m_pair[:, 1:end-1]) ./ T(2),
             m_pair[:, end:end])

    I = [clamp(searchsortedfirst(view(s_old_f, 2:n_nodes), sn), 1, n_nodes - 1)
         for sn in s_new]
    dx = s_old_f[I .+ 1] .- s_old_f[I]
    t_norm = (s_new .- s_old_f[I]) ./ dx
    hh = h_poly(t_norm)

    y_new = transpose(hh[1]) .* y_flat[:, I] .+
            transpose(hh[2]) .* m[:, I] .* transpose(dx) .+
            transpose(hh[3]) .* y_flat[:, I .+ 1] .+
            transpose(hh[4]) .* m[:, I .+ 1] .* transpose(dx)
    return reshape(y_new, n_cases, n_tsteps, n_nodes_out)
end

function regrid_state(data::AbstractArray{T,3}, n_nodes_out::Int=N_REGRID_NODES) where {T}
    n_nodes = size(data, 3)
    s_old = collect(range(zero(T), one(T), length=n_nodes))
    return regrid_with_positions(data, s_old, n_nodes_out)
end

# ---------------------------------------------------------------------------
# Model structs  (copied verbatim from training)
# ---------------------------------------------------------------------------
function build_mlp(input_size::Int, output_size::Int, widths, activation)
    layers = Lux.AbstractLuxLayer[]
    prev = input_size
    for w in widths
        push!(layers, Dense(prev => w, activation))
        prev = w
    end
    push!(layers, Dense(prev => output_size))
    return Chain(layers...)
end

struct FiLMMLP{HT,O,A} <: Lux.AbstractLuxLayer
    hidden::HT
    out::O
    act::A
    film_widths::Vector{Int}
end

function FiLMMLP(input_size::Int, output_size::Int, widths, activation)
    denses = Lux.AbstractLuxLayer[]
    prev = input_size
    for w in widths
        push!(denses, Dense(prev => w))
        prev = w
    end
    hidden = Tuple(denses)
    out = Dense(prev => output_size)
    return FiLMMLP(hidden, out, activation, vcat(collect(widths), output_size))
end

function Lux.initialparameters(rng::AbstractRNG, l::FiLMMLP)
    h = NamedTuple{ntuple(i -> Symbol("layer_$i"), length(l.hidden))}(
        ntuple(i -> Lux.initialparameters(rng, l.hidden[i]), length(l.hidden))
    )
    return (hidden = h, out = Lux.initialparameters(rng, l.out))
end

function Lux.initialstates(rng::AbstractRNG, l::FiLMMLP)
    h = NamedTuple{ntuple(i -> Symbol("layer_$i"), length(l.hidden))}(
        ntuple(i -> Lux.initialstates(rng, l.hidden[i]), length(l.hidden))
    )
    return (hidden = h, out = Lux.initialstates(rng, l.out))
end

Lux.parameterlength(l::FiLMMLP) =
    sum(Lux.parameterlength, l.hidden) + Lux.parameterlength(l.out)
Lux.statelength(l::FiLMMLP) =
    sum(Lux.statelength, l.hidden) + Lux.statelength(l.out)

function (l::FiLMMLP)(input::Tuple, ps, st)
    x, gammas, betas = input
    h = x
    for i in eachindex(l.hidden)
        key = Symbol("layer_$i")
        h, _ = l.hidden[i](h, ps.hidden[key], st.hidden[key])
        h = gammas[i] .* h .+ betas[i]
        h = l.act.(h)
    end
    h, _ = l.out(h, ps.out, st.out)
    h = gammas[end] .* h .+ betas[end]
    return h, st
end

struct FiLMHyper{C} <: Lux.AbstractLuxLayer
    net::C
    target_widths::Vector{Int}
end

function FiLMHyper(nx::Int, target_widths, hyper_widths, hyper_act)
    tw = collect(target_widths)
    total = sum(tw)
    net = build_mlp(nx, 2 * total, hyper_widths, hyper_act)
    return FiLMHyper(net, tw)
end

Lux.initialparameters(rng::AbstractRNG, l::FiLMHyper) =
    (net = Lux.initialparameters(rng, l.net),)
Lux.initialstates(rng::AbstractRNG, l::FiLMHyper) =
    (net = Lux.initialstates(rng, l.net),)
Lux.parameterlength(l::FiLMHyper) = Lux.parameterlength(l.net)
Lux.statelength(l::FiLMHyper) = Lux.statelength(l.net)

function (l::FiLMHyper)(x_b, ps, st)
    raw, _ = l.net(x_b, ps.net, st.net)
    total = sum(l.target_widths)
    raw_g = raw[1:total, :]
    raw_b = raw[total+1:end, :]
    gammas = ntuple(length(l.target_widths)) do i
        offset = sum(l.target_widths[1:i-1])
        w = l.target_widths[i]
        1.0 .+ raw_g[offset+1 : offset+w, :]
    end
    betas = ntuple(length(l.target_widths)) do i
        offset = sum(l.target_widths[1:i-1])
        w = l.target_widths[i]
        raw_b[offset+1 : offset+w, :]
    end
    return (gammas, betas), st
end

struct KQGlobal <: Lux.AbstractLuxLayer
    ncp::Int
    n_real::Int
    nlatent::Int
end

function Lux.initialparameters(rng::AbstractRNG, l::KQGlobal)
    s = 1.0 / sqrt(Float64(l.nlatent))
    runif(n) = (rand(rng, Float64, n) .* 2.0 .- 1.0) .* s
    return (r_raw      = runif(l.ncp),
            theta      = runif(l.ncp),
            K_real_raw = runif(l.n_real),
            Q_raw      = runif(l.nlatent))
end

Lux.initialstates(::AbstractRNG, ::KQGlobal) = NamedTuple()
Lux.parameterlength(l::KQGlobal) = 2 * l.ncp + l.n_real + l.nlatent
Lux.statelength(::KQGlobal) = 0

struct NeuralKoopman{E,D,F,EH,DH,FH,KH,KQ} <:
       Lux.AbstractLuxContainerLayer{(:encoder, :decoder, :force_nn,
                                      :enc_hyper, :dec_hyper, :force_hyper,
                                      :kq_hyper, :kq_global)}
    encoder::E
    decoder::D
    force_nn::F
    enc_hyper::EH
    dec_hyper::DH
    force_hyper::FH
    kq_hyper::KH
    kq_global::KQ
    nstates::Int
    nlatent::Int
    nforce_regrid::Int
    nx::Int
    n_complex_pairs::Int
    n_complex_dims::Int
    n_real_dims::Int
    svd_buffer::Float64
end

function NeuralKoopman(nstates, nlatent, nforce_regrid, nx, n_complex_pairs,
                       ae_widths, ae_act, force_widths, force_act,
                       hyper_widths, hyper_act, svd_buffer)
    M = 2 * n_complex_pairs
    n_real = nlatent - M
    @assert n_real >= 0

    encoder  = FiLMMLP(nstates,       nlatent,       ae_widths,             ae_act)
    decoder  = FiLMMLP(nlatent,       nstates,       reverse(ae_widths),    ae_act)
    force_nn = FiLMMLP(nforce_regrid, nlatent,       force_widths,          force_act)

    enc_hyper   = FiLMHyper(nx, encoder.film_widths,  hyper_widths, hyper_act)
    dec_hyper   = FiLMHyper(nx, decoder.film_widths,  hyper_widths, hyper_act)
    force_hyper = FiLMHyper(nx, force_nn.film_widths, hyper_widths, hyper_act)

    kq_hyper  = FiLMHyper(nx, [n_complex_pairs, n_complex_pairs, n_real, nlatent],
                          hyper_widths, hyper_act)
    kq_global = KQGlobal(n_complex_pairs, n_real, nlatent)

    return NeuralKoopman(encoder, decoder, force_nn,
                         enc_hyper, dec_hyper, force_hyper,
                         kq_hyper, kq_global,
                         nstates, nlatent, nforce_regrid, nx,
                         n_complex_pairs, M, n_real, Float64(svd_buffer))
end

function get_eigs(model::NeuralKoopman, ps, st, x_b)
    c = 1.0 - model.svd_buffer
    (gammas, betas), _ = model.kq_hyper(x_b, ps.kq_hyper, st.kq_hyper)
    g_r, g_th, g_kr, g_q = gammas
    b_r, b_th, b_kr, b_q = betas

    r_raw_b      = g_r  .* ps.kq_global.r_raw      .+ b_r
    theta_b      = g_th .* ps.kq_global.theta      .+ b_th
    K_real_raw_b = g_kr .* ps.kq_global.K_real_raw .+ b_kr
    Q_raw_b      = g_q  .* ps.kq_global.Q_raw      .+ b_q

    r      = sigmoid_fast.(r_raw_b) .* c
    mu     = r .* cos.(theta_b)
    omega  = r .* sin.(theta_b)
    K_real = tanh_fast.(K_real_raw_b) .* c

    r_pair = reshape(r, 1, size(r, 1), size(r, 2))
    eig_mag_complex = reshape(vcat(r_pair, r_pair), 2 * size(r, 1), size(r, 2))
    eig_mag_real    = abs.(K_real)
    eig_mag         = vcat(eig_mag_complex, eig_mag_real)

    q_max = sqrt.(max.(c^2 .- eig_mag .^ 2, 0.0))
    Q     = tanh_fast.(Q_raw_b) .* q_max

    return mu, omega, K_real, Q
end

# ---------------------------------------------------------------------------
# Load data (case 1 only) — same preprocessing path as training
# ---------------------------------------------------------------------------
println("\nLoading dataset...")

function _load_one(path)
    h5open(path, "r") do fid
        (read(fid, "u"), read(fid, "f"), read(fid, "y"), read(fid, "x"))
    end
end

let
    parts = [_load_one(p) for p in DATA_FILES]
    u0, f0, y0, x0 = parts[1]
    println("  file 1: ncase=$(size(u0,1))  ntime=$(size(u0,2))  npoint=$(size(u0,3))")
    global udata_raw = parts[1][1]
    global fdata_raw = parts[1][2]
    global ydata_raw = parts[1][3]
    global xdata_raw = parts[1][4]
end
GC.gc()

udata_raw = udata_raw[1:MAX_CASES, :, :, :]
fdata_raw = fdata_raw[1:MAX_CASES, :, :, :]
ydata_raw = ydata_raw[1:MAX_CASES, :]
xdata_raw = xdata_raw[1:MAX_CASES, :]

const ntimes_full = size(udata_raw, 2)
const skipidx     = 1
const nelem       = size(fdata_raw, 3)
const npoint      = nelem + 1
const dt          = Float64(ydata_raw[1, end])

const fnorm = 1.0e4
const node_r_phys = Float64.(ydata_raw[1, 1:npoint])
const elem_r_phys = 0.5 .* (node_r_phys[1:end-1] .+ node_r_phys[2:end])
const r_min = elem_r_phys[1]
const r_max = elem_r_phys[end]
const s_elem = (elem_r_phys .- r_min) ./ (r_max - r_min)

let
    f_raw = Float64.(fdata_raw[:, 1:ntimes_full, :, :])[:, 1:skipidx:end, :, :] ./ fnorm
    fcomps = [regrid_with_positions(f_raw[:, :, :, k], s_elem, N_REGRID_NODES) for k in 1:6]
    global f_concat = cat(fcomps...; dims=3)
end
const nforce_regrid = 6 * N_REGRID_NODES

const n_comps   = length(STATE_COMPONENT_INDICES)
const scale_vec = Float64[NORM_SCALES[ALL_STATE_TYPES[k ÷ 3 + 1]][k % 3 + 1]
                          for k in STATE_COMPONENT_INDICES]

let
    state_idx_1b = collect(STATE_COMPONENT_INDICES) .+ 1
    u_sel = Float64.(udata_raw[:, 1:ntimes_full, :, state_idx_1b])
    sv = reshape(scale_vec, 1, 1, 1, :)
    u_scaled = u_sel .* sv
    comps = [regrid_state(u_scaled[:, :, :, k], N_REGRID_NODES) for k in 1:n_comps]
    global u_concat = cat(comps...; dims=3)
end

const ncase_total = size(u_concat, 1)
const ntimes      = size(u_concat, 2)
const nstates     = size(u_concat, 3)

u_all = permutedims(u_concat, (3, 2, 1))
f_all = permutedims(f_concat, (3, 2, 1))
x_all = permutedims(Float64.(xdata_raw), (2, 1))
const nx = size(x_all, 1)

u_concat = nothing; f_concat = nothing
udata_raw = nothing; fdata_raw = nothing; ydata_raw = nothing; xdata_raw = nothing
GC.gc()

u_train = u_all
f_train = f_all
x_train = x_all   # seed-only: no normalization (x_mean = 0, x_std = 1)

println("u_train: $(size(u_train))  f_train: $(size(f_train))  x_train: $(size(x_train))")
println("dt = $dt")

# ---------------------------------------------------------------------------
# Build model + load trained weights
# ---------------------------------------------------------------------------
const nlatent         = nstates * LATENT_MULT
const n_complex_pairs = nlatent ÷ 4

model = NeuralKoopman(nstates, nlatent, nforce_regrid, nx, n_complex_pairs,
                      AE_WIDTHS, AE_ACT, FORCE_WIDTHS, FORCE_ACT,
                      HYPER_WIDTHS, HYPER_ACT, SVD_BUFFER)

rng = Random.default_rng(); Random.seed!(rng, 0)
ps, st = Lux.setup(rng, model)
ps = ps |> dev
st = st |> dev

println("\nLoading trained ps from checkpoint...")
let
    ckpt = JLD2.load(MODEL_JLD2)
    sd = ckpt["ps"]
    for k in (:encoder, :decoder, :force_nn,
              :enc_hyper, :dec_hyper, :force_hyper,
              :kq_hyper, :kq_global)
        ps_k = getfield(ps, k)
        sd_k = getfield(sd, k) |> dev
        Lux.fmap((dst, src) -> (dst .= src; dst), ps_k, sd_k)
    end
    println("  ps loaded (encoder/decoder/force/4 hyper-nets + kq_global)")
end

# ---------------------------------------------------------------------------
# AE round-trip on the IC frame of case 1
# ---------------------------------------------------------------------------
println("\n=== AE round-trip on IC frame (case 1) ===")

case_idx  = 1
u_b_train = reshape(u_train[:, 1, case_idx], :, 1)    # (108, 1) scaled, regridded
x_b       = x_train[:, case_idx:case_idx]              # (nx, 1)

(enc_g, enc_b),     _ = model.enc_hyper(x_b,   ps.enc_hyper,   st.enc_hyper)
(dec_g, dec_b),     _ = model.dec_hyper(x_b,   ps.dec_hyper,   st.dec_hyper)
(force_g, force_b), _ = model.force_hyper(x_b, ps.force_hyper, st.force_hyper)
mu, omega, K_real, Q = get_eigs(model, ps, st, x_b)

z0,    _ = model.encoder((u_b_train, enc_g, enc_b), ps.encoder, st.encoder)
u_rec, _ = model.decoder((z0,        dec_g, dec_b), ps.decoder, st.decoder)

res_inf = maximum(abs.(u_rec .- u_b_train))
res_mse = mean((u_rec .- u_b_train) .^ 2)

println(@sprintf("max|u_rec - u_b_train| (scaled units) = %.6e", res_inf))
println(@sprintf("mean((u_rec - u_b_train)^2)           = %.6e", res_mse))
println()

# Indices 1:6 = u_x at 6 regrid nodes,  7:12 = u_y,  13:18 = u_z
# (then theta_xyz, V_xyz, Omega_xyz, F_xyz, M_xyz at 6 nodes each).
function _show_band(label::String, rng::UnitRange{Int})
    println("$label")
    println(@sprintf("  truth : %s", string(round.(vec(u_b_train[rng]); digits=6))))
    println(@sprintf("  pred  : %s", string(round.(vec(u_rec[rng]);     digits=6))))
    println(@sprintf("  Δ     : %s", string(round.(vec(u_rec[rng] .- u_b_train[rng]); digits=6))))
end

_show_band("u_x at 6 nodes (1:6)",     1:6)
_show_band("u_y at 6 nodes (7:12)",    7:12)
_show_band("u_z at 6 nodes (13:18)",  13:18)
_show_band("V_y at 6 nodes (43:48)", 43:48)
_show_band("F_x at 6 nodes (73:78)", 73:78)
_show_band("M_x at 6 nodes (91:96)", 91:96)

# ---------------------------------------------------------------------------
# Dump everything WATT needs for a byte-comparison
# ---------------------------------------------------------------------------
# Saved fields:
#   u_b_train   : (108, 1) scaled+regridded IC frame  (encoder input)
#   x_b         : (nx, 1)  design vector              (FiLM input)
#   enc_g, enc_b: tuples (5×) of (width_i, 1)         (encoder FiLM)
#   dec_g, dec_b: tuples (5×) of (width_i, 1)         (decoder FiLM)
#   force_g, force_b : tuples (5×) of (width_i, 1)    (force_nn FiLM)
#   mu, omega, K_real, Q : Jordan blocks               (rollout)
#   z0          : (nlatent, 1)                         (encoder output)
#   u_rec       : (108, 1)                              (decoder output)
#
# Tuples are saved as plain `Vector` of matrices; the WATT side just iterates.

JLD2.jldsave(OUT_JLD2;
    u_b_train = collect(u_b_train),
    x_b       = collect(x_b),
    enc_g     = collect.(collect(enc_g)),
    enc_b     = collect.(collect(enc_b)),
    dec_g     = collect.(collect(dec_g)),
    dec_b     = collect.(collect(dec_b)),
    force_g   = collect.(collect(force_g)),
    force_b   = collect.(collect(force_b)),
    mu        = collect(mu),
    omega     = collect(omega),
    K_real    = collect(K_real),
    Q         = collect(Q),
    z0        = collect(z0),
    u_rec     = collect(u_rec),
    case_idx  = case_idx,
    nstates   = nstates,
    nlatent   = nlatent,
    nx        = nx,
)

println("\nSaved AE-IC diagnostic tensors → $OUT_JLD2")
println("Load these in the WATT example and compare against your reimplementation's")
println("u_flat, enc_g/b, dec_g/b, z, u_decoded byte-for-byte to localize any divergence.")
