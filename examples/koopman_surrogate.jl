#=
Shared NeuralKoopman structural-surrogate implementations for the WATT
surrogate-driven aeroelastic solvers.

Defines two concrete `AbstractStructuralSurrogate`s that load the *same* trained
Lux NeuralKoopman checkpoint (a plain NamedTuple of matrices, Lux-free):

  * `ConditionedKoopman` — host / single-sim. Byte-compatible with the
    hand-written forward pass in aerostructural_nrel5mw5seg_surrogate.jl.
    Implements `encode_initial` / `step_latent` / `decode` (returns a
    `SurrogateAssemblyState`) for `run_sim_surrogate!`.

  * `BatchedKoopman` — device / batched over `ns` sims, each with its own FiLM
    conditioning vector `x`. Re-expresses the per-component Hermite regrids as
    constant matrices (so regridding is a matmul) and carries FiLM γ/β as
    `(width, ns)` arrays; the Koopman math is otherwise the same matrix ops.
    Implements `encode_initial(surr, ::GPUPointStates)` / `step_latent` /
    `decode!(surr, z, ::GPUPointStates)` for `run_sim_surrogate_gpu!`.

Both share the primitives below; the array eltype/backend is chosen by the
`ArrayType` passed to the batched builder, so the batched model runs unchanged
on the CPU (`Array`) backend (used for validation) or a GPU array type.

Adam Cardoza
=#

using WATT, GXBeam, StaticArrays, JLD2, LinearAlgebra

# ---------------------------------------------------------------------------
# Forward primitives (Lux-free), shared by both models.
# ---------------------------------------------------------------------------
mish(x) = x * tanh(log1p(exp(x)))
sigmoid(x) = one(x) / (one(x) + exp(-x))

f64(x::AbstractArray{<:AbstractFloat}) = Float64.(x)
f64(x::AbstractArray) = x
f64(x::NamedTuple) = NamedTuple{keys(x)}(map(f64, values(x)))
f64(x::Tuple) = map(f64, x)
f64(x::Real) = Float64(x)
f64(x) = x

# Recursive cast of a (possibly nested) NamedTuple parameter tree onto ArrayType.
to_dev(AT, x::AbstractArray) = AT(x)
to_dev(AT, x::NamedTuple) = NamedTuple{keys(x)}(map(v -> to_dev(AT, v), values(x)))
to_dev(AT, x::Tuple) = map(v -> to_dev(AT, v), x)
to_dev(AT, x) = x

# Hermite cubic regrid: y[s_src] → y[s_dst], scalar series (host, for CPU model
# and for building the constant regrid matrices below).
function hermite_regrid(y::AbstractVector{T}, s_src::AbstractVector,
                        s_dst::AbstractVector) where {T}
    n = length(s_src)
    ds = s_src[2:end] .- s_src[1:end-1]
    m_pair = (y[2:end] .- y[1:end-1]) ./ ds
    m = vcat(m_pair[1:1], (m_pair[2:end] .+ m_pair[1:end-1]) ./ T(2), m_pair[end:end])
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

# hermite_regrid is linear in y, so it equals M*y for a constant matrix M.
# Build M (length(s_dst) × length(s_src)) by regridding the unit vectors.
function build_regrid_matrix(s_src::AbstractVector, s_dst::AbstractVector)
    n = length(s_src)
    M = Matrix{Float64}(undef, length(s_dst), n)
    e = zeros(Float64, n)
    for j in 1:n
        fill!(e, 0.0); e[j] = 1.0
        M[:, j] = hermite_regrid(e, s_src, s_dst)
    end
    return M
end

# FiLM-conditioned MLP forward (hidden Denses + out Dense, FiLM on every layer).
# Works for a single column (ns=1) or batched (in, ns); gammas/betas are
# (width, ns) (or (width,1) broadcast). Reused verbatim by both models.
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

# FiLM hypernet forward + split into per-target-layer (gamma, beta) tuples.
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
        offset = sum(target_widths[1:i-1]); w = target_widths[i]
        1.0 .+ raw_g[offset+1:offset+w, :]
    end
    betas = ntuple(length(target_widths)) do i
        offset = sum(target_widths[1:i-1]); w = target_widths[i]
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
    r      = sigmoid.(g_r  .* ps_kq_global.r_raw      .+ b_r) .* c
    theta_b= g_th .* ps_kq_global.theta      .+ b_th
    K_real = tanh.(g_kr .* ps_kq_global.K_real_raw .+ b_kr) .* c
    Q_raw_b= g_q  .* ps_kq_global.Q_raw      .+ b_q
    mu     = r .* cos.(theta_b)
    omega  = r .* sin.(theta_b)
    r_pair = reshape(r, 1, size(r, 1), size(r, 2))
    eig_mag_complex = reshape(vcat(r_pair, r_pair), 2 * size(r, 1), size(r, 2))
    eig_mag = vcat(eig_mag_complex, abs.(K_real))
    q_max = sqrt.(max.(c^2 .- eig_mag .^ 2, 0.0))
    Q     = tanh.(Q_raw_b) .* q_max
    return mu, omega, K_real, Q
end

function apply_jordan_A(z, mu, omega, K_real, M::Int, ncp::Int)
    B = size(z, 2)
    if ncp > 0
        z1 = z[1:2:M, :]; z2 = z[2:2:M, :]
        z1n = mu .* z1 .- omega .* z2
        z2n = omega .* z1 .+ mu .* z2
        stacked = vcat(reshape(z1n, 1, ncp, B), reshape(z2n, 1, ncp, B))
        Az_c = reshape(stacked, M, B)
    else
        Az_c = similar(z, 0, B)
    end
    if size(K_real, 1) > 0
        return vcat(Az_c, K_real .* z[M+1:end, :])
    else
        return Az_c
    end
end

# ---------------------------------------------------------------------------
# Normalization constants (must match the training script).
# ---------------------------------------------------------------------------
const NORM_SCALES_VEC = Float64[
    1e+1, 1e+0, 1e+0,  1e+1, 1e+1, 1e+1,  1e-1, 1e-1, 1e+0,
    1e+0, 1e+0, 1e+0,  1e-6, 1e-6, 1e-6,  1e-5, 1e-7, 1e-7]
const FNORM    = 1e4
const N_REGRID = 6
const N_SEL    = 18
const FLIP_V_OMEGA_AT_BOUNDARY = false

# Grids shared by both models, derived from the assembly.
function koopman_grids(assembly::GXBeam.Assembly)
    nelem = length(assembly.elements)
    elem_r = Float64[norm(assembly.elements[j].x) for j in 1:nelem]
    s_elem = (elem_r .- elem_r[1]) ./ (elem_r[end] - elem_r[1])
    s_regrid = collect(range(0.0, 1.0, length=N_REGRID))
    np = length(assembly.points)
    s_state = collect(range(0.0, 1.0, length=np))
    return s_elem, s_regrid, s_state
end

# ---------------------------------------------------------------------------
# ConditionedKoopman — host / single-sim (SurrogateAssemblyState boundary).
# ---------------------------------------------------------------------------
struct ConditionedKoopman{TG, TB, TPE, TPD, TPF} <: WATT.AbstractStructuralSurrogate
    enc_g::TG;   enc_b::TB
    dec_g::TG;   dec_b::TB
    force_g::TG; force_b::TB
    mu::Matrix{Float64}; omega::Matrix{Float64}; K_real::Matrix{Float64}; Q::Matrix{Float64}
    ps_encoder::TPE; ps_decoder::TPD; ps_force::TPF
    n_complex_pairs::Int; n_complex_dims::Int; nlatent::Int
    s_elem_infer::Vector{Float64}; s_regrid::Vector{Float64}; s_state_infer::Vector{Float64}
end

function build_conditioned_koopman(jld_path::AbstractString, x_norm::AbstractVector,
                                   assembly::GXBeam.Assembly)
    ckpt = JLD2.load(jld_path)
    ps = f64(ckpt["ps"])
    nlatent = ckpt["nlatent"]; ncp = ckpt["n_complex_pairs"]
    nstates = ckpt["nstates"]; nforce_regrid = ckpt["nforce_regrid"]
    svd_buffer = Float64(ckpt["svd_buffer"])
    @assert nstates == N_SEL * N_REGRID
    @assert nforce_regrid == 6 * N_REGRID
    M = 2 * ncp; n_real = nlatent - M
    ae_widths = ckpt["ae_widths"]; force_widths = ckpt["force_widths"]
    enc_widths  = vcat(collect(ae_widths), nlatent)
    dec_widths  = vcat(reverse(collect(ae_widths)), nstates)
    force_widths_full = vcat(collect(force_widths), nlatent)
    x_b = reshape(Float64.(x_norm), :, 1)
    enc_g, enc_b     = hyper_forward(x_b, ps.enc_hyper, enc_widths)
    dec_g, dec_b     = hyper_forward(x_b, ps.dec_hyper, dec_widths)
    force_g, force_b = hyper_forward(x_b, ps.force_hyper, force_widths_full)
    mu, omega, K_real, Q = get_eigs(ps.kq_hyper, ps.kq_global, x_b, ncp, n_real, nlatent, svd_buffer)
    s_elem, s_regrid, s_state = koopman_grids(assembly)
    return ConditionedKoopman(enc_g, enc_b, dec_g, dec_b, force_g, force_b,
        mu, omega, K_real, Q, ps.encoder, ps.decoder, ps.force_nn,
        ncp, M, nlatent, s_elem, s_regrid, s_state)
end

function WATT.encode_initial(surr::ConditionedKoopman, u0::WATT.SurrogateAssemblyState)
    np = length(u0.points)
    u_raw = Matrix{Float64}(undef, np, N_SEL)
    s = FLIP_V_OMEGA_AT_BOUNDARY ? -1.0 : 1.0
    for j in 1:np
        p = u0.points[j]
        u_raw[j, 1:3]=p.u; u_raw[j,4:6]=p.theta; u_raw[j,7:9]=s.*p.V
        u_raw[j,10:12]=s.*p.Omega; u_raw[j,13:15]=p.F; u_raw[j,16:18]=p.M
    end
    u_regrid = Matrix{Float64}(undef, N_REGRID, N_SEL)
    for k in 1:N_SEL
        u_regrid[:, k] = hermite_regrid(u_raw[:, k], surr.s_state_infer, surr.s_regrid) .* NORM_SCALES_VEC[k]
    end
    u_flat = reshape(u_regrid, N_REGRID * N_SEL, 1)
    return film_mlp_forward(u_flat, surr.ps_encoder, surr.enc_g, surr.enc_b, mish)
end

function WATT.step_latent(surr::ConditionedKoopman, z, f_per_element::AbstractMatrix)
    f_regrid = Matrix{Float64}(undef, N_REGRID, 6)
    for k in 1:6
        f_regrid[:, k] = hermite_regrid(Float64.(f_per_element[:, k]) ./ FNORM, surr.s_elem_infer, surr.s_regrid)
    end
    f_flat = reshape(f_regrid, 6 * N_REGRID, 1)
    fi_enc = film_mlp_forward(f_flat, surr.ps_force, surr.force_g, surr.force_b, mish)
    return apply_jordan_A(z, surr.mu, surr.omega, surr.K_real, surr.n_complex_dims, surr.n_complex_pairs) .+ surr.Q .* fi_enc
end

function WATT.decode(surr::ConditionedKoopman, z)
    u_flat = film_mlp_forward(z, surr.ps_decoder, surr.dec_g, surr.dec_b, mish)
    u_scaled = reshape(u_flat, N_REGRID, N_SEL)
    np = length(surr.s_state_infer)
    u_pts = Matrix{Float64}(undef, np, N_SEL)
    for k in 1:N_SEL
        u_pts[:, k] = hermite_regrid(u_scaled[:, k] ./ NORM_SCALES_VEC[k], surr.s_regrid, surr.s_state_infer)
    end
    s = FLIP_V_OMEGA_AT_BOUNDARY ? -1.0 : 1.0
    points = Vector{WATT.SurrogatePointState{Float64}}(undef, np)
    for j in 1:np
        points[j] = WATT.SurrogatePointState{Float64}(
            SVector{3,Float64}(u_pts[j,1],u_pts[j,2],u_pts[j,3]),
            SVector{3,Float64}(u_pts[j,4],u_pts[j,5],u_pts[j,6]),
            SVector{3,Float64}(s*u_pts[j,7],s*u_pts[j,8],s*u_pts[j,9]),
            SVector{3,Float64}(s*u_pts[j,10],s*u_pts[j,11],s*u_pts[j,12]),
            SVector{3,Float64}(u_pts[j,13],u_pts[j,14],u_pts[j,15]),
            SVector{3,Float64}(u_pts[j,16],u_pts[j,17],u_pts[j,18]))
    end
    return WATT.SurrogateAssemblyState{Float64}(points)
end

# ---------------------------------------------------------------------------
# BatchedKoopman — device / batched over ns (each sim its own conditioning x).
# ---------------------------------------------------------------------------
struct BatchedKoopman{TM, TG, TP, TV} <: WATT.AbstractStructuralSurrogate
    M_enc::TM      # (N_REGRID × np)   state regrid np→6
    M_force::TM    # (N_REGRID × nelem) force regrid nelem→6
    M_dec::TM      # (np × N_REGRID)   state regrid 6→np
    enc_g::TG;   enc_b::TG
    dec_g::TG;   dec_b::TG
    force_g::TG; force_b::TG
    mu::TM; omega::TM; K_real::TM; Q::TM         # (·, ns)
    ps_encoder::TP; ps_decoder::TP; ps_force::TP
    scale::TV      # (N_SEL,) NORM_SCALES on device
    ncp::Int; ncd::Int; nlatent::Int; sflip::Float64
end

"""
    build_batched_koopman(jld_path, x_b, assembly; ArrayType=Array{Float64})

Build a [`BatchedKoopman`](@ref) from the trained checkpoint. `x_b` is the
`(nx, ns)` FiLM conditioning matrix — one column per sim. All parameters,
regrid matrices, and FiLM γ/β land on `ArrayType`.
"""
function build_batched_koopman(jld_path::AbstractString, x_b::AbstractMatrix,
                               assembly::GXBeam.Assembly; ArrayType::Type=Array{Float64})
    ckpt = JLD2.load(jld_path)
    ps = f64(ckpt["ps"])
    nlatent = ckpt["nlatent"]; ncp = ckpt["n_complex_pairs"]
    nstates = ckpt["nstates"]; svd_buffer = Float64(ckpt["svd_buffer"])
    M = 2 * ncp; n_real = nlatent - M
    ae_widths = ckpt["ae_widths"]; force_widths = ckpt["force_widths"]
    enc_widths  = vcat(collect(ae_widths), nlatent)
    dec_widths  = vcat(reverse(collect(ae_widths)), nstates)
    force_widths_full = vcat(collect(force_widths), nlatent)
    xb = Float64.(x_b)
    enc_g, enc_b     = hyper_forward(xb, ps.enc_hyper, enc_widths)
    dec_g, dec_b     = hyper_forward(xb, ps.dec_hyper, dec_widths)
    force_g, force_b = hyper_forward(xb, ps.force_hyper, force_widths_full)
    mu, omega, K_real, Q = get_eigs(ps.kq_hyper, ps.kq_global, xb, ncp, n_real, nlatent, svd_buffer)

    s_elem, s_regrid, s_state = koopman_grids(assembly)
    M_enc   = build_regrid_matrix(s_state, s_regrid)   # (6 × np)
    M_force = build_regrid_matrix(s_elem,  s_regrid)   # (6 × nelem)
    M_dec   = build_regrid_matrix(s_regrid, s_state)   # (np × 6)

    dev(x) = ArrayType(Float64.(x))
    devt(t) = map(dev, t)
    sflip = FLIP_V_OMEGA_AT_BOUNDARY ? -1.0 : 1.0
    return BatchedKoopman(
        dev(M_enc), dev(M_force), dev(M_dec),
        devt(enc_g), devt(enc_b), devt(dec_g), devt(dec_b), devt(force_g), devt(force_b),
        dev(mu), dev(omega), dev(K_real), dev(Q),
        to_dev(ArrayType, ps.encoder), to_dev(ArrayType, ps.decoder), to_dev(ArrayType, ps.force_nn),
        dev(NORM_SCALES_VEC), ncp, M, nlatent, sflip)
end

function WATT.encode_initial(surr::BatchedKoopman, ps::GPUPointStates)
    np = size(ps.u, 2); ns = size(ps.u, 3)
    s = surr.sflip
    pt(A) = permutedims(A, (2, 1, 3))                        # (3,np,ns)→(np,3,ns)
    u_raw = cat(pt(ps.u), pt(ps.theta), s .* pt(ps.V), s .* pt(ps.Omega),
                pt(ps.F), pt(ps.M); dims=2)                  # (np,18,ns)
    u_regrid = reshape(surr.M_enc * reshape(u_raw, np, N_SEL * ns), N_REGRID, N_SEL, ns)
    u_regrid = u_regrid .* reshape(surr.scale, 1, N_SEL, 1)
    u_flat = reshape(u_regrid, N_REGRID * N_SEL, ns)
    return film_mlp_forward(u_flat, surr.ps_encoder, surr.enc_g, surr.enc_b, mish)
end

function WATT.step_latent(surr::BatchedKoopman, z, f_elem::AbstractArray)
    nelem = size(f_elem, 1); ns = size(f_elem, 3)
    f_regrid = reshape(surr.M_force * reshape(f_elem, nelem, 6 * ns) ./ FNORM, N_REGRID, 6, ns)
    f_flat = reshape(f_regrid, 6 * N_REGRID, ns)
    fi_enc = film_mlp_forward(f_flat, surr.ps_force, surr.force_g, surr.force_b, mish)
    return apply_jordan_A(z, surr.mu, surr.omega, surr.K_real, surr.ncd, surr.ncp) .+ surr.Q .* fi_enc
end

"""
    HostBridgeKoopman <: AbstractStructuralSurrogate

Single-sim wrapper that runs the Koopman math on `ArrayType` (a device) but
presents the *host* surrogate interface, so the CPU coupled solver
[`run_sim_surrogate!`](@ref) drives it with CPU BEMT + CPU dynamic stall while
the structural step runs on the GPU. Each step uploads the `(nelem,6)` load
matrix, keeps the latent state on device, and downloads the decoded state to a
host `SurrogateAssemblyState`. Use to benchmark "original aero on CPU +
structural surrogate on GPU". Wraps a `ns=1` [`BatchedKoopman`](@ref).
"""
struct HostBridgeKoopman{BK, GPS} <: WATT.AbstractStructuralSurrogate
    bk::BK
    gps_in::GPS
    gps_out::GPS
    AT::Type
    nelem::Int
end

function build_hostbridge_koopman(jld_path::AbstractString, x_norm::AbstractVector,
                                  assembly::GXBeam.Assembly; ArrayType::Type=Array{Float64})
    bk = build_batched_koopman(jld_path, reshape(Float64.(x_norm), :, 1), assembly; ArrayType=ArrayType)
    np = length(assembly.points); nelem = length(assembly.elements)
    return HostBridgeKoopman(bk, GPUPointStates(np, 1; ArrayType=ArrayType),
                             GPUPointStates(np, 1; ArrayType=ArrayType), ArrayType, nelem)
end

function _upload_state!(gps::GPUPointStates, u0::WATT.SurrogateAssemblyState)
    np = length(u0.points)
    h = (u=zeros(3,np,1), theta=zeros(3,np,1), V=zeros(3,np,1),
         Omega=zeros(3,np,1), F=zeros(3,np,1), M=zeros(3,np,1))
    for j in 1:np
        p = u0.points[j]
        h.u[:,j,1]=p.u; h.theta[:,j,1]=p.theta; h.V[:,j,1]=p.V
        h.Omega[:,j,1]=p.Omega; h.F[:,j,1]=p.F; h.M[:,j,1]=p.M
    end
    copyto!(gps.u,h.u); copyto!(gps.theta,h.theta); copyto!(gps.V,h.V)
    copyto!(gps.Omega,h.Omega); copyto!(gps.F,h.F); copyto!(gps.M,h.M)
    return gps
end

function _download_state(gps::GPUPointStates)
    u=Array(gps.u); th=Array(gps.theta); V=Array(gps.V)
    Om=Array(gps.Omega); F=Array(gps.F); M=Array(gps.M)
    np = size(u, 2)
    pts = [WATT.SurrogatePointState{Float64}(
              SVector{3,Float64}(u[:,j,1]),  SVector{3,Float64}(th[:,j,1]),
              SVector{3,Float64}(V[:,j,1]),  SVector{3,Float64}(Om[:,j,1]),
              SVector{3,Float64}(F[:,j,1]),  SVector{3,Float64}(M[:,j,1])) for j in 1:np]
    return WATT.SurrogateAssemblyState{Float64}(pts)
end

function WATT.encode_initial(s::HostBridgeKoopman, u0::WATT.SurrogateAssemblyState)
    _upload_state!(s.gps_in, u0)
    return WATT.encode_initial(s.bk, s.gps_in)     # z stays on device
end

function WATT.step_latent(s::HostBridgeKoopman, z, f::AbstractMatrix)
    f_dev = s.AT(reshape(Float64.(f), s.nelem, 6, 1))
    return WATT.step_latent(s.bk, z, f_dev)
end

function WATT.decode(s::HostBridgeKoopman, z)
    WATT.decode!(s.bk, z, s.gps_out)
    return _download_state(s.gps_out)
end

function WATT.decode!(surr::BatchedKoopman, z, ps::GPUPointStates)
    np = size(ps.u, 2); ns = size(z, 2); s = surr.sflip
    u_flat = film_mlp_forward(z, surr.ps_decoder, surr.dec_g, surr.dec_b, mish)  # (108, ns)
    u_scaled = reshape(u_flat, N_REGRID, N_SEL, ns) ./ reshape(surr.scale, 1, N_SEL, 1)
    u_pts = reshape(surr.M_dec * reshape(u_scaled, N_REGRID, N_SEL * ns), np, N_SEL, ns)  # (np,18,ns)
    ip(A) = permutedims(A, (2, 1, 3))                        # (np,3,ns)→(3,np,ns)
    ps.u     .= ip(view(u_pts, :, 1:3, :))
    ps.theta .= ip(view(u_pts, :, 4:6, :))
    ps.V     .= s .* ip(view(u_pts, :, 7:9, :))
    ps.Omega .= s .* ip(view(u_pts, :, 10:12, :))
    ps.F     .= ip(view(u_pts, :, 13:15, :))
    ps.M     .= ip(view(u_pts, :, 16:18, :))
    return ps
end
