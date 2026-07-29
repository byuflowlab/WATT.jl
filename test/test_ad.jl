#=
AD compatibility tests, gated on ENV["WATT_AD_TESTS"]=="true".

Two layers:

1. BEM-level (cheap): ForwardDiff and ReverseDiff propagate through
   `solve_BEMT` via the ImplicitAD wrapping of the BEM residual. Verified
   against central finite differences to ~1e-10.

2. Coupled-solver level (expensive): ForwardDiff propagates through
   `run_sim!` provided the structural sensitivity parameters
   (per-element compliance and mass matrices) are packed into the
   parameter vector `p` and unpacked inside a `pfunc` that GXBeam calls
   on every step. Aero loads are written into `p` via a `prepp!`
   callback so that dual-typed aero outputs flow into the dual-typed
   distributed_loads. The non-dual `assembly` passed to
   `initialize_sim` is used only as a structural skeleton —
   `pfunc(p, t)` rebuilds the dual-typed Assembly internally on every
   GXBeam call. This is the pattern used in the Cardoza2026 fatigue
   optimization workflow.

Design vars for the gradient test: the full flattened compliance and
mass vector across every element (36·nelem·2 entries). We FD-spot-check
a handful of entries rather than the whole gradient.

Adam Cardoza
=#

using Test
using WATT, ForwardDiff, ReverseDiff, FiniteDifferences
using GXBeam, StaticArrays, StructArrays
using OpenFASTTools, DynamicStallModels, JLD2
using Random

const _of = OpenFASTTools
const _DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)


# ---------------------------------------------------------------------
# Inline NREL 5MW setup (the fixture used by the other test files is a
# plain Float64 skeleton; here we need the assembly to be reconstructible
# from a packed `p`, so we build everything from scratch).
# ---------------------------------------------------------------------
let ofpath = "../data/openfast"
    global blade_ad, assembly_ad, rotor_ad, env_ad, n_ad, nelem_ad

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
    twistvec = adblade["BlTwist"]
    rhub = edfile["HubRad"]
    rvec = adblade["BlSpn"] .+ rhub
    rtip = rvec[end]
    n_ad = length(rvec)

    airfoils = StructArray{_DS.Airfoil}(undef, n_ad)
    xcp = Vector{Float64}(undef, n_ad)
    for i = 1:n_ad
        airfoils[i], xcp[i] = _of.make_dsairfoil(afs[i])
    end

    blade_ad    = WATT.Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils; rhub, rtip)
    rotor_ad    = WATT.Rotor(3, 80.0, true)

    vinf = 10.0
    omega = vinf * 7.55 / rtip
    env_ad = WATT.environment(1.225, 1.464e-5, 343.0, vinf, omega, 0.0)

    assembly_ad = _of.make_assembly(edfile, bdfile, bdblade)
    nelem_ad    = length(assembly_ad.elements)
end


# Reference finite-difference engine (Richardson-extrapolated 5-point
# central stencil). Tolerance for AD-vs-FD comparison is kept wide (rtol
# ~1e-2) because FD on a coupled aeroelastic time-stepping solver is
# inherently noisy.
const FDM = FiniteDifferences.central_fdm(5, 1)
fd_deriv(f, x) = FDM(f, x)


# ---------------------------------------------------------------------
# Helpers for packing/unpacking the structural sensitivity parameter
# vector `p`. WATT's aero arrays have length `n_aero` (number of
# aerodynamic nodes), which generally differs from `nelem` (number of
# structural elements) — interpolation between the two meshes happens
# inside `pfunc` here. Layout:
#     p = [vec(C_1); ...; vec(C_nelem);     <-- 36·nelem compliance entries
#          vec(M_1); ...; vec(M_nelem);     <-- 36·nelem mass entries
#          Fx_1, ..., Fx_n_aero,            <-- aero-node Fx (placeholder; prepp! fills)
#          Fy_1, ..., Fy_n_aero]            <-- aero-node Fy
# ---------------------------------------------------------------------
const NCE = 36          # entries per element-compliance matrix
const NME = 36          # entries per element-mass matrix

n_comp_(nelem)         = NCE * nelem
n_mass_(nelem)         = NME * nelem
n_loads_(n_aero)       = 2 * n_aero
n_p_(nelem, n_aero)    = n_comp_(nelem) + n_mass_(nelem) + n_loads_(n_aero)

function pack_p(assembly, n_aero, ::Type{TF}=Float64) where {TF}
    nelem = length(assembly.elements)
    p = Vector{TF}(undef, n_p_(nelem, n_aero))
    for i in 1:nelem
        p[(i-1)*NCE+1 : i*NCE] = vec(Matrix(assembly.elements[i].compliance))
    end
    off = n_comp_(nelem)
    for i in 1:nelem
        p[off + (i-1)*NME+1 : off + i*NME] = vec(Matrix(assembly.elements[i].mass))
    end
    p[off + n_mass_(nelem) + 1 : end] .= 0.0      # load slot — prepp! fills
    return p
end

# pfunc unpacks compliance + mass to rebuild the assembly, and maps the
# n_aero-length aero loads to nelem-length element loads via a constant
# nearest-node assignment (Fx[i], Fy[i] for the i-th structural
# element). This is intentionally simpler than WATT's update_forces!
# interpolation — the AD test is about derivative correctness through
# the pipeline, not about reproducing the production aero→structure map.
function make_pfunc(assembly, n_aero)
    nelem  = length(assembly.elements)
    points = assembly.points
    np     = length(points)
    start  = 1:nelem
    stop   = 2:np
    xp     = [assembly.elements[i].x for i in 1:nelem]
    nc     = n_comp_(nelem)
    nm     = n_mass_(nelem)

    function pfunc(p, t)
        comp = [SMatrix{6,6}(reshape(p[(i-1)*NCE+1 : i*NCE], 6, 6)) for i in 1:nelem]
        ms   = [SMatrix{6,6}(reshape(p[nc + (i-1)*NME+1 : nc + i*NME], 6, 6)) for i in 1:nelem]
        Fx   = @view p[nc + nm + 1          : nc + nm + n_aero]
        Fy   = @view p[nc + nm + n_aero + 1 : nc + nm + 2*n_aero]

        a = GXBeam.Assembly(points, start, stop; compliance=comp, mass=ms, midpoints=xp)

        TF = eltype(p)
        z  = SVector{3,TF}(0.0, 0.0, 0.0)
        # Constant per-element load: fz = -Fx[i], fy = Fy[i] for i ≤ nelem.
        # (WATT's update_forces! follows the convention GXBeam fz = -Fx.)
        dl = Dict(i => begin
            L = a.elements[i].L
            fy = Fy[i]
            fz = -Fx[i]
            GXBeam.DistributedLoads{TF}(
                SVector{3,TF}(0.0, fy*L/2, fz*L/2),
                SVector{3,TF}(0.0, fy*L/2, fz*L/2),
                z, z, z, z, z, z)
        end for i in 1:nelem)
        return (; assembly=a, distributed_loads=dl)
    end
    return pfunc
end

# Aero→structure callback: copies the n_aero-length aero arrays Fx, Fy
# verbatim into the load slot of p. pfunc handles the n_aero→nelem map.
function make_prepp(nelem, n_aero)
    nc = n_comp_(nelem)
    nm = n_mass_(nelem)
    function prepp!(pp, Fx, Fy, mx)
        pp[nc + nm + 1          : nc + nm + n_aero]   .= Fx
        pp[nc + nm + n_aero + 1 : nc + nm + 2*n_aero] .= Fy
    end
    return prepp!
end


@testset "AD compatibility" begin

    # -----------------------------------------------------------------
    # BEM-level (cheap, no GXBeam) — verifies ImplicitAD wrapping of the
    # BEM residual is correct.
    # -----------------------------------------------------------------
    @testset "ForwardDiff through solve_BEMT" begin
        idx = 10
        Vx = 10.0
        Vy = (10.0 * 7.55 / blade_ad.rtip) * blade_ad.r[idx]

        # chord scale
        bem_phi_scale = function(s)
            o = one(s)
            bld = WATT.Blade(blade_ad.r, blade_ad.c .* s, blade_ad.twist .* o,
                             blade_ad.xcp, blade_ad.airfoils;
                             rhub=blade_ad.rhub, rtip=blade_ad.rtip)
            xv = Vector{typeof(s)}(undef, 11)
            WATT.solve_BEMT(rotor_ad, bld, env_ad, idx, Vx, Vy, 0.0, xv).phi
        end
        d_ad = ForwardDiff.derivative(bem_phi_scale, 1.0)
        d_fd = fd_deriv(bem_phi_scale, 1.0)
        @test isapprox(d_ad, d_fd; rtol=1e-4)

        # pitch
        bem_phi_pitch = function(p)
            xv = Vector{typeof(p)}(undef, 11)
            WATT.solve_BEMT(rotor_ad, blade_ad, env_ad, idx, Vx, Vy, p, xv).phi
        end
        d_ad = ForwardDiff.derivative(bem_phi_pitch, 0.0)
        d_fd = fd_deriv(bem_phi_pitch, 0.0)
        @test isapprox(d_ad, d_fd; rtol=1e-4)
    end


    @testset "ReverseDiff through solve_BEMT" begin
        idx = 10
        Vx = 10.0
        Vy = (10.0 * 7.55 / blade_ad.rtip) * blade_ad.r[idx]

        bem_phi_scale_vec = function(svec)
            s = svec[1]
            o = one(s)
            bld = WATT.Blade(blade_ad.r, blade_ad.c .* s, blade_ad.twist .* o,
                             blade_ad.xcp, blade_ad.airfoils;
                             rhub=blade_ad.rhub, rtip=blade_ad.rtip)
            xv = Vector{typeof(s)}(undef, 11)
            WATT.solve_BEMT(rotor_ad, bld, env_ad, idx, Vx, Vy, 0.0, xv).phi
        end
        g    = ReverseDiff.gradient(bem_phi_scale_vec, [1.0])
        d_fd = fd_deriv(s -> bem_phi_scale_vec([s]), 1.0)
        @test isapprox(g[1], d_fd; rtol=1e-4)
    end


    # -----------------------------------------------------------------
    # Coupled-solver level — ForwardDiff.gradient through run_sim! over
    # the per-element compliance + mass parameter vector.
    # -----------------------------------------------------------------
    @testset "ForwardDiff.gradient through run_sim! (compliance+mass)" begin
        nelem  = nelem_ad
        n_aero = n_ad
        tvec   = collect(range(0.0, 0.01, length=3))    # 3-step, ~3° azimuth

        p0     = pack_p(assembly_ad, n_aero, Float64)
        pfunc  = make_pfunc(assembly_ad, n_aero)
        prepp! = make_prepp(nelem, n_aero)

        # Skeleton assembly handed to initialize_sim — pfunc rebuilds the
        # real (dual-typed) assembly internally on every step.
        a_skel = pfunc(p0, 0.0).assembly

        # Output: integrated tip deflection magnitude across the time
        # window — a structural scalar that exercises the full chain.
        function tip_def_sum(x)
            TF = eltype(x)
            p  = convert(Vector{TF}, x)

            # Promote blade.c and blade.twist to TF so that WATT's
            # `find_inittype(blade.c[1], blade.twist[1])` sees the dual
            # type and allocates dual-typed aerostates + gxhistory. The
            # numeric values are unchanged; this is a type-promotion-only
            # rebuild.
            bld = WATT.Blade(blade_ad.r,
                             convert(Vector{TF}, blade_ad.c),
                             convert(Vector{TF}, blade_ad.twist),
                             blade_ad.xcp, blade_ad.airfoils;
                             rhub=blade_ad.rhub, rtip=blade_ad.rtip)

            aerostates, gxhistory, mesh = WATT.initialize_sim(
                bld, a_skel, tvec; verbose=false, pfunc=pfunc, p=p)
            WATT.run_sim!(rotor_ad, bld, mesh, env_ad, tvec, aerostates, gxhistory;
                          verbose=false, prepp=prepp!, p=p)
            s = zero(TF)
            for i in eachindex(gxhistory)
                s += gxhistory[i].elements[end].u[3]
            end
            return s
        end

        # Non-dual forward call sanity check.
        f0 = tip_def_sum(p0)
        @test isfinite(f0)

        # Full gradient via ForwardDiff (chunked to keep memory sane).
        cfg  = ForwardDiff.GradientConfig(tip_def_sum, p0, ForwardDiff.Chunk{8}())
        grad = ForwardDiff.gradient(tip_def_sum, p0, cfg)

        @test length(grad) == length(p0)
        @test all(isfinite, grad)

        # Directional derivative check: compare ⟨∇f, v⟩ (from AD) to
        # ∂f/∂v (from FD) along a few random directions over the
        # compliance + mass region of p (excluding load slots, which are
        # overwritten by prepp! and have no meaningful gradient). Per-
        # entry spot checks are flaky on this problem — compliance and
        # mass entries span ~10 orders of magnitude, and individual
        # perturbations of the (assumed-)symmetric matrices are
        # ill-conditioned. The directional approach integrates out those
        # pathologies into a single FD evaluation.
        nc = n_comp_(nelem)
        nm = n_mass_(nelem)
        struct_idxs = 1:(nc + nm)            # compliance + mass region only
        rng = MersenneTwister(1234)

        for trial in 1:3
            v = zeros(length(p0))
            v[struct_idxs] .= randn(rng, length(struct_idxs))
            # Scale v by p0's magnitude per entry so the FD step lives
            # on the same scale as each design var (compliance ~1e-10,
            # mass ~1e3).
            v[struct_idxs] .*= abs.(@view p0[struct_idxs])
            v ./= sqrt(sum(abs2, v))         # unit vector

            f1d  = t -> tip_def_sum(p0 .+ t .* v)
            d_fd = fd_deriv(f1d, 0.0)
            d_ad = sum(grad .* v)

            # Wide tolerance: FD over a 3-step coupled aero+structural
            # solve is noisy. rtol=1e-2 is a sanity bound, not precision.
            @test isapprox(d_ad, d_fd; rtol=1e-2, atol=1e-6)
        end
    end


    # -----------------------------------------------------------------
    # Frozen-start windowed sensitivity — ForwardDiff through
    # `initialize_from_state` + `run_from_state!`. Warm up rest→s in
    # Float64, freeze the snapshot at s, then differentiate the windowed
    # march s→s+k w.r.t. the compliance+mass vector p. The snapshot is
    # held constant (zero partials): both AD and FD re-initialize from the
    # same Float64 state each call, so they measure the same frozen-start
    # derivative.
    # -----------------------------------------------------------------
    @testset "ForwardDiff through run_from_state! (frozen-start window)" begin
        nelem    = nelem_ad
        n_aero   = n_ad
        nt       = 6
        tvec     = collect(range(0.0, 0.02, length=nt))
        s        = 3                       # snapshot step
        tvec_win = tvec[s:end]             # window [t_s … t_end]
        t_s      = tvec[s]

        p0     = pack_p(assembly_ad, n_aero, Float64)
        pfunc  = make_pfunc(assembly_ad, n_aero)
        prepp! = make_prepp(nelem, n_aero)
        a_skel = pfunc(p0, 0.0).assembly

        # Warm up rest→s in Float64; freeze the native snapshot at s.
        aw, gw, mw = WATT.initialize_sim(blade_ad, a_skel, tvec; verbose=false, pfunc=pfunc, p=p0)
        WATT.run_sim!(rotor_ad, blade_ad, mw, env_ad, tvec, aw, gw; verbose=false, prepp=prepp!, p=p0)
        state_s   = gw[s]
        xds_s     = copy(aw.xds[s, :])
        azimuth_s = aw.azimuth[s]

        # Frozen-start windowed map: p ↦ sum of tip out-of-plane deflection over the window.
        function window_tip_sum(x)
            TF  = eltype(x)
            p   = convert(Vector{TF}, x)
            bld = WATT.Blade(blade_ad.r,
                             convert(Vector{TF}, blade_ad.c),
                             convert(Vector{TF}, blade_ad.twist),
                             blade_ad.xcp, blade_ad.airfoils;
                             rhub=blade_ad.rhub, rtip=blade_ad.rtip)
            # Dual-typed window mesh (buffers follow eltype(bld) = TF).
            _, _, meshw = WATT.initialize_sim(bld, a_skel, tvec_win; verbose=false, pfunc=pfunc, p=p)
            init = WATT.initialize_from_state(state_s, xds_s, azimuth_s, meshw, bld,
                                              env_ad, tvec_win, t_s; p=p)
            acc = Ref(zero(TF))
            WATT.run_from_state!(init..., meshw, rotor_ad, bld, env_ad, tvec_win;
                                 prepp=prepp!, p=p,
                                 out=(st, j) -> (acc[] += st.elements[end].u[3]))
            return acc[]
        end

        f0 = window_tip_sum(p0)
        @test isfinite(f0)

        cfg  = ForwardDiff.GradientConfig(window_tip_sum, p0, ForwardDiff.Chunk{8}())
        grad = ForwardDiff.gradient(window_tip_sum, p0, cfg)
        @test length(grad) == length(p0)
        @test all(isfinite, grad)

        # Directional AD-vs-FD over the compliance+mass region (frozen start on both sides).
        nc = n_comp_(nelem)
        nm = n_mass_(nelem)
        struct_idxs = 1:(nc + nm)
        rng = MersenneTwister(2468)
        for trial in 1:2
            v = zeros(length(p0))
            v[struct_idxs] .= randn(rng, length(struct_idxs))
            v[struct_idxs] .*= abs.(@view p0[struct_idxs])
            v ./= sqrt(sum(abs2, v))
            d_fd = fd_deriv(t -> window_tip_sum(p0 .+ t .* v), 0.0)
            d_ad = sum(grad .* v)
            @test isapprox(d_ad, d_fd; rtol=1e-2, atol=1e-6)
        end
    end

end #End AD compatibility

nothing
