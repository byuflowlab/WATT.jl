#=
Tests for WATT's BEMT solver.

WATT.solve_BEMT is a thin wrapper around CCBlade.jl — the actual 
blade-element-momentum-theory math lives upstream. What can break
in the wrapper is the *translation* between WATT's data structures 
(Rotor, Blade, Environment,DS.Airfoil polars) and CCBlade's 
(Rotor, Section, OperatingPoint, AlphaAF).
That's what these tests cover: given identical physical inputs, WATT and a
direct CCBlade call must produce numerically identical Outputs.

Adam Cardoza
=#
using Test, WATT, DelimitedFiles, FLOWMath, CCBlade, OpenFASTTools, DynamicStallModels, StructArrays

of = OpenFASTTools
DS = DynamicStallModels

localpath = @__DIR__
cd(localpath)

errfun(x, xt) = 100*(x - xt)/xt

@testset "Blade Element Momentum Theory" begin

    @testset "CCBlade" begin

        ### -------------------------------------------------------------------
        ### Read NREL 5MW rotor geometry from OpenFAST input files
        ### -------------------------------------------------------------------
        ofpath = "../data/openfast"
        adblade = of.read_adblade("sn5_adblade.dat", ofpath)   # AeroDyn blade table (chord, twist, span, airfoil id)
        edfile  = of.read_edfile("sn5_EDfile.dat", ofpath)     # ElastoDyn file (hub radius, etc.)

        ### Read every airfoil polar referenced by the blade table
        aftypes = Array{of.AirfoilInput}(undef, 8)
        aftypes[1] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder1.dat"))
        aftypes[2] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "Cylinder2.dat"))
        aftypes[3] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU40_A17.dat"))
        aftypes[4] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU35_A17.dat"))
        aftypes[5] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU30_A17.dat"))
        aftypes[6] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU25_A17.dat"))
        aftypes[7] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "DU21_A17.dat"))
        aftypes[8] = of.read_airfoilinput(joinpath(ofpath, "airfoils", "NACA64_A17.dat"))

        # Map each blade station to its airfoil polar
        af_idx = Int.(adblade["BlAFID"])
        afs    = aftypes[af_idx]

        ### -------------------------------------------------------------------
        ### Extract blade geometry. Note units: BlTwist is in DEGREES in the
        ### OpenFAST file; WATT.Blade and CCBlade.Section both expect RADIANS.
        ### -------------------------------------------------------------------
        chordvec = adblade["BlChord"]
        twistvec = adblade["BlTwist"]               # degrees (from OpenFAST)
        rhub     = edfile["HubRad"]
        rvec     = adblade["BlSpn"] .+ rhub         # absolute radii: hub-relative span + hub radius
        rtip     = rvec[end]
        precone  = 0.0
        hubht    = 80.0
        n        = length(rvec)

        ### Build the DynamicStallModels airfoil objects. We store them as a
        ### StructArray (not a Vector) because DynamicStallModels overrides
        ### getproperty(::AbstractVector{<:Airfoil}, sym) to do field-wise
        ### broadcasts; that override assumes a StructArray layout and trips
        ### on UndefRef slots in a plain `Vector{Airfoil}(undef, n)`.
        airfoils = StructArray{DS.Airfoil}(undef, n)
        xcp = Vector{Float64}(undef, n)
        for i = 1:n
            airfoils[i], xcp[i] = of.make_dsairfoil(afs[i])
        end

        # Twist must be in radians here.
        blade = Blade(rvec, chordvec, twistvec.*(pi/180), xcp, airfoils)

        ### -------------------------------------------------------------------
        ### Rotor + operating point (NREL 5MW at rated wind speed, TSR=7.55)
        ### -------------------------------------------------------------------
        B       = 3
        turbine = true
        rotor   = WATT.Rotor(B, hubht, turbine)

        vinf   = 10.0
        tsr    = 7.55
        rotorR = rtip*cos(precone)
        omega  = vinf*tsr/rotorR

        rho      = 1.225
        mu       = 1.464e-5
        a        = 343.0
        shearexp = 0.0
        env      = environment(rho, mu, a, vinf, omega, shearexp)

        pitch = 0.0
        Vx    = vinf                # axial inflow at every section (no yaw, no shear)

        xv = zeros(11)   # scratch buffer reused inside solve_BEMT

        ### -------------------------------------------------------------------
        ### Compare WATT.solve_BEMT against a direct CCBlade.solve at each
        ### interior blade station. This is the wrapper-fidelity test: WATT's
        ### data-structure translation must yield the same numbers a direct
        ### CCBlade call would produce for the same physical inputs.
        ###
        ### Endpoints (idx=1 and idx=n) are deliberately excluded:
        ###   rvec[1] == rhub and rvec[end] == rtip (since BlSpn[1]==0 and
        ###   BlSpn[end]==blade length). CCBlade short-circuits any section
        ###   coincident with rhub/rtip to zero Outputs (CCBlade.jl:441),
        ###   while WATT.solve_BEMT computes nonzero loads there. The
        ###   simulation in practice never places aero nodes on the hub or
        ###   tip — see the checkforwarnings() guard in aerostructural.jl —
        ###   so the endpoint disagreement is outside both solvers' contract.
        ### -------------------------------------------------------------------
        @testset "BEM vs CCBlade (interior nodes)" begin
            for idx = 2:n-1
                Vy = omega*rvec[idx]

                # WATT path: pass the WATT-side structs through the wrapper.
                rotorout = WATT.solve_BEMT(rotor, blade, env, idx, Vx, Vy, pitch, xv; npts=10)

                # Direct CCBlade path: rebuild equivalent CCBlade inputs from
                # the *same* underlying polars and geometry, then call solve().
                #
                # Note on units: CCBlade.Section expects twist in RADIANS, so
                # we re-convert twistvec[idx] (which is still in degrees here
                # — we only converted the copy that went into the Blade above).
                # An earlier version of this test omitted the *(pi/180) and
                # CCBlade silently converged to a non-physical root (a≈0.99,
                # phi opposite sign), making every field comparison fail.
                polar    = blade.airfoils[idx].polar
                af       = CCBlade.AlphaAF(polar[:,1], polar[:,2], polar[:,3])
                section  = Section(rvec[idx], chordvec[idx], twistvec[idx]*(pi/180), af)
                rotor_cc = CCBlade.Rotor(blade.rhub, blade.rtip, B; turbine=true, tip=nothing)
                op       = CCBlade.OperatingPoint(Vx, Vy, rho, pitch, mu, a)
                ccout    = solve(rotor_cc, section, op)

                # Type fidelity: WATT must return CCBlade's Outputs struct
                # untouched (not wrap it in its own type).
                @test isa(rotorout, typeof(ccout))

                # Numerical fidelity: every field must agree. A field-named
                # inner @testset gives precise diagnostics on failure
                # (otherwise "Test Failed" alone wouldn't tell us idx or which
                # quantity disagrees).
                @testset "idx=$idx field=$f" for f in fieldnames(typeof(ccout))
                    @test isapprox(getfield(rotorout, f), getfield(ccout, f); rtol=1e-8, atol=1e-10)
                end

                # Physical sanity: a > 1 means the BEM has entered the
                # turbulent-wake regime where standard momentum theory breaks
                # down. If this fires for the NREL 5MW at rated, something
                # upstream (geometry, polars, op point) is broken.
                @test ccout.a < 1
            end
        end


    end #End test CCBlade
end #End test BEMT

nothing

