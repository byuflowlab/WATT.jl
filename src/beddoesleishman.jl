#=
Code to interact with the DynamicStallModels Package, specifically for the Risø model. 

=#



function initializeDSmodel(dsmodel::DS.BeddoesLeishman, dsmodelinit::ModelInit, solver::Solver, turbine::Bool, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

    if isa(dsmodel.detype, DS.Functional)
        error("The Beddoes-Leishman functional form isn't ready yet. ")
    elseif isa(dsmodel.detype, DS.Iterative)
        error("The Beddoes-Leishman iterative implementation is still tying its' shoes.")
    else #Model is an indicial
        if dsmodel.version==1 #Original implementation
            error("Original Beddoes-Leishman model not yet ready.")

        elseif dsmodel.version==2 #AeroDyn Original
            return initializeBLA(dsmodel, dsmodelinit, solver, turbine, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

        elseif dsmodel.version==3 #AeroDyn Gonzalez
            # error("The AreoDyn Gonzalez implementation isn't ready yet.")
            # println("Got to initializing Gonzalez model") #Got here, but I can't figure out where it goes. 
            return initializeBLAG(dsmodel, dsmodelinit, solver, turbine, nt, na, tvec, Vxvec, Vxdotvec, chordvec, twistvec, phivec, pitch, a) 

        elseif dsmodel.version==4 #AeroDyn Minema
            error("The AeroDyn Minema implementation needs work. Try a different model.")
        end
    end
    error("Dynamic stall model failed to initialize. You likely chose a model that is not yet incorporated.")
end




function update_aero_parameters!(dsmodel::DS.BeddoesLeishman, turbine::Bool, pds, na, rvec, Wvec, phivec, twistvec, pitch, env, t)
    if isa(dsmodel.detype, DS.Functional)
        error("The Beddoes-Leishman functional form isn't ready yet. ")
    elseif isa(dsmodel.detype, DS.Iterative)
        error("The Beddoes-Leishman iterative implementation is still tying its' shoes.")
    else #Model is an indicial
        if dsmodel.version==1 #Original implementation
            error("Original Beddoes-Leishman model not yet ready.")
        elseif dsmodel.version==2 #AeroDyn Original
            return update_BLA_parameters!(dsmodel, turbine, pds, na, rvec, Wvec, phivec, twistvec, pitch, env, t)
        elseif dsmodel.version==3 #AeroDyn Gonzalez
            return update_BLAG_parameters!(dsmodel, turbine, pds, na, rvec, Wvec, phivec, twistvec, pitch, env, t)
            # error("The AreoDyn Gonzalez implementation isn't ready yet.")
        elseif dsmodel.version==4 #AeroDyn Minema
            error("The AeroDyn Minema implementation needs work. Try a different model.")
        end
    end
    error("Dynamic stall model failed to initialize. You likely chose a model that is not yet incorporated.")
end


function extractloads(dsmodel::DS.BeddoesLeishman, x, ccout, t, rvec, chordvec, twistvec, pitch, blade::Blade, env::Environment) #TODO: This should probably be an inplace function. 
    if isa(dsmodel.detype, DS.Functional)
        error("The Beddoes-Leishman functional form isn't ready yet. ")
    elseif isa(dsmodel.detype, DS.Iterative)
        error("The Beddoes-Leishman iterative implementation is still tying its' shoes.")
    else #Model is an indicial
        if dsmodel.version==1 #Original implementation
            error("Original Beddoes-Leishman model not yet ready.")
        elseif dsmodel.version==2 #AeroDyn Original
            return extractloads_BLA(dsmodel, x, ccout, t, rvec, chordvec, twistvec, pitch, blade, env)
        elseif dsmodel.version==3 #AeroDyn Gonzalez
            # error("The AreoDyn Gonzalez implementation isn't ready yet.")
            return extractloads_BLAG(dsmodel, x, ccout, t, rvec, chordvec, twistvec, pitch, blade, env)
        elseif dsmodel.version==4 #AeroDyn Minema
            error("The AeroDyn Minema implementation needs work. Try a different model.")
        end
    end
    error("Dynamic stall model failed to initialize. You likely chose a model that is not yet incorporated.")
end


function parsesolution(dsmodel::DS.BeddoesLeishman, xds, W, phi, tvec, chordvec, twistvec, pitch, blade::Blade) #Todo: What is the purpose of this function? 
    if isa(dsmodel.detype, DS.Functional)
        error("The Beddoes-Leishman functional form isn't ready yet. ")
    elseif isa(dsmodel.detype, DS.Iterative)
        error("The Beddoes-Leishman iterative implementation is still tying its' shoes.")
    else #Model is an indicial
        if dsmodel.version==1 #Original implementation
            error("Original Beddoes-Leishman model not yet ready.")
        elseif dsmodel.version==2 #AeroDyn Original
            return intializeBLA()
        elseif dsmodel.version==3 #AeroDyn Gonzalez
            error("The AreoDyn Gonzalez implementation isn't ready yet.")
        elseif dsmodel.version==4 #AeroDyn Minema
            error("The AeroDyn Minema implementation needs work. Try a different model.")
        end
    end
    error("Dynamic stall model failed to initialize. You likely chose a model that is not yet incorporated.")
end
