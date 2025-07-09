include("../simulation/basic_sim.jl")
include("../utils/paths.jl")

using Dates


const ANALYSIS_ROOT_DIR = "work/analysis/"


function getanalysisdir(correlationfunction::String, simparams::SimulationParameters, forcenewdir::Bool=false)::String
    # get abspath of analysis root directory
    analysisrootdir::String = relpath2abspath(ANALYSIS_ROOT_DIR)

    # make correlation function directory if it doesn't already exist
    correlationfnctdir::String = joinpath(analysisrootdir, correlationfunction)
    if !isdir(correlationfnctdir)
        mkdir(correlationfnctdir)
    end
    

    # make interaction type directory if it doesn't already exist
    local interactionstr::String = ""
    if simparams.interaction == nointeraction
        interactionstr = "No Interaction"
    elseif simparams.interaction == alignsimple
        interactionstr = "Align Simple"
    elseif simparams.interaction == antialignsimple
        interactionstr = "Antialign Simple"
    end
   
    interactiondir::String = joinpath(correlationfnctdir, interactionstr)
    if !isdir(interactiondir)
        mkdir(interactiondir)
    end

    # make simparam-specific directory if it doesn't exist already
    local spdirname::String
    if simparams.interaction == nointeraction
        spdirname = joinpath(interactiondir, "N$(simparams.numparticles)-B$(Int64(round(simparams.boxwidth)))-T$(Int64(round(simparams.totaltime)))")
    else
        spdirname = joinpath(interactiondir, "N$(simparams.numparticles)-B$(Int64(round(simparams.boxwidth)))-T$(Int64(round(simparams.totaltime)))-I$(Int64(round(simparams.interactionfliprate)))")
    end

    spdirname_try = spdirname
    i::Int = 0
    # make new dirname if want to make a new one
    if forcenewdir
        while isdir(spdirname_try)
            spdirname_try = "$(spdirname_try)_$(i)"
        end

        spdirname = spdirname_try
    end

    # make dir if it doesn't already exist
    if !isdir(spdirname)
        mkdir(spdirname)
    end

    # return dir name
    return spdirname
end


