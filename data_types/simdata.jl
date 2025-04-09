include("../particle_types/particle.jl")

"""
Simulation Data Structure:

Holds particle and simulation params,
and position data of each particle


"""
struct ParticlePath
    name::String
    times::Array{Real}
    positions::Array{Real}
end

function start_t(path::ParticlePath)
    return path.times[1]
end

function end_t(path::ParticlePath)
    return path.times[end]
end


mutable struct ParticleSimData
    particlepaths::Array{ParticlePath}
    totaltime::Real
    dt::Real
    boxwidth::Real
end

function npaths(simdata::ParticleSimData)
    return length(simdata.particlepaths)
end

function periodicboundary_paths(simdata::ParticleSimData)
    # dictionary from index of particle path
    # to exttra offset needed for plot
    addedlines = Dict{Int, Tuple{Int, Int}}()

    for i in 1:length(simdata.particlepaths)
        currpath = simdata.particlepaths[i].positions

        # check how much the path goes over the positive edge
        overboundary = (maximum(currpath) - (simdata.boxwidth / 2)) / simdata.boxwidth
        if overboundary < 0
            overboundary = 0
        end

        # check how much the path goes under the negative edge
        underboundary = (minimum(currpath) + (simdata.boxwidth / 2)) / simdata.boxwidth
        if underboundary > 0
            underboundary = 0
        end
        # create tuple of how many other lines need to be plotted
        linereach = tuple(ceil(abs(underboundary)), ceil(overboundary))

        # add entry if the path crossed either boundary
        if linereach[1] != 0 || linereach[2] != 0 
            addedlines[i] = linereach
        end

    end

    return addedlines
end


mutable struct EnsembleSimData
    xmat::Matrix{Float64}
    smat::Matrix{Float64}
    totaltime::Real
    dt::Real
    boxwidth::Real
end

function timesteps(enssim::EnsembleSimData)
    return 0:enssim.dt:enssim.totaltime
end



"""
Methods to calculate heuristics from data
"""


