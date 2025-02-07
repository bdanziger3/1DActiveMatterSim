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


mutable struct SimData
    particlepaths::Array{ParticlePath}
    totaltime::Real
    dt::Real
end

function npaths(simdata::SimData)
    return length(simdata.particlepaths)
end





"""
Methods to calculate heuristics from data
"""


