
"""
single_time_correlation_functions.jl
Blake Danziger
1D Active Matter Sim
MSc Theoretical Physics Dissertation (2025)

File containing correlation functions and helper functions
for computing quantities and distributions of interest that are
limited to a single frame in time.

Functions usually take in only a single particle state array and compute
the values of interest from that line only. For values/distributions that change over time.
"""

include("./histogram_helper.jl")
include("../simulation/sim_structs.jl")

# binwidth with respect to the minimum space step `dx`
dx_per_bin::Float64 = 10
n_density_bins::Int64 = 100


"""
Takes a list of positions and bins them based on the `binsettings`,
then counts the number of particles in each bin and divides by the width of the bin.
Returns an array of particle densities for each bin with the first bin starting at -boxwidth/2 and the final bin ending at boxwidth/2.
"""
function posdensityhist(positions::Array{Float64}, simparams::SimulationParameters)::Array{Float64}

    posbinsettings::BinSettings = positionbinsettings(simparams.boxwidth, dx_per_bin * dx(simparams))

    # bin the particle positions
    posbincounts::Array{Int64} = bincounts(positions, posbinsettings)

    # convert to particle density (particles / unit of length)
    particledensities = posbincounts ./ posbinsettings.binwidth

    # now bin the particle densities
    densitiesbinsettings = BinSettings(maximum(particledensities), n_density_bins)
    posdensitycounts::Array{Int64} = bincounts(particledensities, densitiesbinsettings)

    return posdensitycounts    
end