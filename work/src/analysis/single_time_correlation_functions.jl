"""
single_time_correlation_functions.jl
Blake Danziger
1D Active Solids
MSc Theoretical Physics Dissertation (2025)

File containing correlation functions and helper functions
for computing quantities and distributions of interest that are
limited to a single frame in time.

Functions usually take in only a single particle state array and compute
the values of interest from that line only. For values/distributions that change over time.
"""

include("./histogram_helper.jl")


function posdensity(positions::Array{Float64}, spins::Array{Int8}, binsettings::BinSettings)
    

end