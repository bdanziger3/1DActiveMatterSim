"""
histogram_helper.jl
Blake Danziger
1D Active Solids
MSc Theoretical Physics Dissertation (2025)

File containing helper functions for building histograms and distribution data.
"""

struct BinSettings
    boxwidth::Float64
    nbins::Int64
    binwidth::Float64
    binboundaries::Array{Float64}


    function # Inner constructor with just nbins input
        BinSettings(boxwidth::Float64, nbins::Int64)
        calculatedbinwidth::Float64 = boxwidth / nbins
        calculatedbinboundaries = collect(range(-boxwidth/2, boxwidth/2, length=nbins+1))
        return new(boxwidth, nbins, calculatedbinwidth, calculatedbinboundaries)
    end
    function # Inner constructor with just binwidth input
        BinSettings(boxwidth::Float64, binwidth::Float64)
        # rounds binwidth to nearest value that evenly divides the `boxwidth`
        calculatednbins::Int64 = Int64(round(boxwidth / binwidth))

        return BinSettings(boxwidth, calculatednbins)
    end
    function # Inner constructor with just binboundaries
        BinSettings(boxwidth::Float64, binboundaries::Array{Float64})
        nbins::Int64 = length(binboundaries) - 1
        calculatedbinboundaries = collect(range(-boxwidth/2, boxwidth/2, length=nbins+1))
        if calculatedbinboundaries != binboundaries
            println("Warning: inputted `binboundaries` parameters are not evenly spaced. New binboundaries calculated.")
        end    

        return new(boxwidth, nbins, binwidth, calculatedbinboundaries)
    end
end