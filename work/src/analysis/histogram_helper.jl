########################
# histogram_helper.jl
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# File containing helper functions for building histograms and distribution data.
########################

using StatsBase


struct BinSettings
    valuesrange::Float64
    nbins::Int64
    binwidth::Float64
    binboundaries::Vector{Float64}

    # Inner constructors assume data in the range [0, valuesrange]
    function # Inner constructor with all fields
        BinSettings(valuesrange::Float64, nbins::Int64, binwidth::Float64, binboundaries::Vector{Float64})
        return new(valuesrange, nbins, binwidth, binboundaries)
    end
    function # Inner constructor with just nbins input for particle positions
        BinSettings(valuesrange::Float64, nbins::Int64)
        calculatedbinwidth::Float64 = valuesrange / nbins
        calculatedbinboundaries::Vector{Float64} = collect(range(0, valuesrange, length=nbins+1))
        return BinSettings(valuesrange, nbins, calculatedbinwidth, calculatedbinboundaries)
    end
    function # Inner constructor with just binwidth input
        BinSettings(valuesrange::Float64, binwidth::Float64)
        # rounds binwidth to nearest value that evenly divides the `valuesrange`
        calculatednbins::Int64 = Int64(round(valuesrange / binwidth))

        return BinSettings(valuesrange, calculatednbins)
    end
    function # Inner constructor with just binboundaries
        BinSettings(valuesrange::Float64, binboundaries::Array{Float64})
        nbins::Int64 = length(binboundaries) - 1
        calculatedbinboundaries::Array{Float64} = collect(range(-valuesrange/2, valuesrange/2, length=nbins+1))
        if calculatedbinboundaries != binboundaries
            println("Warning: inputted `binboundaries` parameters are not evenly spaced. New binboundaries calculated.")
        end    

        return BinSettings(valuesrange, nbins, binwidth, calculatedbinboundaries)
    end


end

# Position BinSettings constructors assume values in the range [-boxwidth/2, boxwidth/2]
function # outer constructor with just nbins input for particle positions
    positionbinsettings(boxwidth::Float64, nbins::Int64)
    calculatedbinwidth::Float64 = boxwidth / nbins
    calculatedbinboundaries = collect(range(-boxwidth/2, boxwidth/2, length=nbins+1))
    return BinSettings(boxwidth, nbins, calculatedbinwidth, calculatedbinboundaries)
end
function # Outer constructor with just binwidth input
    positionbinsettings(boxwidth::Float64, binwidth::Float64)
    # rounds binwidth to nearest value that evenly divides the `boxwidth`
    calculatednbins::Int64 = Int64(round(boxwidth / binwidth))

    return positionbinsettings(boxwidth, calculatednbins)
end
function # Outer constructor with just binboundaries
    positionbinsettings(boxwidth::Float64, binboundaries::Array{Float64})
    nbins::Int64 = length(binboundaries) - 1
    calculatedbinboundaries = collect(range(-boxwidth/2, boxwidth/2, length=nbins+1))
    if calculatedbinboundaries != binboundaries
        println("Warning: inputted `binboundaries` parameters are not evenly spaced. New binboundaries calculated.")
    end    

    return BinSettings(boxwidth, nbins, binwidth, calculatedbinboundaries)
end

# get midpoints of bins
function midpoints(binsettings::BinSettings)::Array{Float64}
    return binsettings.binboundaries[1:end-1] .+ (binsettings.binwidth / 2)
end



"""
Takes a list of values and bins them based on the `binsettings`,
then counts the number of occurences in each bin.
Returns an array of particle counts for each bin.
"""
function bincounts(vals::Array{Float64}, binsettings::BinSettings)::Array{Int64}

    @assert minimum(vals) >= binsettings.binboundaries[1] && maximum(vals) <= binsettings.binboundaries[end] "Values are not all within the range provided by `binsettings`."

    # sort positions into bins and count
    shiftedvalues = vals .- binsettings.binboundaries[1]
    binplacements = div.(shiftedvalues, binsettings.binwidth) .+ 1

    countsdict = countmap(binplacements)

    # get array with number of particles in each bin
    nparticlesinbin = (x -> x in keys(countsdict) ? countsdict[x] : 0).(collect(1:binsettings.nbins))

    return nparticlesinbin
end