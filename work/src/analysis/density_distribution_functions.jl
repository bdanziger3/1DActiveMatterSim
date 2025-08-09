using CSV
using Statistics
using StatsBase
using ProgressMeter

include("../simulation/basic_sim.jl")
include("./histogram_helper.jl")



"""
Density Distribution

Returns bin counts of particles in simulation with `nbins` bins
"""
function densitydist(positions::Array, nbins::Int64, simparams::SimulationParameters)::Matrix{<:Real}
    # get nsaves
    nsaves = getnsaves(simparams)
    
    binwidth::Float64 = simparams.boxwidth / nbins
    
    bincounts = zeros(nsaves, nbins)
    binpositions = Int64.(div.(positions .+ (simparams.boxwidth / 2), binwidth)) .+ 1
    println(binpositions[1,1:10])
    println(minimum(binpositions))
    println(maximum(binpositions))

    # populate density matrix
    @showprogress 1 "calculating densities..." for t_i in 1:nsaves
        for binpos in binpositions[t_i, :]
            bincounts[t_i, binpos] += 1
        end
    end

    return bincounts ./ simparams.numparticles
end 
