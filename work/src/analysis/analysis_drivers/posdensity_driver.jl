########################
# posdensity_driver.jl
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run position density function on variou data files
########################

using PyPlot
using PyCall

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../../file_management/analysis_files.jl")



"""
Produces a histogram of all particle densities of the bins of particle positions.

Provide `outputfile` to save the particle density data to that location.

Set `savetxt` to save a .txt file with the densities of each bin for all times in the simulation
Set `logplot` to make the y-axis of the plot logarithmic (Log10)
Set `show` to display the plot as well as saving it.


Set `weighting` to `"particle"` to weigh density counts by the number of particles in the bin.
                Or `"bins"` to weigh density counts by each position bin having an equal weight.
"""
function posdensity_hist(filename::String, posbinwidth::Float64, serialized::Bool=false, saveplot::Bool=true, savetxt::Bool=true, logplot::Bool=false, show::Bool=false, weighting::String="particle")
    
    PARTICLE_DENSITY_DIRNAME = "Particle Density Histogram"
    particledensity_nbins = 30

    FILE_NAME_PREFIX = "particle_density_hist"

    
 
    # get the `SimulationData` from the file
    local sd::SimulationData
    if serialized
        sd = loadcompressedfile(filename)
    else
        sd = loadsim(filename, rowwisetxt)
    end


    # posdensities::Array{Float64} = posdensity_fixedbins(sd)
    posdensities::Array{Float64} = particledensity_alltimes(sd, posbinwidth)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(PARTICLE_DENSITY_DIRNAME, sd.simparams), "$(FILE_NAME_PREFIX)_binwidth-$(posbinwidth).txt")

        open(outputtextfilefullpath, "w") do io
            println(io, "Position Density  Data")
            writedlm(io, posdensities, ",")
        end
    end



    # posdensities = posdensities[end-1000:end, :]

    # use parameters to determine labels
    local logstr::String = ""
    local logstr1::String = raw"$P_{bin}(n)$"
    local logstr2::String = ""
    if logplot
        logstr = "Log of "
        logstr1 = raw"$\log{(P_{bin}(n))}$"
        logstr2 = "-logplot"
    end

    local interactionstr::String = ""
    if sd.simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(sd.simparams.interaction) I=$(Int64(round(sd.simparams.interactionfliprate)))  "
    end

    # determine weighting of densities
    local FILE_NAME_PREFIX
    local binweights
    local weightingdescription
    if lowercase(weighting) == "bins"
        binweights = ones(length(vec(posdensities)))
        weightingdescription = "(Weighted by number of bins with each density)"
        FILE_NAME_PREFIX = "particle_density_hist_binweighted"
    else
        # weighting == particles
        binweights = vec(posdensities)
        weightingdescription = "(Weighted by number of particles that are in each bin)"
        FILE_NAME_PREFIX = "particle_density_hist_particleweighted"
    end



    # more file naming controls
    FILE_NAME_PREFIX_DBINS = "$(FILE_NAME_PREFIX)_dbins-$(particledensity_nbins)"


    ### Plotting

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)
    plt.hist(vec(posdensities), bins=particledensity_nbins, weights=binweights, log=logplot, edgecolor="black", linewidth=.5, zorder=3, density=true)
    plt.xlabel("$(raw"Particle Density $n$")\n(particles per unit length)")
    plt.ylabel("$(logstr)$(raw"Probability Density of Bin with Particle Density $n$")\n$(logstr1)")
    # plt.title(raw"Occurence of Particle Densities")
    # plt.ylabel(raw"Log of Probability Density (Log$P(n)$) of Bin with Particle Density $n$")
    plt.title("Occurence of Particle Densities $(weightingdescription)\n N=$(sd.simparams.numparticles)  Boxwidth=$(sd.simparams.boxwidth)  t=$(Int64(sd.simparams.totaltime))  $(interactionstr)\nBin Width=$(posbinwidth)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(PARTICLE_DENSITY_DIRNAME, sd.simparams)
        plt.savefig(joinpath(datadirname, "$(FILE_NAME_PREFIX_DBINS)_binwidth-$(posbinwidth)$(logstr2).pdf"), bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end


end

datafile = fixpath("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01.txt")
datafile2 = fixpath("work/data/26-6/N5000-B100-alignsimple-100-t4-sn0.01.txt")
datafile_nointeraction = fixpath("work/data/22-6/N10000-nointeraction-t100-sn0.01.txt")
datafile_nointeraction_2 = fixpath("work/data/8-7/N1000-B100.0-nointeraction-100-T100.txt")



serialized::Bool = false
weighting::String = "particles"
# getposdensitydata(datafile, false, datafile_out)


for logplot in [false]
    for binwidth = [1.0, 2.0, 5.0, 8.5, 10.0, 20]
        # posdensity_hist_binweighted(datafile, binwidth, false, true, true, logplot, false)
        posdensity_hist(datafile, binwidth, serialized, true, true, logplot, false, weighting)
    end
end