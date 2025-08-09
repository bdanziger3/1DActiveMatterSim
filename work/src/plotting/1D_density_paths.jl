########################
# 1D_density_paths.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing code to plot simulation data as 1D paths
########################

using PyPlot
using PyCall
using ColorSchemes
using FFTW

include("../file_management/simdata_files.jl")
include("../utils/paths.jl")
include("../file_management/analysis_files.jl")
include("../analysis/density_distribution_functions.jl")


PARTICLE_PATH_DIRNAME = "1D Partilce Paths"
FILE_NAME_PREFIX = "1d_particle_paths"

D_PLOT_XLABEL = "Position"
D_PLOT_YLABEL = "Time"
D_PLOT_TITLE_TOPLINE = "Particle Density Trajectories"

D_LINE_COLORMAP = ColorSchemes.Blues_8

PARTICLE_PATH_COLOR = "mediumblue"



"""
Produces a 1D path of particle densities in a simulation

Set `savetxt=true` to save a .txt file with the mean spin data results
Set `show` to display the plot.
"""
function plot_1d_density_path(filename::String, saveplot::Bool=false, show::Bool=true)
    
    # get the spin data from the file
    positions, simparams = loadsimpositions(filename)

    # calculate needed nbins
    # nbins = simparams.numparticles
    nbins = 10000
    binwidth = simparams.boxwidth / nbins

    # run density calculation on data
    densitymat::Matrix{Float64} = densitydist(positions, nbins, simparams)


    savetimes = getsavetimes(simparams)
    bincenters = (-simparams.boxwidth / 2) + (binwidth / 2) .+ (binwidth .* (collect(0:nbins-1)))


    #################
    ### Plotting ####
    #################
    # prepare strings for labels used in plot
    local interactionstr::String = ""
    if simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(simparams.interaction) I=$(Int64(round(simparams.interactionfliprate)))  "
    end


    # more file naming controls

    # clear plot
    plt.clf()

    maxdensity = maximum(densitymat)
    
    for t in 1:size(savetimes)[1]
        densitycolors = get(D_LINE_COLORMAP, densitymat[t,:] / maxdensity)
        # plt.plot(bincenters, repeat([savetimes[t]], nbins), marker="o", s=densitymat[t,:], c=densitycolors)
        alphaarray = densitymat[t,:] / maxdensity
        # plt.scatter(bincenters, repeat([savetimes[t]], nbins), c=densitymat[t,:] .^ 2, alpha=1.0, cmap="Blues", marker="s", s=.1)
        plt.scatter(bincenters, repeat([savetimes[t]], nbins), c="blue", alpha=alphaarray, marker="s", s=1)
        # plt.scatter(bincenters, repeat([savetimes[t]], nbins), c=densitymat[t,:], cmap="Blues", marker="o", s=30)
    end

    # plt.plot(densitymat[end,:] / maxdensity)

    gca().invert_yaxis()

    plt.xlabel(D_PLOT_XLABEL)
    plt.ylabel(D_PLOT_YLABEL)
    plt.title("$(D_PLOT_TITLE_TOPLINE)\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")
    

    # if saveplot
    #     # get dir path to save plot
    #     datadirname = getanalysisdir(MEAN_SPIN_DIRNAME, simparams)
    #     outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf")
    #     outputfilename = checkfilename(outputfilename)
    #     plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    # end
    
    if show
        plt.show()
    end
end

"""
Plots the lines of particle positoins.

"""
function plot_1d_path_lines(filename::String, saveplot::Bool=false, show::Bool=true, minindex::Int64=1, maxindex::Int64=0, highlight::Array{Int64}=[0])    
    # get the spin data from the file
    positions, simparams = loadsimpositions(filename)

    savetimes = getsavetimes(simparams)

    # by default max index is the last time index
    if maxindex == 0
        maxindex = size(positions)[1]
    end 


    #################
    ### Plotting ####
    #################
    # prepare strings for labels used in plot
    local interactionstr::String = ""
    if simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(simparams.interaction) I=$(Int64(round(simparams.interactionfliprate)))  "
    end

    # more file naming controls
    
    unwrappedpos = unwrappositions(positions, simparams)

    pathalpha =  .01

    # clear plot
    plt.clf()
    plt.xlim([-simparams.boxwidth / 2, simparams.boxwidth / 2])
    
    @showprogress 1 "plotting paths..." for p in 1:simparams.numparticles
        maxscrens = round(maximum(unwrappedpos[:,p]) / simparams.boxwidth)
        minscrens = round(minimum(unwrappedpos[:,p]) / simparams.boxwidth)
        for s in collect(minscrens:maxscrens)
            plt.plot(unwrappedpos[minindex:maxindex ,p] .- (s .* simparams.boxwidth), savetimes[minindex:maxindex], c=PARTICLE_PATH_COLOR, alpha=pathalpha, linewidth=.1)
        end
    end

    # now highlight certain paths
    highlightcolors = ["orangered", "mediumvioletred", "darkgreen"]
    numhighlightpaths = length(highlight)
    for p_i in 1:numhighlightpaths
        if highlight[p_i] != 0
            maxscrens = round(maximum(unwrappedpos[:,p_i]) / simparams.boxwidth)
            minscrens = round(minimum(unwrappedpos[:,p_i]) / simparams.boxwidth)
            for s in collect(minscrens:maxscrens)
                println("$(minimum(unwrappedpos[minindex:maxindex ,p_i] .+ (s * simparams.boxwidth))), $(maximum(unwrappedpos[minindex:maxindex ,p_i] .+ (s * simparams.boxwidth)))")
                plt.plot(unwrappedpos[minindex:maxindex ,p_i] .- (s * simparams.boxwidth), savetimes[minindex:maxindex], c=highlightcolors[Int64(mod(p_i, numhighlightpaths)) + 1], alpha=1.0, linewidth=1)
            end
        end
    end
        
    gca().invert_yaxis()
    
    plt.xlabel("Position")
    plt.ylabel("Time")
    plt.title("$(D_PLOT_TITLE_TOPLINE)\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")


    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(PARTICLE_PATH_DIRNAME, simparams)
        outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX)_thin.pdf")
        outputfilename = checkfilename(outputfilename)
        println("Saving figure...")
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end

# file = fixpath("work/data/sweeps/alignsimple/boxwidthsweep/d1000/N10000-B10.0-alignsimple-100.0_t10.txt")
# plot_1d_path_lines(file, true, false)
file = fixpath("work/data/sweeps/alignsimple/interactionsweep/N1000-sweep-t1000/extended/N1000-B100.0-alignsimple-500_t1000.txt")
plot_1d_path_lines(file, true, false)
# file = fixpath("work/data/sweeps/alignsimple/interactionsweep/N1000-sweep-t1000/extended/N1000-B100.0-alignsimple-2_t1000.txt")
# plot_1d_path_lines(file, true, false, 98000)
# file = fixpath("work/data/sweeps/turnaway/boxwidthsweep/small_int/d100/N200-B2.0-turnaway-0.1.txt")
# plot_1d_path_lines(file, true, false)

