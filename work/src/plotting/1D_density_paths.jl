########################
# 1D_density_paths.jl
# Blake Danziger
# 1D Active Matter Sim
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

PARTICLE_PATH_COLOR = "teal"

LR = "lightcoral"
LB = "lightskyblue"


FONT = "Times New Roman"



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
    
    # for t in 1:size(savetimes)[1]
    @showprogress 1 "plotting densities..." for t in 1:size(savetimes)[1]
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
function plot_1d_path_lines(filename::String, saveplot::Bool=false, show::Bool=true, mintime::Real=0, maxtime::Real=-1, highlight::Array{Int64}=[0], title::Bool=true, alpha::Real=0.05, lw::Real=1)    
    # get the spin data from the file
    positions, simparams = loadsimpositions(filename)

    savetimes = getsavetimes(simparams)

    # by default max index is the last time index
    if maxtime < 0
        maxtime = simparams.totaltime
    end 

    # convert times to indices
    minindex = Int64(round(mintime / simparams.snapshot_dt)) + 1
    maxindex = Int64(round(maxtime / simparams.snapshot_dt)) + 1


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

    pathalpha =  alpha

    # clear plot
    plt.clf()
    plt.xlim([-simparams.boxwidth / 2, simparams.boxwidth / 2])
    plt.xlim([-25, -15])
    plt.ylim([mintime, maxtime])
    
    @showprogress 1 "plotting paths..." for p in 1:simparams.numparticles
        maxscrens = round(maximum(unwrappedpos[:,p]) / simparams.boxwidth)
        minscrens = round(minimum(unwrappedpos[:,p]) / simparams.boxwidth)
        for s in collect(minscrens:maxscrens)
            plt.plot(unwrappedpos[minindex:maxindex ,p] .- (s .* simparams.boxwidth), savetimes[minindex:maxindex], c=PARTICLE_PATH_COLOR, alpha=pathalpha, linewidth=lw)
        end
    end

    # now highlight certain paths
    # highlightcolors = [ "orangered", "mediumvioletred", "darkgreen"]
    highlightcolors = ["k"]
    numhighlightpaths = length(highlight)
    for p_i in 1:numhighlightpaths
        if highlight[p_i] != 0
            p = highlight[p_i]
            maxscrens = round(maximum(unwrappedpos[:,p]) / simparams.boxwidth)
            minscrens = round(minimum(unwrappedpos[:,p]) / simparams.boxwidth)
            for s in collect(minscrens:maxscrens)
                plt.plot(unwrappedpos[minindex:maxindex ,p] .- (s * simparams.boxwidth), savetimes[minindex:maxindex], c=highlightcolors[Int64(mod(p_i, length(highlightcolors))) + 1], alpha=1.0, linewidth=1)
            end
        end
    end
    
    ax = gca()
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", which="both", direction="out")
    
    plt.xlabel(raw"Position", fontsize=16, fontname=FONT)
    plt.ylabel(raw"$\longleftarrow$ Time", fontsize=16, fontname=FONT)

    for label in cat(ax.get_xticklabels(), ax.get_yticklabels(), dims=1)
        label.set_fontname(FONT)  # font family
    end

    # xticks = Int64.(round.(collect(-simparams.boxwidth / 2:5:simparams.boxwidth / 2)))
    # xticks = Int64.(round.(collect(minindex / simparams.snapshot_dt:5:maxindex / simparams.snapshot_dt)))
    # plt.xticks([-50, -25, 0, 25, 50], [0, 25, 50, 75, 100])
    # plt.xticks([-10, -5, 0, 5, 10], [0, 5, 10, 15, 20])
    # plt.xticks([-25, -15, -5, 5, 15, 25], [0, 10, 20, 30, 40, 50])
    plt.xticks([-25, -15], [0, 10])
    # plt.xticks([-10, -5, 0, 5, 10], [0, 5, 10, 15, 20])
    # plt.yticks([80, 85, 90, 95, 100])
    plt.yticks(collect(mintime:5:maxtime))



    if title
        plt.title("$(D_PLOT_TITLE_TOPLINE)\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")
    # else
        
        # ax.spines["bottom"].set_visible(false)
        # ax.spines["top"].set_visible(false)
    end


    timestr = ""
    if mintime != 0 || maxtime != simparams.totaltime
        timestr = "_times_$(mintime)_to_$(maxtime)"
    end

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(PARTICLE_PATH_DIRNAME, simparams)
        outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX)$(timestr).pdf")
        outputfilename = checkfilename(outputfilename)
        println("Saving figure...")

        # calculate appropriate plot dimensions
        # fig = plt.gcf()
        # figdims = fig.get_size_inches()
        # println(size)

        # fig.set_size_inches(figdims[1], 3 * figdims[1])
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end


# path = fixpath("work/data/sweeps/alignsimple/interactionsweep/Aug15-B50-Isweep/N500-B50.0-alignsimple-0.5_rep2.txt")
# path = fixpath("work/data/13-8/N2000-B20.0-alignsimple-10-T100.txt")
# path = fixpath("work/data/sweeps/alignsimple/densitysweep/Aug13-density-sweep-B50/N500-B50.0-alignsimple-100.0_rep3.txt")
# path = fixpath("work/data/sweeps/alignsimple/densitysweep/Aug13-density-sweep-B50/N2500-B50.0-alignsimple-100.0_rep5.txt")
# plot_1d_path_lines(path, true, false, 480, -1, [0], false)
# path = fixpath("work/data/19-8/N2000-B20.0-nointeraction-0-T100.txt")

# for i in [2, 5, 10, 20, 50, 100]
#     path = fixpath("work/data/sweeps/turnaway/interactionsweep/N500-B50/interactionsweep/N500-B50.0-turnaway-$i.0_rep5.txt")
#     plot_1d_path_lines(path, true, false, 480, 500, [1, 70, 400], false, .1, 1)
# end
for i in [2]
    for j in 1:5
        path = fixpath("work/data/sweeps/turnaway/interactionsweep/N500-B50/interactionsweep/N500-B50.0-turnaway-$i.0_rep5.txt")
        plot_1d_path_lines(path, true, false, 480, 500, (20 * j) .+ [10,11,12,13,14,15,16,17,18,19,20], false, .8, .4)
    end
end
# plot_1d_path_lines(path, true, false, 480, 500, [1, 5, 70], false, 0.2, .1)
