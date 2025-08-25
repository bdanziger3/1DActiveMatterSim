########################
# density_plots.jl
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# Functions to plot a density histogram of simulation data
########################


using CSV
using PyPlot
using PyCall
using ColorSchemes

include("../simulation/sim_structs.jl")
include("../file_management/simdata_files.jl")


function getpatchpolarity(positions, spins, binleftedge, binrightedge)
    # find average spin spins of all particles within given bounds

    # get indices of all particles in the bin
    inbin = findall(x -> x >= binleftedge && x <= binrightedge, positions)

    return sum((spins[inbin]))  / length(inbin)
end


function plotparticledensityline(simdata::SimulationData, timeindex::Int64, particledensity_nbins::Int64=-1)

    positions = simdata.positions[timeindex, :]
    spins = simdata.spins[timeindex, :]

    if particledensity_nbins == -1
        particledensity_nbins = Int64(round(simdata.simparams.boxwidth))

        # except if it's 1
        if particledensity_nbins == 1
            particledensity_nbins = 100
        end
    end

    local interactionstr::String = ""
    if simdata.simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(simdata.simparams.interaction) I=$(Int64(round(simdata.simparams.interactionfliprate)))  "
    end

    barcolors = repeat(["red", "darkviolet"], 5)

    ### Plotting

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)
    n, bins, patches = plt.hist(positions, bins=particledensity_nbins, edgecolor="black", linewidth=.5, zorder=3, density=true)
    plt.xlabel("Particle Position")
    plt.ylabel("Particle Density")
    # plt.title(raw"Occurence of Particle Densities")
    # plt.ylabel(raw"Log of Probability Density (Log$P(n)$) of Bin with Particle Density $n$")
    plt.title("Particle Densities\n N=$(simdata.simparams.numparticles)  Boxwidth=$(simdata.simparams.boxwidth)  t=$(Int64(simdata.simparams.totaltime))  $(interactionstr)")
    

    cs = ColorSchemes.coolwarm

    # Calculate polarities of bins
    for i in 1:particledensity_nbins
        binleftedge = bins[i]
        binrightedge = bins[i + 1]

        binpolarity = (getpatchpolarity(positions, spins, binleftedge, binrightedge) / 2) + .5

        bincolor = get(cs, binpolarity)
        println(binpolarity)
        println(bincolor)
        println(fieldnames(typeof(bincolor)))
        patches[i].set_facecolor((bincolor.r, bincolor.g, bincolor.b))

    end

    # Create a ScalarMappable for the colorbar
    norm = mcolors.Normalize(vmin=0, vmax=100)
    cmap = cm.viridis
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # required for older matplotlib

    # Add the colorbar
    cbar = plt.colorbar(sm)
    cbar.set_label("Custom Data")
    cbar.set_ticks([0, 100])  # labels on the ends
    cbar.ax.set_yticklabels(["Low (0)", "High (100)"])

    plt.colorbar()

    println(n)
    println(bins)
    println(patches)
    plt.show()


end




file_align = fixpath("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01.txt")

sd = loadsim(file_align, rowwisetxt)

plotparticledensityline(sd, 100, 20)

