########################
# meanclusterlength_driver.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to mean cluster length calculations on simulation data

using PyPlot
using PyCall
using ColorSchemes
using FFTW

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../../file_management/analysis_files.jl")
include("../analysis_functions.jl")


MS_LINE_DEFAULT_COLOR = "darkviolet"
MEAN_SPIN_DIRNAME = "Mean Spin"
FILE_NAME_PREFIX = "mean_spin"

MS_PLOT_XLABEL = "Time"
MS_PLOT_YLABEL = "Mean Spin \n$(raw"$\langle S_i \rangle_{i} = \frac{1}{N}\sum^N_i S_i$")"
MS_PLOT_TITLE_TOPLINE = "Mean Spin of Active Particle Simulation"

MS_LINE_DENSITY_COLORMAP = ColorSchemes.Reds_8
MS_LINE_INTERACTION_COLORMAP = ColorSchemes.Blues_8
MS_LINE_BOXWIDTH_COLORMAP = ColorSchemes.Greens_8


PO_LINE_DEFAULT_COLOR = "darkviolet"
POLAR_ORDER_DIRNAME = "Polar Order"
PO_FILE_NAME_PREFIX = "polar_order"

PO_PLOT_YLABEL = "Polar Order \n$(raw"$\frac{1}{N}\sum^N_i \left| S_i \right| $")"
PO_PLOT_TITLE_TOPLINE = "Polar Order of Active Particle Simulation"







"""
Produces a line plot of the mean spin as a function of the time

Set `savetxt=true` to save a .txt file with the mean spin data results
Set `show` to display the plot.
"""
function mcl_plot(filename::String, saveplot::Bool=false, savetxt::Bool=false, show::Bool=true)
    
    # get the spin data from the file
    nsegments, serialized = getsimfiletype(filename)


    positions, simparams = loadsimpositions(filename)
    spins, simparams = loadsimspins(filename)
   
    local simdata
    if serialized
        simdata = loadcompressedfile(filename)
    else
        simdata = loadsim(filename, rowwisetxt)
    end

    mclmat = meanclusterlength(positions, spins, simparams)
    # meanclusterlength(simdata)

    # # save as a datafile if requested
    # if savetxt
    #     outputtextfilefullpath = joinpath(getanalysisdir(MEAN_SPIN_DIRNAME, simparams), "$(FILE_NAME_PREFIX).txt")
    #     outputtextfilefullpath = checkfilename(outputtextfilefullpath)
    #     open(outputtextfilefullpath, "w") do io
    #         println(io, "Mean Spin Data")
    #         println(io, "Simulation Parameters")
    #         println(io, csv_serialize(simparams))
    #         println(io, "Plot Data")
    #         writedlm(io, pomatrix, ",")
    #     end
    # end

    #################
    ### Plotting ####
    #################
    # prepare strings for labels used in plot
    # local interactionstr::String = ""
    # if simparams.interaction == nointeraction
    #     interactionstr = "no-interaction  "
    # else
    #     interactionstr = "$(simparams.interaction) I=$(Int64(round(simparams.interactionfliprate)))  "
    # end


    # more file naming controls

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)
    plt.plot(mclmat[1, :], mclmat[2, :])
    # plt.xlabel(MS_PLOT_XLABEL)
    # plt.ylabel(MS_PLOT_YLABEL)
    # plt.title("$(MS_PLOT_TITLE_TOPLINE)\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")
    

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

p1 = fixpath("work/data/sweeps/alignsimple/interactionsweep/N1000-sweep-t1000/N1000-B100.0-alignsimple-50.txt")
p2 = fixpath("work/data/sweeps/alignsimple/interactionsweep/N1000-sweep-t1000/N1000-B100.0-alignsimple-100.txt")
p3 = fixpath("work/data/sweeps/alignsimple/interactionsweep/N1000-sweep-t1000/N1000-B100.0-alignsimple-300.txt")

mcl_plot(p1)

# """
# Plots the data previously saved to a `.txt` file.

# Can specify `minindex` and `maxindex` to only plot that many data points along the horizontal axis.
# """
# function ms_saveddata_plot(txtfilename::String, minindex::Int64=1, maxindex::Int64=0, saveplot::Bool=false, show::Bool=true)
    

#     # get the plot data from the file
#     xydatamat, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

#     # handle default `maxindex` behavior
#     if maxindex == 0
#         maxindex = size(xydatamat)[2]
#     end

#     ### Plotting
#     # prepare plot
#     plt.clf()
#     plt.grid(true, zorder=0)
#     ###

#     local paramsstr
#     if sweeptype != nosweep
#         filenameprefix = "mean_spin_plot_$(sweeptype)"
#         local legendlabels::Array{String}
#         local legendtitle::String
#         local sweepcm
#         if sweeptype == densitysweep
#             legendtitle = "Number of Particles"
#             sweepcm = MS_LINE_DENSITY_COLORMAP
#         elseif sweeptype == boxwidthsweep
#             legendtitle = "Box Width"
#             sweepcm = MS_LINE_BOXWIDTH_COLORMAP
#         elseif sweeptype == interactionstrengthsweep
#             legendtitle = "Interaction Fliprate"
#             sweepcm = MS_LINE_INTERACTION_COLORMAP
#         end
        
#         legendlabels = (val -> "$val").(legendvalues)
        
#         # get params string for plot
#         paramsstr = simparamsout
#         nshots = length(legendvalues)

#         for i in 1:nshots            
#             # plot each line one at a time
#             linecolor = get(sweepcm, (i / nshots))
#             plt.plot(xydatamat[1, minindex:maxindex], xydatamat[i+1, minindex:maxindex], color=(linecolor.r, linecolor.g, linecolor.b))
#         end
        
#         # add legends and labels to plot
#         plt.legend(legendlabels, title=legendtitle)

#     else
#         filenameprefix = "mean_spin_plot"

#         # prepare strings for labels used in plot
#         local interactionstr::String = ""
#         if simparamsout.interaction == nointeraction
#             interactionstr = "no-interaction  "
#         else
#             interactionstr = "$(simparamsout.interaction) I=$(Int64(round(simparamsout.interactionfliprate)))  "
#         end
        
#         # get params string for plot
#         paramsstr = "N=$(simparamsout.numparticles)  Boxwidth=$(simparamsout.boxwidth)  t=$(Int64(simparamsout.totaltime))  $(interactionstr)"

#         # plot the line
#         plt.plot(xydatamat[1, minindex:maxindex], xydatamat[2, minindex:maxindex], color=MS_LINE_DEFAULT_COLOR)
#     end

#     plt.xlabel(MS_PLOT_XLABEL)
#     plt.ylabel(MS_PLOT_YLABEL)
#     plt.title("$(MS_PLOT_TITLE_TOPLINE)\n$(paramsstr)")        
    
#     if saveplot
#         # get dir path to save plot
#         datadirname = dirname(txtfilename)

#         ndatapoints = size(xydatamat)[2]
#         local rangestr::String = ""
#         if minindex != 0 && maxindex != ndatapoints
#             rangestr = "_indexrange_$minindex-$maxindex"
#         elseif minindex == 1 && maxindex != ndatapoints
#             rangestr = "_maxindex-$maxindex"
#         elseif minindex != 1 && maxindex == ndatapoints
#             rangestr = "_minindex-$minindex"
#         end

#         outputfilename = joinpath(datadirname, "$(filenameprefix)_from_txt$rangestr.pdf")
#         outputfilename = checkfilename(outputfilename)
#         plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
#     end
    
#     if show
#         plt.show()
#     end
# end


# """
# Loads the data previously saved to a `.txt` file.

# Returns a tuple of 4 objects:


#     `xydatamatrix`:   The matrix holding the analysis data. Row 0 holds the x-axis data.
#                         The remaining lines hold the y-axis data
                        
#     `sweeptype`:        The `SweepType` of the data. Is `nosweep` if not a sweep.
    
#     `simparamsout`:     The simulation parameters used. Returns a `SimulationParameters` object if not a sweep.
#                         Returns a string to use in the plot title if it is a sweep

#     `legendvalues`:     The list of values of the swept parameter in the sweep. If the file is not for a sweep, contains `[nothing]`.
# """
# function ms_saveddata_load(txtfilename::String)
    
#     # get the plot data from the file
#     file = open(txtfilename, "r")

#     line = readline(file)
#     # Check that data file header is correct
#     @assert contains(lowercase(line), "mean spin data") "MS data file does not have correct first line. File may be corrupted."
    
#     # check if it is a sweep
#     line = readline(file)

#     local paramsstr
#     local podata::Array{Float64} = zeros(0, 0)
#     local tdata
#     local simparamsout
#     local legendvalues

#     if contains(lowercase(line), "sweep type")
#         sweeptype = eval(Symbol(line[length("sweep type: ")+1:end]))

#         local legendlabels::Array{String}
#         local legendtitle::String
#         # local paramsstr::String

#         # now get sweep values
#         line = readline(file)
#         line = readline(file)
        
#         legendlabels = split(line, ",")
#         legendvalues = parse.(Float64, legendlabels)
        
#         # get params string for plot
#         line = readline(file)
#         paramsstr = line[length("Params String: ")+1:end]
#         simparamsout = paramsstr
       

#     else  
#         sweeptype = nosweep      
#         # get params
#         # read next line to get number of position data points
#         simparaminfo_str = readline(file)
#         simparams = SimulationParameters(simparaminfo_str)
#         simparamsout = simparams

#         legendvalues = [nothing]
#     end

#     # skip to data matrix
#     while !contains(lowercase(line), "data matrix") && !contains(lowercase(line), "plot data")
#         line = readline(file)
#     end

     
#     # read x data (t)
#     dataline = readline(file)
#     tdata = parse.(Float64, split(dataline, ","))
    
#     # read y data (MS)
#     nshots = length(legendvalues)
#     podata = zeros(nshots, length(tdata))

#     # load each line
#     for i in 1:nshots
#         dataline = readline(file)
#         podata[i, :] = parse.(Float64, split(dataline, ","))
#     end

#     xydatamatrix = zeros(size(podata)[1] + 1, size(podata)[2])
#     xydatamatrix[1, :] = tdata
#     xydatamatrix[2:end, :] = podata

#     return xydatamatrix, sweeptype, simparamsout, legendvalues
# end


# """
# Produces a line plot of the mean spin as a function of the time interval `t`
# of all the files in a given directory, and plots the lines on the same axes.

# Set `savetxt` to save a .txt file with the orientation self-correlation function results
# Set `show` to display the plot as well as saving it.
# """
# function ms_plot_dir(dirname::String, sweepname::String, sweeptype::SweepType, saveplot::Bool=true, savetxt::Bool=true, individualplots::Bool=false, show::Bool=false)
    
#     filename_prefix = "mean_spin_sweep_plot"
    
#     # initialize data matrix
#     local full_pomatrix::Matrix{Float64} = zeros(0,0)
    
#     # load every file in the directory
#     datafiles_all = readdir(dirname)
#     datafiles = []

#     # only look at .txt files
#     for file in datafiles_all
#         if endswith(file, ".txt")
#             push!(datafiles, file)
#         end
#     end

#     sp_array::Array{SimulationParameters} = []
#     for (i, filename) in enumerate(datafiles)
#         inputfilename = joinpath(dirname, filename)
        
#         # run correlation function and add to plot
#         # get the spins from the file
        
#         # load the spins from the file
#         spins, simparams = loadsimspins(inputfilename)
#         # add simparams to array
#         push!(sp_array, simparams)

#         # run mean spin calculation on data
#         pomatrix_i::Matrix{Float64} = meanspin(spins, simparams)

#         # if `full_pomatrix` is empty, add x-axis data and first line of y-axis data
#         if length(full_pomatrix) == 0
#             full_pomatrix = zeros(length(datafiles) + 1, size(pomatrix_i)[2])

#             full_pomatrix[1:2, :] = pomatrix_i
#         else
#             # else just add y-axis data
#             full_pomatrix[i+1, :] = pomatrix_i[2,:]
#         end
#     end

    
#     # generate legend labels and title

#     local interactionstr::String = ""
#     if sp_array[1].interaction == nointeraction
#         interactionstr = "no-interaction  "
#     else
#         interactionstr = "$(sp_array[1].interaction) I=$(sp_array[1].interactionfliprate)"
#     end
    
#     local sweepvalues::Array{Real}
#     local legendtitle::String
#     local paramsstr::String
#     local sortedorder::Vector{Int}
#     local sweepcm
#     if sweeptype == densitysweep
#         legendtitle = "Number of Particles"
#         sortedorder = sortperm([sp.numparticles for sp in sp_array])
#         sweepvalues = (sp -> sp.numparticles).(sp_array[sortedorder])
#         paramsstr = "t=$(Int64(sp_array[1].totaltime))  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
#         sweepcm = MS_LINE_DENSITY_COLORMAP
#     elseif sweeptype == boxwidthsweep
#         legendtitle = "Box Width"
#         sortedorder = sortperm([sp.boxwidth for sp in sp_array])
#         sweepvalues = (sp -> sp.boxwidth).(sp_array[sortedorder])
#         paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  $(interactionstr)"
#         sweepcm = MS_LINE_BOXWIDTH_COLORMAP
#     elseif sweeptype == interactionstrengthsweep
#         interactionstr = "$(sp_array[1].interaction)"
#         legendtitle = "Interaction Fliprate"
#         sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
#         sweepvalues = (sp -> sp.interactionfliprate).(sp_array[sortedorder])
#         paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
#         sweepcm = MS_LINE_INTERACTION_COLORMAP
#     end
#     legendlabels::Array{String} = (val -> "$(val)").(sweepvalues)


#     # save as a datafile if requested
#     if savetxt
#         outputtextfilefullpath = joinpath(getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array[sortedorder], sweepname), "$(filename_prefix).txt")
#         outputtextfilefullpath = checkfilename(outputtextfilefullpath)

#         open(outputtextfilefullpath, "w") do io
#             println(io, "Mean Spin Data")
#             println(io, "Sweep Type: $(sweeptype)")
#             println(io, "Sweep Values")
#             writedlm(io, permutedims(legendlabels), ",")
#             println(io, "Params String: $(paramsstr)")
#             println(io, "Data Matrix")
#             writedlm(io, permutedims(full_pomatrix[1, :]), ",")
#             for i in sortedorder
#                 writedlm(io, permutedims(full_pomatrix[i+1, :]), ",")
#             end
#         end
#     end





#     #################
#     ### Plotting ####
#     #################

#     # if `individualplots`, plot the data seperately
#     if individualplots
#         for (i, s) in enumerate(sortedorder)        
#             # clear plot
#             plt.clf()

#             plt.grid(true, zorder=0)
#             plt.plot(full_pomatrix[1, :], full_pomatrix[s+1, :])
#             plt.xlabel(MS_PLOT_XLABEL)
#             plt.ylabel(MS_PLOT_YLABEL)
#             plt.title("$(MS_PLOT_TITLE_TOPLINE)\n N=$(sp_array[s].numparticles)  Boxwidth=$(sp_array[s].boxwidth)  t=$(Int64(sp_array[s].totaltime))  $(sp_array[1].interaction) I=$(sp_array[s].interactionfliprate)")

#             # figure out swept value and var and add to names
#             datadirname = getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array, sweepname)
#             outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX)_$(sweepvalues[i]).pdf")
#             outputfilename = checkfilename(outputfilename)
#             plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
#         end
#     end

#     # clear plot
#     plt.clf()
    
#     plt.grid(true, zorder=0)

#     for i in 1:length(sortedorder)
#         linecolor = get(sweepcm, sweepvalues[i] / sweepvalues[end])
#         plt.plot(full_pomatrix[1,:], full_pomatrix[i+1,:], color=(linecolor.r, linecolor.g, linecolor.b))
#     end

#     plt.legend(legendlabels, title=legendtitle)
#     plt.xlabel(MS_PLOT_XLABEL)
#     plt.ylabel(MS_PLOT_YLABEL)
#     plt.title("$(MS_PLOT_TITLE_TOPLINE)\n$(paramsstr)")
    
#     if saveplot
#         # get dir path to save plot
#         datadirname = getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array, sweepname)
#         outputfilename = joinpath(datadirname, "$(filename_prefix).pdf")
#         outputfilename = checkfilename(outputfilename)
#         plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
#     end
    
#     if show
#         plt.show()
#     end
# end


# """
# Plots the data previously saved to a `.txt` file.

# Can specify `minindex` and `maxindex` to only plot that many data points along the horizontal axis.
# """
# function po_saveddata_plot(txtfilename::String, minindex::Int64=1, maxindex::Int64=0, saveplot::Bool=false, show::Bool=true)
    

#     # get the plot data from the file
#     xydatamat, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

#     # handle default `maxindex` behavior
#     if maxindex == 0
#         maxindex = size(xydatamat)[2]
#     end

#     ### Plotting
#     # prepare plot
#     plt.clf()
#     plt.grid(true, zorder=0)
#     ###

#     local paramsstr
#     local xaxistitle
#     if sweeptype != nosweep
#         filenameprefix = "polar_order_plot_$(sweeptype)"
#         if sweeptype == densitysweep
#             xaxistitle = "Number of Particles"
#             sweepcm = MS_LINE_DENSITY_COLORMAP
#         elseif sweeptype == boxwidthsweep
#             xaxistitle = "Box Width"
#             sweepcm = MS_LINE_BOXWIDTH_COLORMAP
#         elseif sweeptype == interactionstrengthsweep
#             xaxistitle = "Interaction Fliprate"
#             sweepcm = MS_LINE_INTERACTION_COLORMAP
#         end
                
#         # get params string for plot
#         paramsstr = simparamsout
#         nshots = length(legendvalues)

#         avgpo::Vector{Float64} = zeros(length(legendvalues))
#         for i in 1:nshots
#             # calculate average PO for the given range
#             avgpo[i] = sum(abs.(xydatamat[i+1, minindex:maxindex]))[1] / (maxindex - minindex + 1)
#         end

#         # plot results
#         linecolor = get(sweepcm, 1.0)
#         plt.plot(legendvalues, avgpo, "o")
#         axes = plt.gca()
#         axes.loglog()

#     else
#         filenameprefix = "polar_order_plot"

#         # prepare strings for labels used in plot
#         local interactionstr::String = ""
#         if simparamsout.interaction == nointeraction
#             interactionstr = "no-interaction  "
#         else
#             interactionstr = "$(simparamsout.interaction) I=$(Int64(round(simparamsout.interactionfliprate)))  "
#         end

#         absxymat = abs.(xydatamat[2, minindex:maxindex])
#         po_sum = sum(absxymat)
#         avgpo = po_sum / (maxindex - minindex + 1)
#         println(avgpo)

        
        
#     end


#     plt.xlabel(MS_PLOT_XLABEL)
#     plt.ylabel(PO_PLOT_YLABEL)
#     plt.title("$(PO_PLOT_TITLE_TOPLINE)\n$(paramsstr)")        
    
#     # if saveplot
#     #     # get dir path to save plot
#     #     datadirname = dirname(txtfilename)

#     #     ndatapoints = size(xydatamat)[2]
#     #     local rangestr::String = ""
#     #     if minindex != 0 && maxindex != ndatapoints
#     #         rangestr = "_indexrange_$minindex-$maxindex"
#     #     elseif minindex == 1 && maxindex != ndatapoints
#     #         rangestr = "_maxindex-$maxindex"
#     #     elseif minindex != 1 && maxindex == ndatapoints
#     #         rangestr = "_minindex-$minindex"
#     #     end

#     #     outputfilename = joinpath(datadirname, "$(filenameprefix)_from_txt$rangestr.pdf")
#     #     outputfilename = checkfilename(outputfilename)
#     #     plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
#     # end
    
#     if show
#         plt.show()
#     end
# end

# file = fixpath("work/analysis/Mean Spin/Align Simple/small_interaction2_3-8/mean_spin_sweep_plot.txt")
# # data = ms_saveddata_load(file)

# po_saveddata_plot(file)
