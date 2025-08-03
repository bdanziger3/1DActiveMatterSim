########################
# polarorder_driver.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run polar order analysis functions on various data files and plot results
########################

using PyPlot
using PyCall
using ColorSchemes
using FFTW

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../../file_management/analysis_files.jl")
include("../analysis_functions.jl")


@enum SweepType densitysweep boxwidthsweep interactionstrengthsweep nosweep


PO_LINE_DEFAULT_COLOR = "darkviolet"
POLAR_ORDER_DIRNAME = "Polar Order"
FILE_NAME_PREFIX = "polar_order"

PO_PLOT_XLABEL = "Time"
PO_PLOT_YLABEL = "Polar Order \n$(raw"$\langle S_i \rangle_{i} = \frac{1}{N}\sum^N_i S_i$")"
PO_PLOT_TITLE_TOPLINE = "Polar Order of Active Particle Simulation"

PO_LINE_DENSITY_COLORMAP = ColorSchemes.Reds_8
PO_LINE_INTERACTION_COLORMAP = ColorSchemes.Blues_8
PO_LINE_BOXWIDTH_COLORMAP = ColorSchemes.Greens_8






"""
Produces a line plot of the polar order as a function of the time

Set `savetxt=true` to save a .txt file with the polar order data results
Set `show` to display the plot.
"""
function po_plot(filename::String, saveplot::Bool=false, savetxt::Bool=false, show::Bool=true)
    
    # get file info
    nsegments, serialized = getsimfiletype(filename)

    @assert nsegments == 1 "More than 1 segment in file $filename. Found $nsegments segments"

    # get the spin data from the file
    spins, simparams = loadsimspins(inputfilename)

    # run polar order calculation on data
    pomatrix::Matrix{Float64} = polarorder(spins, simparams)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(POLAR_ORDER_DIRNAME, simparams), "$(FILE_NAME_PREFIX).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)
        open(outputtextfilefullpath, "w") do io
            println(io, "Polar Order Data")
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simparams))
            println(io, "Plot Data")
            writedlm(io, pomatrix, ",")
        end
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

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)
    plt.plot(pomatrix[1, :], pomatrix[2, :])
    plt.xlabel(PO_PLOT_XLABEL)
    plt.ylabel(PO_PLOT_YLABEL)
    plt.title("$(PO_PLOT_TITLE_TOPLINE)\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(POLAR_ORDER_DIRNAME, simparams)
        outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end

"""
Plots the data previously saved to a `.txt` file.

Can specify `minindex` and `maxindex` to only plot that many data points along the horizontal axis.
"""
function po_saveddata_plot(txtfilename::String, minindex::Int64=1, maxindex::Int64=0, saveplot::Bool=false, show::Bool=true)
    

    # get the plot data from the file
    xydatamat, sweeptype, simparamsout, legendvalues = po_saveddata_load(txtfilename)

    # handle default `maxindex` behavior
    if maxindex == 0
        maxindex = size(xydatamat)[2]
    end

    ### Plotting
    # prepare plot
    plt.clf()
    plt.grid(true, zorder=0)
    ###

    local paramsstr
    if sweeptype != nosweep
        filenameprefix = "polar_order_plot_$(sweeptype)"
        local legendlabels::Array{String}
        local legendtitle::String
        local sweepcm
        if sweeptype == densitysweep
            legendtitle = "Number of Particles"
            sweepcm = PO_LINE_DENSITY_COLORMAP
        elseif sweeptype == boxwidthsweep
            legendtitle = "Box Width"
            sweepcm = PO_LINE_BOXWIDTH_COLORMAP
        elseif sweeptype == interactionstrengthsweep
            legendtitle = "Interaction Fliprate"
            sweepcm = PO_LINE_INTERACTION_COLORMAP
        end
        
        legendlabels = (val -> "$val").(legendvalues)
        
        # get params string for plot
        paramsstr = simparamsout
        nshots = length(legendvalues)

        for i in 1:nshots            
            # plot each line one at a time
            linecolor = get(sweepcm, (i / nshots))
            plt.plot(xydatamat[1, minindex:maxindex], xydatamat[i+1, minindex:maxindex], color=(linecolor.r, linecolor.g, linecolor.b))
        end
        
        # add legends and labels to plot
        plt.legend(legendlabels, title=legendtitle)

    else
        filenameprefix = "polar_order_plot"

        # prepare strings for labels used in plot
        local interactionstr::String = ""
        if simparamsout.interaction == nointeraction
            interactionstr = "no-interaction  "
        else
            interactionstr = "$(simparamsout.interaction) I=$(Int64(round(simparamsout.interactionfliprate)))  "
        end
        
        # get params string for plot
        paramsstr = "N=$(simparamsout.numparticles)  Boxwidth=$(simparamsout.boxwidth)  t=$(Int64(simparamsout.totaltime))  $(interactionstr)"

        # plot the line
        plt.plot(xydatamat[1, minindex:maxindex], xydatamat[2, minindex:maxindex], color=PO_LINE_DEFAULT_COLOR)
    end

    plt.xlabel(PO_PLOT_XLABEL)
    plt.ylabel(PO_PLOT_YLABEL)
    plt.title("$(PO_PLOT_TITLE_TOPLINE)\n$(paramsstr)")        
    
    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)

        ndatapoints = size(xydatamat)[2]
        local rangestr::String = ""
        if minindex != 0 && maxindex != ndatapoints
            rangestr = "_indexrange_$minindex-$maxindex"
        elseif minindex == 1 && maxindex != ndatapoints
            rangestr = "_maxindex-$maxindex"
        elseif minindex != 1 && maxindex == ndatapoints
            rangestr = "_minindex-$minindex"
        end

        outputfilename = joinpath(datadirname, "$(filenameprefix)_from_txt$rangestr.pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end


"""
Loads the data previously saved to a `.txt` file.

Returns a tuple of 4 objects:


    `xydatamatrix`:   The matrix holding the analysis data. Row 0 holds the x-axis data.
                        The remaining lines hold the y-axis data
                        
    `sweeptype`:        The `SweepType` of the data. Is `nosweep` if not a sweep.
    
    `simparamsout`:     The simulation parameters used. Returns a `SimulationParameters` object if not a sweep.
                        Returns a string to use in the plot title if it is a sweep

    `legendvalues`:     The list of values of the swept parameter in the sweep. If the file is not for a sweep, contains `[nothing]`.
"""
function po_saveddata_load(txtfilename::String)
    
    # get the plot data from the file
    file = open(txtfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert contains(lowercase(line), "polar order data") "PO data file does not have correct first line. File may be corrupted."
    
    # check if it is a sweep
    line = readline(file)

    local paramsstr
    local podata::Array{Float64} = zeros(0, 0)
    local tdata
    local simparamsout
    local legendvalues

    if contains(lowercase(line), "sweep type")
        sweeptype = eval(Symbol(line[length("sweep type: ")+1:end]))

        local legendlabels::Array{String}
        local legendtitle::String
        # local paramsstr::String

        # now get sweep values
        line = readline(file)
        line = readline(file)
        
        legendlabels = split(line, ",")
        legendvalues = parse.(Float64, legendlabels)
        
        # get params string for plot
        line = readline(file)
        paramsstr = line[length("Params String: ")+1:end]
        simparamsout = paramsstr
       

    else  
        sweeptype = nosweep      
        # get params
        # read next line to get number of position data points
        simparaminfo_str = readline(file)
        simparams = SimulationParameters(simparaminfo_str)
        simparamsout = simparams

        legendvalues = [nothing]
    end

    # skip to data matrix
    while !contains(lowercase(line), "data matrix") && !contains(lowercase(line), "plot data")
        line = readline(file)
    end

     
    # read x data (t)
    dataline = readline(file)
    tdata = parse.(Float64, split(dataline, ","))
    
    # read y data (PO)
    nshots = length(legendvalues)
    podata = zeros(nshots, length(tdata))

    # load each line
    for i in 1:nshots
        dataline = readline(file)
        podata[i, :] = parse.(Float64, split(dataline, ","))
    end

    xydatamatrix = zeros(size(podata)[1] + 1, size(podata)[2])
    xydatamatrix[1, :] = tdata
    xydatamatrix[2:end, :] = podata

    return xydatamatrix, sweeptype, simparamsout, legendvalues
end


"""
Produces a line plot of the polar order as a function of the time interval `t`
of all the files in a given directory, and plots the lines on the same axes.

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function po_plot_dir(dirname::String, sweepname::String, sweeptype::SweepType, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
    
    filename_prefix = "polar_order_sweep_plot"
    
    # initialize data matrix
    local full_pomatrix::Matrix{Float64} = zeros(0,0)
    
    # load every file in the directory
    datafiles_all = readdir(dirname)
    datafiles = []

    # only look at .txt files
    for file in datafiles_all
        if endswith(file, ".txt")
            push!(datafiles, file)
        end
    end

    sp_array::Array{SimulationParameters} = []
    for (i, filename) in enumerate(datafiles)
        inputfilename = joinpath(dirname, filename)

        # add simparams to array
        initialstate = loadsim_nlines(inputfilename, 1, 1, rowwisetxt, true)
        push!(sp_array, initialstate.simparams)

        # run correlation function and add to plot
        # get the `SimulationData` from the file

        # get file info
        nsegments, serialized = getsimfiletype(inputfilename)

        if nsegments != 1
            println("More than 1 segment in file $filename. Found $nsegments segments")
            continue
        end

        # load the spins from the file
        spins, simparams = loadsimspins(inputfilename)

        # run polar order calculation on data
        pomatrix_i::Matrix{Float64} = polarorder(spins, simparams)

        # if `full_pomatrix` is empty, add x-axis data and first line of y-axis data
        if length(full_pomatrix) == 0
            full_pomatrix = zeros(length(datafiles) + 1, size(pomatrix_i)[2])

            full_pomatrix[1:2, :] = pomatrix_i
        else
            # else just add y-axis data
            full_pomatrix[i+1, :] = pomatrix_i[2,:]
        end
    end

    
    # generate legend labels and title

    local interactionstr::String = ""
    if sp_array[1].interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(sp_array[1].interaction) I=$(sp_array[1].interactionfliprate)"
    end
    
    local sweepvalues::Array{Real}
    local legendtitle::String
    local paramsstr::String
    local sortedorder::Vector{Int}
    local sweepcm
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
        sortedorder = sortperm([sp.numparticles for sp in sp_array])
        sweepvalues = (sp -> sp.numparticles).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
        sweepcm = PO_LINE_DENSITY_COLORMAP
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sortedorder = sortperm([sp.boxwidth for sp in sp_array])
        sweepvalues = (sp -> sp.boxwidth).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  $(interactionstr)"
        sweepcm = PO_LINE_BOXWIDTH_COLORMAP
    elseif sweeptype == interactionstrengthsweep
        interactionstr = "$(sp_array[1].interaction)"
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        sweepvalues = (sp -> sp.interactionfliprate).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
        sweepcm = PO_LINE_INTERACTION_COLORMAP
    end
    legendlabels::Array{String} = (val -> "$(val)").(sweepvalues)


    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysissweepdir(POLAR_ORDER_DIRNAME, sp_array[sortedorder], sweepname), "$(filename_prefix).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)

        open(outputtextfilefullpath, "w") do io
            println(io, "Polar Order Data")
            println(io, "Sweep Type: $(sweeptype)")
            println(io, "Sweep Values")
            writedlm(io, permutedims(legendlabels), ",")
            println(io, "Params String: $(paramsstr)")
            println(io, "Data Matrix")
            writedlm(io, permutedims(full_pomatrix[1, :]), ",")
            for i in sortedorder
                writedlm(io, permutedims(full_pomatrix[i+1, :]), ",")
            end
        end
    end





    #################
    ### Plotting ####
    #################

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)

    for i in 1:length(sortedorder)
        linecolor = get(sweepcm, sweepvalues[i] / sweepvalues[end])
        plt.plot(full_pomatrix[1,:], full_pomatrix[i+1,:], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    plt.legend(legendlabels, title=legendtitle)
    plt.xlabel(PO_PLOT_XLABEL)
    plt.ylabel(PO_PLOT_YLABEL)
    plt.title("$(PO_PLOT_TITLE_TOPLINE)\n$(paramsstr)")
    
    if saveplot
        # get dir path to save plot
        datadirname = getanalysissweepdir(POLAR_ORDER_DIRNAME, sp_array, sweepname)
        outputfilename = joinpath(datadirname, "$(filename_prefix).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end
