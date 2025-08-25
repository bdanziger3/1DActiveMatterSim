########################
# polarorder_driver.jl
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run polar order and mean spin as a function of time analysis on various data files and plot results
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


FONT = "Times New Roman"

# Use serif font in math text
rc("mathtext", fontset="cm")
PyPlot.rc("font", family="serif")           
PyPlot.rc("font", family=FONT)






"""
Produces a line plot of the mean spin as a function of the time

Set `savetxt=true` to save a .txt file with the mean spin data results
Set `show` to display the plot.
"""
function ms_plot(filename::String, saveplot::Bool=false, savetxt::Bool=false, show::Bool=true)
    
    # get the spin data from the file
    spins, simparams = loadsimspins(filename)

    # run mean spin calculation on data
    pomatrix::Matrix{Float64} = meanspin(spins, simparams)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(MEAN_SPIN_DIRNAME, simparams), "$(FILE_NAME_PREFIX).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)
        open(outputtextfilefullpath, "w") do io
            println(io, "Mean Spin Data")
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
    plt.xlabel(MS_PLOT_XLABEL)
    plt.ylabel(MS_PLOT_YLABEL)
    plt.title("$(MS_PLOT_TITLE_TOPLINE)\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(MEAN_SPIN_DIRNAME, simparams)
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
function ms_saveddata_plot(txtfilename::String, mintime::Real=0, maxtime::Real=-1, saveplot::Bool=false, show::Bool=true)

    # get the plot data from the file
    xydatamat, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

    # handle `maxtime` default value
    if maxtime == -1
        maxtime = xydatamat[1, end]
    end

    ### Plotting
    # prepare plot
    plt.clf()
    ###

    local paramsstr
    if sweeptype != nosweep
        filenameprefix = "mean_spin_plot_$(sweeptype)"
        local legendlabels::Array{String}
        local legendtitle::String
        local sweepcm
        if sweeptype == densitysweep
            legendtitle = "Number of Particles"
            sweepcm = MS_LINE_DENSITY_COLORMAP
        elseif sweeptype == boxwidthsweep
            legendtitle = "Box Width"
            sweepcm = MS_LINE_BOXWIDTH_COLORMAP
        elseif sweeptype == interactionstrengthsweep
            legendtitle = "Interaction Fliprate"
            sweepcm = MS_LINE_INTERACTION_COLORMAP
        end
        
        legendlabels = (val -> "$val").(legendvalues)
        
        # get params string for plot
        paramsstr = simparamsout
        nshots = length(legendvalues)

        for i in 1:nshots            
            # plot each line one at a time
            linecolor = get(sweepcm, (i / nshots))
            plt.plot(xydatamat[1, :], xydatamat[i+1, :], color=(linecolor.r, linecolor.g, linecolor.b))
        end
        
        # add legends and labels to plot
        plt.legend(legendlabels, title=legendtitle)

    else
        filenameprefix = "mean_spin_plot"

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
        plt.plot(xydatamat[1, :], xydatamat[2, :], color=MS_LINE_DEFAULT_COLOR)
    end

    plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
    plt.ylabel(L"$\langle s_i \rangle $", fontsize=16, fontname=FONT)

    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")
    # plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(simparamsout)")
    
    ax = plt.gca()
    for label in cat(ax.get_xticklabels(), ax.get_yticklabels(), dims=1)
        label.set_fontname(FONT)  # font family
    end
    
    plt.xlim([mintime, maxtime])
    # plt.title("$(MS_PLOT_TITLE_TOPLINE)\n$(paramsstr)")        
    
    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)

        actualmaxtime = xydatamat[1, end]
        local rangestr::String = ""
        if mintime != 0 && maxtime != actualmaxtime
            rangestr = "_trange_$(Int64(round(mintime)))-$(Int64(round(maxtime)))"
        elseif mintime == 1 && maxtime != actualmaxtime
            rangestr = "_maxtime-$(Int64(round(maxtime)))"
        elseif mintime != 1 && maxtime == actualmaxtime
            rangestr = "_mintime-$(Int64(round(mintime)))"
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
Plots the data previously saved to a `.txt` file to a different plot for each shot

Can specify `minindex` and `maxindex` to only plot that many data points along the horizontal axis.
"""
function ms_saveddata_plot_individual(txtfilename::String, mintime::Real=0, maxtime::Real=-1, saveplot::Bool=false, show::Bool=true)

    # get the plot data from the file
    xydatamat, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

    # handle `maxtime` default value
    if maxtime == -1
        maxtime = xydatamat[1, end]
    end

    ### Plotting
    # prepare plot
    plt.clf()
    ###

    local paramsstr
    filenameprefix = "mean_spin_plot_$(sweeptype)"
    local legendlabels::Array{String}
    local legendtitle::String
    local sweepcm
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
        sweepcm = MS_LINE_DENSITY_COLORMAP
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sweepcm = MS_LINE_BOXWIDTH_COLORMAP
    elseif sweeptype == interactionstrengthsweep
        legendtitle = "Interaction Fliprate"
        sweepcm = MS_LINE_INTERACTION_COLORMAP
    end
    
    legendlabels = (val -> "$val").(legendvalues)
    
    # get params string for plot
    paramsstr = simparamsout
    nshots = length(legendvalues)

    for i in 1:nshots  
        plt.clf()
        plt.plot(xydatamat[1, :], xydatamat[i+1, :], color=MS_LINE_DEFAULT_COLOR)
    
    
        # plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
        plt.ylabel(L"$\langle s_i \rangle $", fontsize=16, fontname=FONT)

        plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both", labelsize=16, pad = 10)
        plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both", labelsize=16, pad = 10)
        # plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(simparamsout)")
    
        plt.xticks([450, 500])
        plt.yticks([-1, 0, 1])

        ax = plt.gca()
        for label in cat(ax.get_xticklabels(), ax.get_yticklabels(), dims=1)
            label.set_fontname(FONT)  # font family
        end
    
        plt.xlim([mintime, maxtime])
        plt.ylim([-1, 1])
        # plt.title("$(MS_PLOT_TITLE_TOPLINE)\n$(paramsstr)")        
        
        if saveplot
            # get dir path to save plot
            datadirname = dirname(txtfilename)

            actualmaxtime = xydatamat[1, end]
            local rangestr::String = ""
            if mintime != 0 && maxtime != actualmaxtime
                rangestr = "_trange_$(Int64(round(mintime)))-$(Int64(round(maxtime)))"
            elseif mintime == 1 && maxtime != actualmaxtime
                rangestr = "_maxtime-$(Int64(round(maxtime)))"
            elseif mintime != 1 && maxtime == actualmaxtime
                rangestr = "_mintime-$(Int64(round(mintime)))"
            end

            outputfilename = joinpath(datadirname, "$(filenameprefix)_from_txt_$(legendlabels[i])$rangestr.pdf")
            outputfilename = checkfilename(outputfilename)
            plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
        end
        
        if show
            plt.show()
        end
    end
end


"""
Plots the data previously saved to a `.txt` file to a different plot for each shot

Can specify `minindex` and `maxindex` to only plot that many data points along the horizontal axis.
"""
function ms_saveddata_plot_individual_subplots(txtfilename::String, mintime::Real=0, maxtime::Real=-1, saveplot::Bool=false, show::Bool=true)

    # get the plot data from the file
    xydatamat, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

    # handle `maxtime` default value
    if maxtime == -1
        maxtime = xydatamat[1, end]
    end

    ### Plotting
    # prepare plot
    ###

    local paramsstr
    filenameprefix = "mean_spin_plot_$(sweeptype)"
    local legendlabels::Array{String}
    local legendtitle::String
    local sweepcm
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
        sweepcm = MS_LINE_DENSITY_COLORMAP
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sweepcm = MS_LINE_BOXWIDTH_COLORMAP
    elseif sweeptype == interactionstrengthsweep
        legendtitle = "Interaction Fliprate"
        sweepcm = MS_LINE_INTERACTION_COLORMAP
    end
    
    legendlabels = (val -> "$val").(legendvalues)
    
    # get params string for plot
    paramsstr = simparamsout
    nshots = length(legendvalues)

    fig, axs = plt.subplots(2, 3, figsize=(9, 6), sharex="col", sharey="row")

    for i in 1:2
        for j in 1:3
            ax = axs[i, j]
            realindex = j + (3 * (i-1))
            ax.plot(xydatamat[1, :], xydatamat[1+realindex, :], color=MS_LINE_DEFAULT_COLOR, linewidth=1)
            # ax.text(400, 0.8, raw"$\r_{\text{align}} = 100$")
            ax[:text](0.95, 0.95, "$(raw"$r_{align} = $") $(Int64(legendvalues[realindex]))", transform=ax[:transAxes], ha="right", va="top", fontsize=12)

            ax.set_xlim(mintime, maxtime)
            ax.set_ylim(-1, 1)

            ax.set_xticks([mintime, maxtime])
            ax.set_yticks([-1, 0, 1])

            ax[:tick_params](axis="x", top=true, bottom=true, direction="in", which = "both", labelsize=12, pad = 10)
            ax[:tick_params](axis="y", left=true, right=true, direction="in", which = "both", labelsize=12, pad = 10)

            for label in cat(ax.get_xticklabels(), ax.get_yticklabels(), dims=1)
                label.set_fontname(FONT)  # font family
            end

            # Only label the leftmost y-axis
            if j == 1
                ax.set_ylabel(L"$\langle s_i \rangle $", fontsize=16, fontname=FONT)
            end
            
            # # Only label the bottom x-axis
            if i == 2
                ax.set_xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
            end
        end
    end      
        
    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)

        actualmaxtime = xydatamat[1, end]
        local rangestr::String = ""
        if mintime != 0 && maxtime != actualmaxtime
            rangestr = "_trange_$(Int64(round(mintime)))-$(Int64(round(maxtime)))"
        elseif mintime == 1 && maxtime != actualmaxtime
            rangestr = "_maxtime-$(Int64(round(maxtime)))"
        elseif mintime != 1 && maxtime == actualmaxtime
            rangestr = "_mintime-$(Int64(round(mintime)))"
        end

        outputfilename = joinpath(datadirname, "$(filenameprefix)_subplots_from_txt$rangestr.pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end

    fig[:tight_layout]()
        
    if show
        plt.show()
    end
end


"""
Calculate Mean Spin of all files in a directory and save data of all of them in a single .txt file
"""
function ms_savetxt_dir(dirpath::String, sweepname::String, sweeptype::SweepType)
    filename_prefix = "ms_sweep_plot"
    
    # initialize data matrix
    local full_msmatrix::Matrix{Float64} = zeros(0,0)
    
    # load every file in the directory
    datafiles_all = readdir(dirpath)
    datafiles = []

    # only look at .txt files
    for file in datafiles_all
        if endswith(file, ".txt")
            push!(datafiles, file)
        end
    end

    sp_array::Array{SimulationParameters} = []
    for (i, filename) in enumerate(datafiles)
        inputfilename = joinpath(dirpath, filename)
        
        # find mean spins function and add to plot
        # get the position data from the file
        spins, simparams = loadsimspins(inputfilename)
        
        # add simparams to array
        push!(sp_array, simparams)
        
        # run mean spin calculation on data
        msmat_i::Matrix{Float64} = meanspin(spins, simparams)

        # if `full_msmatrix` is empty, add x-axis data and first line of y-axis data
        if length(full_msmatrix) == 0
            full_msmatrix = zeros(length(datafiles) + 1, size(msmat_i)[2])

            full_msmatrix[1:2, :] = msmat_i
        else
            # else just add y-axis data
            full_msmatrix[i+1, :] = msmat_i[2,:]
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
        sweepcm = MS_LINE_DENSITY_COLORMAP
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sortedorder = sortperm([sp.boxwidth for sp in sp_array])
        sweepvalues = (sp -> sp.boxwidth).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  $(interactionstr)"
        sweepcm = MS_LINE_BOXWIDTH_COLORMAP
    elseif sweeptype == interactionstrengthsweep
        interactionstr = "$(sp_array[1].interaction)"
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        sweepvalues = (sp -> sp.interactionfliprate).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
        sweepcm = MS_LINE_INTERACTION_COLORMAP
    end
    legendlabels::Array{String} = (val -> "$(val)").(sweepvalues)


    # save as a datafile if requested
    outputtextfilefullpath = joinpath(getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array[sortedorder], sweepname), "$(filename_prefix).txt")
    outputtextfilefullpath = checkfilename(outputtextfilefullpath)

    open(outputtextfilefullpath, "w") do io
        println(io, "Mean Spin Data")
        println(io, "Sweep Type: $(sweeptype)")
        println(io, "Sweep Values")
        writedlm(io, permutedims(legendlabels), ",")
        println(io, "Params String: $(paramsstr)")
        println(io, "Data Matrix")
        writedlm(io, permutedims(full_msmatrix[1, :]), ",")
        for i in sortedorder
            writedlm(io, permutedims(full_msmatrix[i+1, :]), ",")
        end
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
function ms_saveddata_load(txtfilename::String)
    
    # get the plot data from the file
    file = open(txtfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert contains(lowercase(line), "mean spin data") "MS data file does not have correct first line. File may be corrupted."
    
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
    
    # read y data (MS)
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
Produces a line plot of the mean spin as a function of the time interval `t`
of all the files in a given directory, and plots the lines on the same axes.

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function ms_plot_dir(dirname::String, sweepname::String, sweeptype::SweepType, saveplot::Bool=true, savetxt::Bool=true, individualplots::Bool=false, show::Bool=false)
    
    filename_prefix = "mean_spin_sweep_plot"
    
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
        
        # run correlation function and add to plot
        # get the spins from the file
        
        # load the spins from the file
        spins, simparams = loadsimspins(inputfilename)
        # add simparams to array
        push!(sp_array, simparams)

        # run mean spin calculation on data
        pomatrix_i::Matrix{Float64} = meanspin(spins, simparams)

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
        sweepcm = MS_LINE_DENSITY_COLORMAP
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sortedorder = sortperm([sp.boxwidth for sp in sp_array])
        sweepvalues = (sp -> sp.boxwidth).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  $(interactionstr)"
        sweepcm = MS_LINE_BOXWIDTH_COLORMAP
    elseif sweeptype == interactionstrengthsweep
        interactionstr = "$(sp_array[1].interaction)"
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        sweepvalues = (sp -> sp.interactionfliprate).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
        sweepcm = MS_LINE_INTERACTION_COLORMAP
    end
    legendlabels::Array{String} = (val -> "$(val)").(sweepvalues)


    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array[sortedorder], sweepname), "$(filename_prefix).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)

        open(outputtextfilefullpath, "w") do io
            println(io, "Mean Spin Data")
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

    # if `individualplots`, plot the data seperately
    if individualplots
        for (i, s) in enumerate(sortedorder)        
            # clear plot
            plt.clf()

            plt.grid(true, zorder=0)
            plt.plot(full_pomatrix[1, :], full_pomatrix[s+1, :])
            plt.xlabel(MS_PLOT_XLABEL)
            plt.ylabel(MS_PLOT_YLABEL)
            plt.title("$(MS_PLOT_TITLE_TOPLINE)\n N=$(sp_array[s].numparticles)  Boxwidth=$(sp_array[s].boxwidth)  t=$(Int64(sp_array[s].totaltime))  $(sp_array[1].interaction) I=$(sp_array[s].interactionfliprate)")

            # figure out swept value and var and add to names
            datadirname = getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array, sweepname)
            outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX)_$(sweepvalues[i]).pdf")
            outputfilename = checkfilename(outputfilename)
            plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
        end
    end

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)

    for i in 1:length(sortedorder)
        linecolor = get(sweepcm, sweepvalues[i] / sweepvalues[end])
        plt.plot(full_pomatrix[1,:], full_pomatrix[i+1,:], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    plt.legend(legendlabels, title=legendtitle)
    plt.xlabel(MS_PLOT_XLABEL)
    plt.ylabel(MS_PLOT_YLABEL)
    plt.title("$(MS_PLOT_TITLE_TOPLINE)\n$(paramsstr)")
    
    if saveplot
        # get dir path to save plot
        datadirname = getanalysissweepdir(MEAN_SPIN_DIRNAME, sp_array, sweepname)
        outputfilename = joinpath(datadirname, "$(filename_prefix).pdf")
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
    xydatamat, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

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
    local xaxistitle
    if sweeptype != nosweep
        filenameprefix = "polar_order_plot_$(sweeptype)"
        if sweeptype == densitysweep
            xaxistitle = "Number of Particles"
            sweepcm = MS_LINE_DENSITY_COLORMAP
        elseif sweeptype == boxwidthsweep
            xaxistitle = "Box Width"
            sweepcm = MS_LINE_BOXWIDTH_COLORMAP
        elseif sweeptype == interactionstrengthsweep
            xaxistitle = "Interaction Fliprate"
            sweepcm = MS_LINE_INTERACTION_COLORMAP
        end
                
        # get params string for plot
        paramsstr = simparamsout
        nshots = length(legendvalues)

        avgpo::Vector{Float64} = zeros(length(legendvalues))
        for i in 1:nshots
            # calculate average PO for the given range
            avgpo[i] = sum(abs.(xydatamat[i+1, minindex:maxindex]))[1] / (maxindex - minindex + 1)
        end

        # plot results
        linecolor = get(sweepcm, 1.0)
        plt.plot(legendvalues, avgpo, "o")
        axes = plt.gca()
        axes.loglog()

    else
        filenameprefix = "polar_order_plot"

        # prepare strings for labels used in plot
        local interactionstr::String = ""
        if simparamsout.interaction == nointeraction
            interactionstr = "no-interaction  "
        else
            interactionstr = "$(simparamsout.interaction) I=$(Int64(round(simparamsout.interactionfliprate)))  "
        end

        absxymat = abs.(xydatamat[2, minindex:maxindex])
        po_sum = sum(absxymat)
        avgpo = po_sum / (maxindex - minindex + 1)
        println(avgpo)

        
        
    end


    plt.xlabel(MS_PLOT_XLABEL)
    plt.ylabel(PO_PLOT_YLABEL)
    plt.title("$(PO_PLOT_TITLE_TOPLINE)\n$(paramsstr)")        
    
    # if saveplot
    #     # get dir path to save plot
    #     datadirname = dirname(txtfilename)

    #     ndatapoints = size(xydatamat)[2]
    #     local rangestr::String = ""
    #     if minindex != 0 && maxindex != ndatapoints
    #         rangestr = "_indexrange_$minindex-$maxindex"
    #     elseif minindex == 1 && maxindex != ndatapoints
    #         rangestr = "_maxindex-$maxindex"
    #     elseif minindex != 1 && maxindex == ndatapoints
    #         rangestr = "_minindex-$minindex"
    #     end

    #     outputfilename = joinpath(datadirname, "$(filenameprefix)_from_txt$rangestr.pdf")
    #     outputfilename = checkfilename(outputfilename)
    #     plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    # end
    
    if show
        plt.show()
    end
end

function po_avg_sweep_txt(txtfilename::String, settletime::Real, saveplot::Bool=false, show::Bool=true)
    # load data
    datamatrix, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)

    # handle settletime
    settleindex = Int64(round(settletime / datamatrix[1,2])) + 1

    # average over cases with same sweep value
    # find the unique sweep values
    sweepvaluecounts = countmap(legendvalues)
    uniquesweepvalues = sort(collect(keys(sweepvaluecounts)))
    nvalues = length(uniquesweepvalues)

    indivpo = [[] for _=1:nvalues]

    averagedpo::Array{Float64} = zeros(nvalues)

    for (i, val) in enumerate(legendvalues)
        legendvalue_uniqueindex = findfirst(uniquesweepvalues .== val)
        meanpo = mean(abs.(datamatrix[i+1, settleindex:end]))
        averagedpo[legendvalue_uniqueindex] += (meanpo ./ sweepvaluecounts[val])
        
        push!(indivpo[legendvalue_uniqueindex], meanpo) 
    end


    plt.plot(uniquesweepvalues, averagedpo, "-o")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.show()
end


function ms_plot_logrevtimes(txtfilename::String, settletime::Real=100, saveplot::Bool=false, show::Bool=true)
    datamatrix, sweeptype, simparamsout, legendvalues = ms_saveddata_load(txtfilename)
    
    revtimes2, revtimes = getfftreversaltimes(datamatrix, settletime)


    # average data with same sweep value
    # find the unique sweep values
    sweepvaluecounts = countmap(legendvalues)
    uniquesweepvalues = sort(collect(keys(sweepvaluecounts)))
    nvalues = length(uniquesweepvalues)

    infivrevs = [[] for _=1:nvalues]

    averagedrevtimes::Array{Float64} = zeros(nvalues)

    for (i, val) in enumerate(legendvalues)
        legendvalue_uniqueindex = findfirst(uniquesweepvalues .== val)
        averagedrevtimes[legendvalue_uniqueindex] += (revtimes[i] ./ sweepvaluecounts[val])


        
        push!(infivrevs[legendvalue_uniqueindex], revtimes[i]) 
    end


    
    avgrevs = (l -> mean(l)).(infivrevs)
    stderrrev = (l -> std(l) / sqrt(length(l))).(infivrevs)

    interval95p = 1.96

    # find linear fit of loglog plot
    criticalval = 2
    criticalindex = findmin(abs.(uniquesweepvalues .- criticalval))[2]
    
    # params, stderrs, linearmodel = logloglinearfit(uniquesweepvalues[1:criticalindex], averagedrevtimes[1:criticalindex])
    params2, stderrs2, linearmodel2 = logloglinearfit(uniquesweepvalues[criticalindex:end], averagedrevtimes[criticalindex:end])

    # println(params)
    # println(stderrs)

    println(params2)
    println(stderrs2)

    # plot line
    # remove 0 from data
    cleanedsweepvals = uniquesweepvalues[uniquesweepvalues .!= 0]
    cleanedrevtimes = averagedrevtimes[uniquesweepvalues .!= 0]
    cleanedstderrrev = stderrrev[uniquesweepvalues .!= 0]
    linex = collect(minimum(cleanedsweepvals):1:maximum(cleanedsweepvals))
    # liney = linearmodel(linex)
    liney2 = linearmodel2(linex)
    
    # plot fit lines
    plt.clf()
    # plt.plot(linex, liney, "--", c="k", linewidth=0.5)
    plt.plot(linex, liney2, "--", c="k", linewidth=0.5)
   
    # Plotting

    marker = "s"
    markercolor = "teal"

    plt.plot(cleanedsweepvals, cleanedrevtimes, marker, markeredgecolor=markercolor,  markerfacecolor="none", markersize=5)
    plt.errorbar(cleanedsweepvals, cleanedrevtimes, yerr=interval95p.*cleanedstderrrev, fmt="none", elinewidth=0.25, ecolor="k", barsabove=false, capsize=1.5)
    
    plt.yscale("log")
    plt.xscale("log")

    plt.xlabel(L"$r_{\text{align}}$", fontsize=24, fontname=FONT)
    plt.ylabel(L"$\tau_{\text{rev}}$", fontsize=28, fontname=FONT)

    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")
    


    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        outputfilename = joinpath(datadirname, "ms_revtime_loglog_plot_from_txt.pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end

    if show
        plt.show()
    end

    return uniquesweepvalues, averagedrevtimes

end


# datatxt = fixpath("work/analysis/Mean Spin/Align Simple/compare_ms_plots/ms_sweep_plot.txt")

# dir1 = fixpath("work/data/sweeps/alignsimple/interactionsweep/Aug15-B50-Isweep/compare")

# ms_savetxt_dir(dir1, "compare_ms_plots", interactionstrengthsweep)
# dir2 = fixpath("work/data/sweeps/alignsimple/interactionsweep/Aug15-B50-Isweep/firsthalf")
# ms_savetxt_dir(dir2, "Aug15_rsweep_half1", interactionstrengthsweep)

# file = fixpath("work/analysis/Mean Spin/Align Simple/Aug15_rsweep_half2/ms_sweep_plot_test.txt")

# ms_saveddata_plot_individual(datatxt, 450, 500, true, false)
# ms_saveddata_plot_individual_subplots(datatxt, 450, 500, true, true)

# filetxt1 = fixpath("work/analysis/Mean Spin/Align Simple/Aug15_rsweep/Aug15_rsweep_full/ms_sweep_plot.txt")
# filetxt2 = fixpath("work/analysis/Mean Spin/Align Simple/small_interaction2_3-8/mean_spin_sweep_plot.txt")


# d1x, d1y = ms_plot_logrevtimes(filetxt1, 100, true, true)
# d2x, d2y = ms_plot_logrevtimes(filetxt2, 100, true, true)

f1 = fixpath("work/analysis/Mean Spin/Align Simple/Aug15_rsweep/Aug15_rsweep_full/ms_sweep_plot.txt")
po_avg_sweep_txt(f1, 100)

# marker = "s"
# plt.plot(d1x, d1y, marker, markeredgecolor="firebrick",  markerfacecolor="none", markersize=5)
# plt.plot(d2x, d2y, marker, markeredgecolor="forestgreen",  markerfacecolor="none", markersize=5)

# plt.yscale("log")
# plt.xscale("log")

# plt.xlabel(L"$\rho$", fontsize=24, fontname=FONT)
# plt.ylabel(L"$\tau_{\text{rev}}$", fontsize=28, fontname=FONT)

# plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
# plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")

# plt.show()
