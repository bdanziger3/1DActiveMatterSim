########################
# orientationcorrelation_driver.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run orientation self-correlation function on variou data files and plot results
########################

using PyPlot
using PyCall
using LsqFit
using ColorSchemes
using StatsBase

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../../file_management/analysis_files.jl")

MSD_DIRNAME = "Mean Squared Displacement"

FILE_NAME_PREFIX = "msd"

MSD_LINE_DENSITY_COLORMAP = ColorSchemes.Reds_8
MSD_LINE_INTERACTION_COLORMAP = ColorSchemes.Blues_8
MSD_LINE_BOXWIDTH_COLORMAP = ColorSchemes.Greens_8

# MSD_LINE_COLORMAP = ColorSchemes.vik100
# MSD_LINE_COLORMAP = ColorSchemes.Blues_8
MSD_LINE_COLORMAP = ColorSchemes.bam50

FONT = "Times New Roman"

# Use serif font in math text
rc("mathtext", fontset="cm")
PyPlot.rc("font", family="serif")           
PyPlot.rc("font", family=FONT)


"""
Calculate MSD of a single simulation and save data as a .txt file
"""
function msd_savetxt(filename::String, settletime::Real=-1.0)
    # get the position data from the file
    positions, simparams = loadsimpositions(filename)

    msdmat::Matrix{Float64} = meansqdisp(positions, simparams, settletime)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(MSD_DIRNAME, simparams), "$(FILE_NAME_PREFIX)-$(settletime == -1 ? "all" : settletime).txt")

        open(outputtextfilefullpath, "w") do io
            println(io, "MSD Data")
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simparams))
            println(io, "Plot Data")
            writedlm(io, msdmat, ",")
        end
    end
end

"""
Calculate MSD of all files in a directory and save data of all of them in a single .txt file

"""
function msd_savetxt_dir(dirpath::String, sweepname::String, sweeptype::SweepType, settletime::Real=-1.0, maxdt::Real=20)
    filename_prefix = "msd_sweep_plot"
    
    # initialize data matrix
    local full_msdmatrix::Matrix{Float64} = zeros(0,0)
    
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
        
        # run MSD function and add to plot
        # get the position data from the file
        positions, simparams = loadsimpositions(inputfilename)
        
        # add simparams to array
        push!(sp_array, simparams)
        
        # run mean spin calculation on data
        msdmat_i::Matrix{Float64} = meansqdisp(positions, simparams, settletime, maxdt)

        # if `full_msdmatrix` is empty, add x-axis data and first line of y-axis data
        if length(full_msdmatrix) == 0
            full_msdmatrix = zeros(length(datafiles) + 1, size(msdmat_i)[2])

            full_msdmatrix[1:2, :] = msdmat_i
        else
            # else just add y-axis data
            full_msdmatrix[i+1, :] = msdmat_i[2,:]
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
        sweepcm = MSD_LINE_DENSITY_COLORMAP
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sortedorder = sortperm([sp.boxwidth for sp in sp_array])
        sweepvalues = (sp -> sp.boxwidth).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  $(interactionstr)"
        sweepcm = MSD_LINE_BOXWIDTH_COLORMAP
    elseif sweeptype == interactionstrengthsweep
        interactionstr = "$(sp_array[1].interaction)"
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        sweepvalues = (sp -> sp.interactionfliprate).(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
        sweepcm = MSD_LINE_INTERACTION_COLORMAP
    end
    legendlabels::Array{String} = (val -> "$(val)").(sweepvalues)


    # save as a datafile if requested
    outputtextfilefullpath = joinpath(getanalysissweepdir(MSD_DIRNAME, sp_array[sortedorder], sweepname), "$(filename_prefix).txt")
    outputtextfilefullpath = checkfilename(outputtextfilefullpath)

    open(outputtextfilefullpath, "w") do io
        println(io, "MSD Data")
        println(io, "Sweep Type: $(sweeptype)")
        println(io, "Sweep Values")
        writedlm(io, permutedims(legendlabels), ",")
        println(io, "Params String: $(paramsstr)")
        println(io, "Data Matrix")
        writedlm(io, permutedims(full_msdmatrix[1, :]), ",")
        for i in sortedorder
            writedlm(io, permutedims(full_msdmatrix[i+1, :]), ",")
        end
    end
end

"""
Produces a line plot of the mean squared displacement as a function of the time interval `dt`

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function msd_plot(filename::String, settletime::Real=-1.0, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
    
    # get the position data from the file
    positions, simparams = loadsimpositions(filename)

    msdmat::Matrix{Float64} = meansqdisp(positions, simparams, settletime)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(MSD_DIRNAME, simparams), "$(FILE_NAME_PREFIX)-$(settletime == -1 ? "all" : settletime).txt")

        open(outputtextfilefullpath, "w") do io
            println(io, "Mean Squared Displacement Data")
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simparams))
            println(io, "Plot Data")
            writedlm(io, msdmat, ",")
        end
    end


    local interactionstr::String = ""
    if simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(simparams.interaction) I=$(Int64(round(simparams.interactionfliprate)))  "
    end


    # more file naming controls


    ### Plotting

    # clear plot
    plt.clf()
    
    plt.plot(msdmat[1,:], msdmat[2,:])
    plt.xlabel("$(raw"Time Interval $dt$")")
    plt.ylabel("Mean Squared Displacement")
    plt.title("Mean Squared Displacement of Active Particle Simulation\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(round(simparams.totaltime)))  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(MSD_DIRNAME, simparams)
        plt.savefig(joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf"), bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end

"""
loads msd data from saved .txt data file
"""
function msd_loadtxt(txtfilename::String)
    
    # get the plot data from the file
    file = open(txtfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert contains(lowercase(line), "msd") "MSD data file does not have correct first line. File may be corrupted."
    
    # check if it is a sweep
    line = readline(file)

    local simparamsout
    local nshots
    if contains(lowercase(line), "sweep type")
        sweeptype = eval(Symbol(line[length("sweep type: ")+1:end]))
        local legendlabels::Array{String}
        local legendtitle::String
        local paramsstr::String
        if sweeptype == densitysweep
            legendtitle = "Number of Particles"
        elseif sweeptype == boxwidthsweep
            legendtitle = "Box Width"
        elseif sweeptype == interactionstrengthsweep
            legendtitle = "Interaction Fliprate"
        end

        # now get sweep values
        line = readline(file)
        line = readline(file)
        
        legendlabels = split(line, ",")
        legendvalues = parse.(Float64, legendlabels)
        nshots = length(legendlabels)
        
        # get params string for plot
        line = readline(file)
        paramsstr = line[length("Params String: ")+1:end]

        simparamsout = paramsstr
    else
        sweeptype = nosweep
        nshots = 1   
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


    # read x data (dt)
    dataline = readline(file)
    dtdata = parse.(Float64, split(dataline, ","))

    # initialize full matrix
    xydatamatrix::Matrix{Float64} = zeros(nshots+1, length(dtdata))
    xydatamatrix[1, :] = dtdata

    # read y data (OC)
    for i in 1:nshots
        dataline = readline(file)
        xydatamatrix[i+1, :] = parse.(Float64, split(dataline, ","))
    end
    
    return xydatamatrix, sweeptype, legendvalues, simparamsout
end


"""
Plots the data previously saved to a txt file.

Can specify `maxdt` to only plot up to that value of dt.
"""
function msd_plot_txt(txtfilename::String, maxdt::Real=-1., saveplot::Bool=false, show::Bool=true, valuestoplot=nothing)
    xydatamatrix, sweeptype, legendvalues, paramsstr = msd_loadtxt(txtfilename)

    # if no `maxdt` provided, defaults to maximum in data
    if maxdt == -1
        maxdt = xydatamatrix[1, end]
    end

    # find the unique sweep values
    sweepvaluecounts = countmap(legendvalues)
    uniquesweepvalues = sort(collect(keys(sweepvaluecounts)))
    nvalues = length(uniquesweepvalues)


    # average data with same sweep value
    averagedxydatamat::Matrix{Float64} = zeros(nvalues+1, size(xydatamatrix)[2])
    averagedxydatamat[1, :] = xydatamatrix[1, :]


    # TODO Calculate characteristic times for each individual sim
    # indivtaus = [[] for _=1:nvalues]


    for (i, val) in enumerate(legendvalues)
        legendvalue_uniqueindex = findfirst(uniquesweepvalues .== val)
        averagedxydatamat[legendvalue_uniqueindex+1, :] += (xydatamatrix[i+1, :] ./ sweepvaluecounts[val])


        # tau_i, _ = fitdecay(xydatamatrix[1, 1:fit_samples], xydatamatrix[i+1, 1:fit_samples])
        # push!(indivtaus[legendvalue_uniqueindex], tau_i) 
    end

    # handle`valuestoplot`
    if isnothing(valuestoplot)
        valuestoplot = uniquesweepvalues
    end

    ### Plotting
    # prepare plot
    plt.clf()
    ###

    linecolors = []

    linesplotted = 0
    for (i, val) in enumerate(sort(intersect(uniquesweepvalues, valuestoplot)))
        # find index of run in big matrix
        indexinmat = findfirst(uniquesweepvalues .== val)
        
        # plot each line one at a time
        positivehalf = false     # using positive half of the color scheme?
        # if uniquesweepvalues[i] in valuestoplot
        linesplotted += 1
        if positivehalf
            # linecolor = get(OC_LINE_COLORMAP, 0.5 + (0.5 * (val / maximum(intersect(uniquesweepvalues, valuestoplot)))))
            linecolor = get(MSD_LINE_COLORMAP, 0.5 + (0.5 * (i / length(intersect(uniquesweepvalues, valuestoplot)))))
        else
            # linecolor = get(OC_LINE_COLORMAP, 0.5 - (0.5 * (i / nvalues)))
            linecolor = get(MSD_LINE_COLORMAP, 0.5 - (0.5 * (i / length(intersect(uniquesweepvalues, valuestoplot)))))
        end

        push!(linecolors, linecolor)

        # TODO Calculate characteristic times for averaged data
            # tau, tau_stderr = fitdecay(averagedxydatamat[1, 1:fit_samples], averagedxydatamat[indexinmat+1, 1:fit_samples])

        plt.plot(averagedxydatamat[1, :], averagedxydatamat[1+indexinmat, :], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    if maxdt > 0
        plt.xlim([0, maxdt])
    end

    # add legends and labels to plot
    local legendlabels
    local alllegendvalues
    if sweeptype == densitysweep
        # convert from N to \rho
        legendtitle = L"$\rho$"
        boxwidth = 50
        legendlabels = round.(intersect(uniquesweepvalues, valuestoplot) ./ boxwidth, digits=1)
        alllegendvalues = uniquesweepvalues ./ boxwidth
    elseif sweeptype == interactionstrengthsweep
        legendtitle = raw"$r_{\text{align}}$"
        legendlabels = intersect(uniquesweepvalues, valuestoplot)
        alllegendvalues = uniquesweepvalues
    end

    plt.legend(legendlabels, title=legendtitle, title_fontsize=16, prop=Dict("family" => FONT, "size" => 8), frameon=false, handlelength=1.0, ncol=3)
    plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
    plt.ylabel(raw"MSD($t$)", fontsize=16, fontname=FONT)
    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")
    

    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        plotfilename = checkfilename(joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf"))
        plt.savefig(plotfilename, bbox_inches = "tight", pad_inches=0.1)
    end

    if show
        plt.show()
    end
    
    # plot on loglog axes data
    # clean data for log plot by removing 0s
    cleandtdata = averagedxydatamat[1, averagedxydatamat[1, :] .!= 0]
    cleanmsddata = averagedxydatamat[2:end, averagedxydatamat[1, :] .!= 0]

    plt.clf()


    linesplotted = 0
    for (i, val) in enumerate(sort(intersect(uniquesweepvalues, valuestoplot)))
        linesplotted += 1
        # find index of run in big matrix
        indexinmat = findfirst(uniquesweepvalues .== val)
        linecolor = linecolors[linesplotted]
        plt.plot(cleandtdata, cleanmsddata[i, :], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    plt.legend(legendlabels, title=legendtitle, title_fontsize=16, prop=Dict("family" => FONT, "size" => 8), frameon=false, handlelength=1.0, ncol=3)
    plt.xlim([cleandtdata[1], maxdt])
    axes = plt.gca()
    axes.loglog()
    plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
    plt.ylabel(raw"MSD($t$)", fontsize=16, fontname=FONT)
    
    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")
    
    ax = plt.gca()
    for label in cat(ax.get_xticklabels(), ax.get_yticklabels(), dims=1)
        label.set_fontname(FONT)  # font family
    end

    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        plotfilename = checkfilename(joinpath(datadirname, "$(FILE_NAME_PREFIX)_loglog.pdf"))
        plt.savefig(plotfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end


    # calculate slope of logplots
    derivlog = (log10.(cleanmsddata[:, 2:end]) .- log10.(cleanmsddata[:, 1:end-1])) ./ permutedims((log10.(cleandtdata[2:end]) .- log10.(cleandtdata[1:end-1])))
    dtmidpoints = (cleandtdata[1:end-1] .+ cleandtdata[2:end]) ./ 2

    ### Plotting

    # clear plot
    plt.clf()

    linesplotted = 0
    for (i, val) in enumerate(sort(intersect(uniquesweepvalues, valuestoplot)))
        linesplotted += 1
        # find index of run in big matrix
        indexinmat = findfirst(uniquesweepvalues .== val)
        linecolor = linecolors[linesplotted]
        plt.plot(dtmidpoints, derivlog[i, :], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    plt.legend(legendlabels, title=legendtitle, title_fontsize=16, prop=Dict("family" => FONT, "size" => 8), frameon=false, handlelength=1.0, ncol=3)
    plt.xlim([cleandtdata[1], maxdt])
    plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
    plt.ylabel(raw"Power Law Exponent of MSD($t$)", fontsize=16, fontname=FONT)
    
    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")

    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        plotfilename = checkfilename(joinpath(datadirname, "$(FILE_NAME_PREFIX)_loglog_slope.pdf"))
        plt.savefig(plotfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end


"""
Plots the data previously saved to a txt file.

Can specify `maxindex` to only plot that many data points along the horizontal axis.
"""
function msd_plot_txt_logslope(txtfilename::String, maxdt::Real=-1, saveplot::Bool=false, show::Bool=true)
    # load data from file
    xydatamatrix, sweeptype, legendvalues, paramsstr = msd_loadtxt(txtfilename)

    nshots = length(legendvalues)

    # if no `maxdt` provided, defaults to maximum in data
    if maxdt == -1
        maxdt = xydatamatrix[1, end]
    end

    # clean data for log plot by removing 0s
    cleandtdata = xydatamatrix[1, xydatamatrix[1, :] .!= 0]
    cleanmsddata = xydatamatrix[2:end, xydatamatrix[1, :] .!= 0]

    # calculate slope of logplots
    derivlog = (log10.(cleanmsddata[:, 2:end]) .- log10.(cleanmsddata[:, 1:end-1])) ./ permutedims((log10.(cleandtdata[2:end]) .- log10.(cleandtdata[1:end-1])))
    dtmidpoints = (cleandtdata[1:end-1] .+ cleandtdata[2:end]) ./ 2

    ### Plotting

    # clear plot
    plt.clf()

    for i in 1:nshots
        plt.plot(dtmidpoints, derivlog[i, :])
    end

    plt.xlim([cleandtdata[1], maxdt])
    plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
    plt.ylabel(raw"Slope of the LogLog MSD($t$) (Power Law Exponent)", fontsize=16, fontname=FONT)
    
    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")

    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        plotfilename = checkfilename(joinpath(datadirname, "$(FILE_NAME_PREFIX)_loglog_slope.pdf"))
        plt.savefig(plotfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end



"""
Produces a line plot of the mean squared displacement as a function of the time interval `dt`
of all the files in a given directory, and plots the lines on the same axes.

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function msd_plot_dir(dirname::String, sweepname::String, sweeptype::SweepType, settletime::Float64=-1.0, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
        
    # initialize data matrix
    local full_msdmat::Matrix{Float64} = zeros(0,0)
    
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

        # get the position data from the file
        positions, simparams = loadsimpositions(filename)
        push!(sp_array, simparams)

        # run correlation function and add to plot
        msdmat_i::Matrix{Float64} = meansqdisp(positions, simparams, settletime)


        # if `full_msdmat` is empty, add x-axis data and first line of y-axis data
        if length(full_msdmat) == 0
            full_msdmat = zeros(length(datafiles) + 1, size(msdmat_i)[2])

            full_msdmat[1:2, :] = msdmat_i
            
        else
            # else just add y-axis data
            full_msdmat[i+1, :] = msdmat_i[2,:]
        end
    end

        # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysissweepdir(MSD_DIRNAME, sp_array, sweepname), "$(FILE_NAME_PREFIX)-settletime$(settletime).txt")

        open(outputtextfilefullpath, "w") do io
            println(io, "Mean Squared Displacement Data")
            println(io, "Simulation Parameters")
            println(io, csv_serialize(sd.simparams))
            println(io, "Plot Data")
            writedlm(io, full_msdmat, ",")
        end
    end

    
    # generate legend labels and title

    local interactionstr::String = ""
    if sp_array[1].interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(sp_array[1].interaction) I=$(sp_array[1].interactionfliprate)"
    end
    
    local legendlabels::Array{String}
    local legendtitle::String
    local paramsstr::String
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
        legendlabels = (sp -> "$(sp.numparticles)").(sp_array)
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        legendlabels = (sp -> "$(sp.boxwidth)").(sp_array)
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sd.simparams.numparticles)   Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
    elseif sweeptype == interactionstrengthsweep
        legendtitle = "Interaction Fliprate"
        legendlabels = (sp -> "$(sp.interactionfliprate)").(sp_array)
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sd.simparams.numparticles)   $(interactionstr)"
    end



    

    ### Plotting

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)

    for i in 2:size(full_msdmat)[1]
        plt.plot(full_msdmat[1,:], full_msdmat[i,:])
    end

    plt.legend(legendlabels, title=legendtitle)
    plt.xlabel("$(raw"Time Interval $dt$")")
    plt.ylabel("Mean Squared Displacement")
    plt.title("Mean Squared Displacement of Active Particle Simulation\n$(paramsstr)")
    
    if saveplot
        # get dir path to save plot
        datadirname = getanalysissweepdir(MSD_DIRNAME, sp_array, sweepname)
        plt.savefig(joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf"), bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end



# dirpath = fixpath("work/data/sweeps/alignsimple/interactionsweep/Aug15-B50-Isweep/firsthalf")
# msd_savetxt_dir(dirpath, "Aug15_rsweep_half2", interactionstrengthsweep, 100, 20)


txtfile = fixpath("work/analysis/Mean Squared Displacement/Align Simple/Aug15_rsweep/msd_sweep_plot.txt")
# txtfile2t1 = fixpath("work/analysis/Mean Squared Displacement/Align Simple/rsweep_t1/mean_spin_sweep_plot_0.txt")

msd_plot_txt(txtfile, 10, true, true, [1, 2, 5, 10, 20, 50, 100, 200, 500])
# msd_plot_txt_logslope(txtfile2t1, -1, false, true)