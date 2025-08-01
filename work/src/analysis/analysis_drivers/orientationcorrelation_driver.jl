########################
# orientationcorrelation_driver.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run orientation self-correlation function on variou data files and plot results
########################

using PyPlot
using PyCall
using ColorSchemes

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../../file_management/analysis_files.jl")


@enum SweepType densitysweep boxwidthsweep interactionstrengthsweep


OC_LINE_COLORMAP = ColorSchemes.vik100



"""
Produces a line plot of the orientation self-correlation as a function of the time interval `dt`

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function oc_plot(filename::String, settletime::Float64=-1.0, serialized::Bool=false, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
    
    ORIENTATION_SELFCORRELATION_DIRNAME = "Orienation Self-Correlation"

    FILE_NAME_PREFIX = "orientation_selfcorrelation"

    
    # get the `SimulationData` from the file
    local sd::SimulationData
    if serialized
        sd = loadcompressedfile(filename)
    else
        sd = loadsim(filename, rowwisetxt)
    end

    ocmat::Matrix{Float64} = orientationcorrelation(sd, settletime)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(ORIENTATION_SELFCORRELATION_DIRNAME, sd.simparams), "$(FILE_NAME_PREFIX)-$(settletime).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)
        open(outputtextfilefullpath, "w") do io
            println(io, "Orientation Self-Correlation  Data")
            writedlm(io, ocmat, ",")
        end
    end



    local interactionstr::String = ""
    if sd.simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(sd.simparams.interaction) I=$(Int64(round(sd.simparams.interactionfliprate)))  "
    end


    # more file naming controls


    ### Plotting

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)
    plt.plot(ocmat[1,:], ocmat[2,:])
    plt.xlabel("$(raw"Time Interval $dt$")")
    plt.ylabel("Orientation Self-Correlation\n$(raw"$\langle u_i (t + \Delta t) u_t(t)\rangle_{t,i}$")")
    plt.title("Orientation Self-Correlation of Active Particle Simulation\n N=$(sd.simparams.numparticles)  Boxwidth=$(sd.simparams.boxwidth)  t=$(Int64(sd.simparams.totaltime))  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(ORIENTATION_SELFCORRELATION_DIRNAME, sd.simparams)
        outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end


"""
Plots the data previously saved to a txt file.

Can specify `maxindex` to only plot that many data points along the horizontal axis.
"""
function oc_plot_txt(txtfilename::String, maxindex::Int64=1000, saveplot::Bool=false, show::Bool=true)
    
    # get the plot data from the file
    file = open(txtfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert contains(lowercase(line), "orientation self-correlation data") "OC data file does not have correct first line. File may be corrupted."
    
    # check if it is a sweep
    line = readline(file)

    ### Plotting
    # prepare plot
    plt.clf()
    plt.grid(true, zorder=0)
    ###

    
    if contains(lowercase(line), "sweep type")
        FILE_NAME_PREFIX = "orientation_selfcorrelation_plot"
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
        println(legendvalues)
        
        # get params string for plot
        line = readline(file)
        paramsstr = line[length("Params String: ")+1:end]
        
        # read x data (dt)
        line = readline(file)
        dataline = readline(file)
        dtdata::Array{Float64} = parse.(Float64, split(dataline, ","))
        # read y data (OC)
        ocdata_i::Array{Float64} = zeros(1, length(dtdata))
        nshots = length(legendlabels)

        for i in 1:nshots
            dataline = readline(file)
            ocdata_i = parse.(Float64, split(dataline, ","))
            
            # plot each line one at a time
            positivehalf = true     # using positive half of the color scheme?
            if positivehalf
                linecolor = get(OC_LINE_COLORMAP, 0.5 + (0.5 * (i / nshots)))
            else
                linecolor = get(OC_LINE_COLORMAP, 0.5 - (0.5 * (i / nshots)))
            end
            # linecolor = get(OC_LINE_COLORMAP, legendvalues[i] / maximum(legendvalues))
            println(linecolor)
            println((linecolor.r, linecolor.g, linecolor.b))
            plt.plot(dtdata[1:maxindex], ocdata_i[1:maxindex], color=(linecolor.r, linecolor.g, linecolor.b))
        end
        
        # add legends and labels to plot
        plt.legend(legendlabels, title=legendtitle)
        plt.xlabel("$(raw"Time Interval $dt$")")
        plt.ylabel("Orientation Self-Correlation\n$(raw"$\langle u_i (t + \Delta t) u_t(t)\rangle_{t,i}$")")
        plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(paramsstr)")
        
        if saveplot
            # get dir path to save plot
            datadirname = dirname(txtfilename)
            outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX)_from_txt_maxindex-$(maxindex).pdf")
            outputfilename = checkfilename(outputfilename)
            plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
        end
        
        if show
            plt.show()
        end

    else
        FILE_NAME_PREFIX = "orientation_selfcorrelation"
        @assert false "Cannot plot files that are not sweeps yet."
    end
end


"""
Produces a line plot of the orientation self-correlation as a function of the time interval `dt`
of all the files in a given directory, and plots the lines on the same axes.

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function oc_plot_dir(dirname::String, sweepname::String, sweeptype::SweepType, settletime::Float64=-1.0, serialized::Bool=false, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
    
    ORIENTATION_SELFCORRELATION_DIRNAME = "Orienation Self-Correlation"

    FILE_NAME_PREFIX = "orientation_selfcorrelation_plot"
    
    # initialize data matrix
    local full_ocmat::Matrix{Float64} = zeros(0,0)
    
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
        local sd::SimulationData
        if serialized
            sd = loadcompressedfile(inputfilename)
        else
            sd = loadsim(inputfilename, rowwisetxt)
        end

        ocmat_i::Matrix{Float64} = orientationcorrelation(sd, settletime)


        # if `full_ocmat` is empty, add x-axis data and first line of y-axis data
        if length(full_ocmat) == 0
            full_ocmat = zeros(length(datafiles) + 1, size(ocmat_i)[2])

            full_ocmat[1:2, :] = ocmat_i
            
        else
            # else just add y-axis data
            full_ocmat[i+1, :] = ocmat_i[2,:]
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
    local sortedorder::Vector{Int}
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
        sortedorder = sortperm([sp.numparticles for sp in sp_array])
        legendlabels = (sp -> "$(sp.numparticles)").(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  Boxwidth=$(sp_array[1].boxwidth)  $(interactionstr)"
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
        sortedorder = sortperm([sp.boxwidth for sp in sp_array])
        legendlabels = (sp -> "$(sp.boxwidth)").(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  $(interactionstr)"
    elseif sweeptype == interactionstrengthsweep
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        legendlabels = (sp -> "$(sp.interactionfliprate)").(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth)"
    end


    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysissweepdir(ORIENTATION_SELFCORRELATION_DIRNAME, sp_array, sweepname), "$(FILE_NAME_PREFIX)-settletime$(settletime).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)

        open(outputtextfilefullpath, "w") do io
            println(io, "Orientation Self-Correlation Data")
            println(io, "Sweep Type: $(sweeptype)")
            println(io, "Sweep Values")
            writedlm(io, permutedims(legendlabels), ",")
            println(io, "Params String: $(paramsstr)")
            println(io, "Data Matrix")
            writedlm(io, permutedims(full_ocmat[1, :]), ",")
            for i in sortedorder
                writedlm(io, permutedims(full_ocmat[i+1, :]), ",")
            end
        end
    end





    ### Plotting

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)

    for i in sortedorder
        linecolor = get(OC_LINE_COLORMAP, sp_array[i].numparticles / sp_array[sortedorder[end]].numparticles)
        plt.plot(full_ocmat[1,:], full_ocmat[i+1,:], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    plt.legend(legendlabels, title=legendtitle)
    plt.xlabel("$(raw"Time Interval $dt$")")
    plt.ylabel("Orientation Self-Correlation\n$(raw"$\langle u_i (t + \Delta t) u_t(t)\rangle_{t,i}$")")
    plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(paramsstr)")
    
    if saveplot
        # get dir path to save plot
        datadirname = getanalysissweepdir(ORIENTATION_SELFCORRELATION_DIRNAME, sp_array, sweepname)
        outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end

datafile = fixpath("work/data/sweeps/densitysweep/N500-B100.0-alignsimple-100.txt")
datafile2 = fixpath("work/data/26-6/N5000-B100-alignsimple-100-t4-sn0.01.txt")
datafile_nointeraction = fixpath("work/data/22-6/N10000-nointeraction-t100-sn0.01.txt")
datafile_nointeraction_2 = fixpath("work/data/8-7/N1000-B100.0-nointeraction-100-T100.txt")



serialized::Bool = true
# getposdensitydata(datafile, false, datafile_out)


# oc_plot(datafile, -1.0, serialized, false, true, true)
    

# sweepdir = fixpath("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/sweeps/densitysweep/collapsed")
# sweepdir = fixpath("work/data/sweeps/alignsimple/interactionsweep/N1000-sweep")
# oc_plot_dir(sweepdir, "Interaction Sweep N1000", interactionstrengthsweep, -1.0, serialized, false, false, true)

dsweep_txt = fixpath("work/analysis/Orienation Self-Correlation/Align Simple/Interaction Sweep N1000/orientation_selfcorrelation_plot-settletime-1.0_0.txt")
dsweep_txt2 = fixpath("work/analysis/Orienation Self-Correlation/Align Simple/Density Sweep 10-7/orientation_selfcorrelation_plot-settletime-1.0.txt")
# oc_plot_txt(dsweep_txt2, 1000, false, true)
oc_plot_txt(dsweep_txt2, 1000, true, false)


