########################
# orientationcorrelation_driver.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run orientation self-correlation function on variou data files and plot results
########################

using PyPlot
using PyCall

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../../file_management/analysis_files.jl")


@enum SweepType densitysweep boxwidthsweep interactionstrengthsweep

MSD_DIRNAME = "Mean Squared Displacement"

FILE_NAME_PREFIX = "msd"

"""
Produces a line plot of the mean squared displacement as a function of the time interval `dt`

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function msd_plot(filename::String, settletime::Float64=-1.0, serialized::Bool=false, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
    

    
    # get the `SimulationData` from the file
    local sd::SimulationData
    if serialized
        sd = loadcompressedfile(filename)
    else
        sd = loadsim(filename, rowwisetxt)
    end

    msdmat::Matrix{Float64} = meansqdisp(sd, settletime)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(MSD_DIRNAME, sd.simparams), "$(FILE_NAME_PREFIX)-$(settletime).txt")

        open(outputtextfilefullpath, "w") do io
            println(io, "Mean Squared Displacement Data")
            println(io, "Simulation Parameters")
            println(io, csv_serialize(sd.simparams))
            println(io, "Plot Data")
            writedlm(io, msdmat, ",")
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
    plt.plot(msdmat[1,:], msdmat[2,:])
    plt.xlabel("$(raw"Time Interval $dt$")")
    plt.ylabel("Mean Squared Displacement")
    plt.title("Mean Squared Displacement of Active Particle Simulation\n N=$(sd.simparams.numparticles)  Boxwidth=$(sd.simparams.boxwidth)  t=$(Int64(round(sd.simparams.totaltime)))  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(MSD_DIRNAME, sd.simparams)
        plt.savefig(joinpath(datadirname, "$(FILE_NAME_PREFIX).pdf"), bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end


function msd_plot_txt(txtfilename::String, maxindex::Int64=1000, saveplot::Bool=true, show::Bool=false)
    
    # get the plot data from the file
    file = open(txtfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert contains(lowercase(line), "mean squared displacement data") "MSD data file does not have correct first line. File may be corrupted."
    
    line = readline(file)
    @assert contains(lowercase(line), "simulation parameters") "MSD data file does not have correct first line. File may be corrupted."
    
    # read simparams
    # read next line to get number of position data points
    simparaminfo_str = readline(file)
    simparams = SimulationParameters(simparaminfo_str)
    
    line = readline(file)
    @assert contains(lowercase(line), "plot data") "MSD data file does not have correct first line. File may be corrupted."
    
    # read x data (dt)
    dataline = readline(file)
    dtdata::Array{Float64} = parse.(Float64, split(dataline, ","))
    
    # read y data (MSD)
    dataline = readline(file)
    msddata::Array{Float64} = parse.(Float64, split(dataline, ","))



    local interactionstr::String = ""
    if simparams.interaction == nointeraction
        interactionstr = "no-interaction  "
    else
        interactionstr = "$(simparams.interaction) I=$(Int64(round(simparams.interactionfliprate)))  "
    end

    # calculate slope of logplots
    derivlog = (log10.(msddata[2:maxindex+1]) .- log10.(msddata[1:maxindex])) ./ (log10.(dtdata[2:maxindex+1]) .- log10.(dtdata[1:maxindex]))

    # more file naming controls


    ### Plotting

    # clear plot
    plt.clf()
    
    plt.grid(true, zorder=0)
    # plt.plot(dtdata[1:maxindex], msddata[1:maxindex])
    # plt.plot(log10.(dtdata[1:maxindex]), log10.(msddata[1:maxindex]))
    # plt.xlabel("$(raw"Time Interval $dt$")")
    # plt.ylabel("Mean Squared Displacement")
    plt.plot(log10.(dtdata[1:maxindex]), derivlog)
    plt.xlabel("$(raw"Log of Time Interval $\log(dt)$")")
    plt.ylabel("$(raw"Slope of the Log of the Mean Squared Displacement \n $\frac{d\log(\langle\Delta^2\rangle)}{d\log(dt)}$")")
    plt.title("Slope of Log-Log Plot of Mean Squared Displacement\nof Active Particle Simulation\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(simparams.totaltime)  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        plt.savefig(joinpath(datadirname, "$(FILE_NAME_PREFIX)_loglog_slope_all.pdf"), bbox_inches = "tight", pad_inches=0.1)
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
function msd_plot_dir(dirname::String, sweepname::String, sweeptype::SweepType, settletime::Float64=-1.0, serialized::Bool=false, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
        
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

        msdmat_i::Matrix{Float64} = meansqdisp(sd, settletime)


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

datafile_100 = fixpath("work/data/sweeps/densitysweep/collapsed/N100-B100.0-alignsimple-100.txt")
datafile_100_txt = fixpath("work/analysis/Mean Squared Displacement/Align Simple/N100-B100-T100-I100/msd--1.0.txt")
datafile_nointeraction = fixpath("work/data/22-6/N10000-nointeraction-t100-sn0.01.txt")
datafile_nointeraction_2 = fixpath("work/data/8-7/N1000-B100.0-nointeraction-100-T100.txt")
datafile_nointeraction_2_txt = fixpath("work/analysis/Mean Squared Displacement/No Interaction/N1000-B100-T100/msd--1.0.txt")



serialized::Bool = true
# msd_plot(datafile_100, -1.0, serialized, true, true, true)
# msd_plot(datafile_nointeraction_2, -1.0, serialized, true, true, true)
    
msd_plot_txt(datafile_nointeraction_2_txt, 4999, false, true)
# testdir = fixpath("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/sweeps/densitysweep/collapsed")
# oc_plot_dir(testdir, "Density Sweep 10-7", densitysweep, -1.0, serialized, true, true, false)
