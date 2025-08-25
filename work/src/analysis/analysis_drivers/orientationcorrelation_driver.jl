########################
# orientationcorrelation_driver.jl
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to run orientation self-correlation function on variou data files and plot results
########################

using PyPlot
using PyCall
using ColorSchemes
using StatsBase

include("../correlation_functions.jl")
include("../../file_management/simdata_files.jl")
include("../../utils/paths.jl")
include("../analysis_functions.jl")
include("../../file_management/analysis_files.jl")


# OC_LINE_COLORMAP = ColorSchemes.vik100
# OC_LINE_COLORMAP = ColorSchemes.Blues_8
OC_LINE_COLORMAP = ColorSchemes.bam50



FONT = "Times New Roman"


# Use serif font in math text
rc("mathtext", fontset="cm")
PyPlot.rc("font", family="serif")           
PyPlot.rc("font", family=FONT)



"""
Produces a line plot of the orientation self-correlation as a function of the time interval `dt`

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function oc_plot(filename::String, settletime::Real=-1, maxdt::Real=-1, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false)
    
    ORIENTATION_SELFCORRELATION_DIRNAME = "Orienation Self-Correlation"

    FILE_NAME_PREFIX = "orientation_selfcorrelation"

    spins, simparams = loadsimspins(filename)


    ocmat::Matrix{Float64} = orientationcorrelation(spins, simparams, settletime, maxdt)

    # save as a datafile if requested
    if savetxt
        outputtextfilefullpath = joinpath(getanalysisdir(ORIENTATION_SELFCORRELATION_DIRNAME, simparams), "$(FILE_NAME_PREFIX)-$(settletime).txt")
        outputtextfilefullpath = checkfilename(outputtextfilefullpath)
        open(outputtextfilefullpath, "w") do io
            println(io, "Orientation Self-Correlation  Data")
            writedlm(io, ocmat, ",")
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
    
    plt.grid(true, zorder=0)
    plt.plot(ocmat[1,:], ocmat[2,:])
    plt.xlabel("$(raw"Time Interval $dt$")")
    plt.ylabel("Orientation Self-Correlation\n$(raw"$\langle u_i (t + \Delta t) u_t(t)\rangle_{t,i}$")")
    plt.title("Orientation Self-Correlation of Active Particle Simulation\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(interactionstr)")
    

    if saveplot
        # get dir path to save plot
        datadirname = getanalysisdir(ORIENTATION_SELFCORRELATION_DIRNAME, simparams)
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

Can specify `maxdt` to only plot up to that many data points along the horizontal axis.
"""
function oc_plot_txt(txtfilename::String, maxdt::Real=-1, saveplot::Bool=false, show::Bool=true)
    
    # get the plot data from the file
    xydatamat, sweeptype, legendvalues, simparamsout = oc_load_txt(txtfilename)

    # handle `maxdt`
    # local maxindex::Int64 = 0
    # if maxdt == -1
    #     maxindex = size(xydatamat)[2]
    # else
    #     maxindex = min(Int64(round(maxdt / (xydatamat[1,2]))) + 1, size(xydatamat)[2])
    # end

    nshots::Int64 = size(xydatamat)[1] - 1

    # make legend title
    local legendtitle
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
    elseif sweeptype == interactionstrengthsweep
        legendtitle = "Interaction Fliprate"
    end

    ### Plotting
    # prepare plot
    plt.clf()
    plt.grid(true, zorder=0)
    ###

    for i in 1:nshots
        # plot each line one at a time
        positivehalf = true     # using positive half of the color scheme?
        if positivehalf
            linecolor = get(OC_LINE_COLORMAP, 0.5 + (0.5 * (i / nshots)))
        else
            linecolor = get(OC_LINE_COLORMAP, 0.5 - (0.5 * (i / nshots)))
        end

        plt.plot(xydatamat[1, :], xydatamat[1+i, :], color=(linecolor.r, linecolor.g, linecolor.b))
    end

    if maxdt > 0
        plt.xlim([0, maxdt])
    end

    # add legends and labels to plot
    plt.legend(legendvalues, title=legendtitle)
    plt.xlabel("$(raw"Time Interval $t$")")
    plt.ylabel("Orientation Self-Correlation\n$(raw"$C_p(t)$")")
    plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(simparamsout)")
    
    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        outputfilename = joinpath(datadirname, "$(FILE_NAME_PREFIX)_from_txt_maxdt-$(maxdt).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end
end

"""
Loads the data previously saved to a txt file.

Returns a tuple of 4 objects:


    `xydatamatrix`:   The matrix holding the analysis data. Row 0 holds the x-axis data.
                        The remaining lines hold the y-axis data
                        
    `sweeptype`:        The `SweepType` of the data. Is `nosweep` if not a sweep.

    `legendvalues`:     The list of values of the swept parameter in the sweep. If the file is not for a sweep, contains `[nothing]`.

    `simparamsout`:     The simulation parameters used. Returns a `SimulationParameters` object if not a sweep.
                        Returns a string to use in the plot title if it is a sweep
"""
function oc_load_txt(txtfilename::String)
    
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

    local simparamsout
    local nshots
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
        interactionstr = "$(sp_array[1].interaction)"
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        legendlabels = (sp -> "$(sp.interactionfliprate)").(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth). $(interactionstr)"
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

"""
Averages the OC data over multiple runs that are saved to a file
"""
function oc_average_reps_txt(txtfilename::String, maxdt::Real=-1, saveplot::Bool=false, show::Bool=true, valuestoplot=nothing)
    oc_plot_name_prefix = "oc_plot"
    persistence_time_plot_name_prefix = "ptime_plot"

    
    # get the plot data from the file
    xydatamat, sweeptype, legendvalues, simparamsout = oc_load_txt(txtfilename)
    
    # for fitting to decay
    # fit_length = 8
    fit_length = maxdt
    fit_samples = Int64(round(fit_length / xydatamat[1,2]))

    # find the unique sweep values
    sweepvaluecounts = countmap(legendvalues)
    uniquesweepvalues = sort(collect(keys(sweepvaluecounts)))
    nvalues = length(uniquesweepvalues)


    # average data with same sweep value
    averagedxydatamat::Matrix{Float64} = zeros(nvalues+1, size(xydatamat)[2])
    averagedxydatamat[1, :] = xydatamat[1, :]

    indivtaus = [[] for _=1:nvalues]


    for (i, val) in enumerate(legendvalues)
        legendvalue_uniqueindex = findfirst(uniquesweepvalues .== val)
        averagedxydatamat[legendvalue_uniqueindex+1, :] += (xydatamat[i+1, :] ./ sweepvaluecounts[val])


        tau_i, _ = fitdecay(xydatamat[1, 1:fit_samples], xydatamat[i+1, 1:fit_samples])
        push!(indivtaus[legendvalue_uniqueindex], tau_i) 
    end


    # handle `maxdt`
    # local maxindex::Int64 = 0
    # if maxdt == -1
    #     maxindex = size(xydatamat)[2]
    # else
    #     maxindex = min(Int64(round(maxdt / (xydatamat[1,2]))) + 1, size(xydatamat)[2])
    # end

    # make legend title
    local legendtitle
    if sweeptype == densitysweep
        legendtitle = "Number of Particles"
    elseif sweeptype == boxwidthsweep
        legendtitle = "Box Width"
    elseif sweeptype == interactionstrengthsweep
        legendtitle = "Interaction Fliprate"
    end

    # handle`valuestoplot`
    if isnothing(valuestoplot)
        valuestoplot = uniquesweepvalues
    end

    ### Plotting
    # prepare plot
    plt.clf()
    ###

    linesplotted = 0
    println(uniquesweepvalues)
    for (i, val) in enumerate(sort(intersect(uniquesweepvalues, valuestoplot)))
        # find index of run in big matrix
        indexinmat = findfirst(uniquesweepvalues .== val)
        
        # plot each line one at a time
        positivehalf = true     # using positive half of the color scheme?
        # if uniquesweepvalues[i] in valuestoplot
        linesplotted += 1
        if positivehalf
            # linecolor = get(OC_LINE_COLORMAP, 0.5 + (0.5 * (val / maximum(intersect(uniquesweepvalues, valuestoplot)))))
            linecolor = get(OC_LINE_COLORMAP, 0.5 + (0.5 * (i / length(intersect(uniquesweepvalues, valuestoplot)))))
        else
            # linecolor = get(OC_LINE_COLORMAP, 0.5 - (0.5 * (i / nvalues)))
            linecolor = get(OC_LINE_COLORMAP, 0.5 - (0.5 * (i / length(intersect(uniquesweepvalues, valuestoplot)))))
        end

        # t0x = collect(0:10:100*(20-fit_length))
        # tauoft0 = zeros(length(t0x))
        # for (j, t0) in enumerate(t0x)
        #     tau, tau_stderr = fitdecay(averagedxydatamat[1, t0+1:t0+1+fit_samples], averagedxydatamat[i+1, t0+1:t0+1+fit_samples])
        #     tauoft0[j] = tau
        # end
        tau, tau_stderr = fitdecay(averagedxydatamat[1, 1:fit_samples], averagedxydatamat[indexinmat+1, 1:fit_samples])
        # println(val)
        # plt.plot(t0x./100, tauoft0, "o-")
        # plt.show()
        # plt.plot(averagedxydatamat[1, :], averagedxydatamat[1+indexinmat, :], color=(linecolor.r, linecolor.g, linecolor.b))
        plt.plot(averagedxydatamat[1, :], averagedxydatamat[1+indexinmat, :], color="dodgerblue", linewidth=1.5)

        # if val == 100 || val == 250 || val == 500
        #     xtemp = collect(0:.01:10)
        #     ytemp = exp.(- xtemp ./ tau)
        #     plt.plot(xtemp, ytemp, c="blue", alpha=val/500)
        # end

        # end
    end

    if maxdt > 0
        plt.xlim([0, maxdt])
    end
    # plt.ylim([.1, 1])

        
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

    xpts = collect(0:.01:maxdt)
    ypts = exp.(-2 .* xpts)
    plt.plot(xpts, ypts, "--", linewidth=1, c="k")


    # plt.legend(legendlabels, title=legendtitle, title_fontsize=16, prop=Dict("family" => FONT, "size" => 8), frameon=false, handlelength=1.0, ncol=3)
    plt.xlabel("$(raw"$t$")", fontsize=16, fontname=FONT)
    plt.ylabel(L"$C_p(t)$", fontsize=16, fontname=FONT)

    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")
    # plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(simparamsout)")
    
    ax = plt.gca()
    for label in cat(ax.get_xticklabels(), ax.get_yticklabels(), dims=1)
        label.set_fontname(FONT)  # font family
    end

    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        outputfilename = joinpath(datadirname, "$(oc_plot_name_prefix)_from_txt_maxdt-$(maxdt).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end


    # Now find the fits of the averaged data
    avgtaus = zeros(nvalues)
    avgtau_stderrs = zeros(nvalues)
    for i in 1:nvalues
        tau, tau_stderr = fitdecay(averagedxydatamat[1, 1:fit_samples], averagedxydatamat[i+1, 1:fit_samples])
        avgtaus[i] = tau
        avgtau_stderrs[i] = tau_stderr
    end

    avgtau2 = (l -> mean(l)).(indivtaus)
    stderrvtau = (l -> std(l) / sqrt(length(l))).(indivtaus)

    interval95p = 1.96
    # stderrvtau = zeros(length(avgtaus))
    # tau5 = [5.976437180610176,5.346917249487149,8.403178089540821,7.6934739567549615,9.793882120555722,6.9753952800077625,7.817340007708796,6.585548150809208,8.59011475485227,4.476300102908587,4.168320392984825,4.122124432879653,4.263805772019667,3.9625558733727653,4.015041357650596,3.955493678615965,3.6141829308799878,3.7715173366825168,3.1567060029064145,3.3160528608834503,3.143363787242488,3.7271000820603146,3.2282443768922398,3.159912556884323,2.650382050562255,2.5535647095944043,2.6335971461453034,2.56991767831472,2.6645820315323165,2.4858820852292705,2.264381922776225,2.293948082592277,2.519060605707039,2.481232941965637,2.3619390408997827,2.2792973052058874,2.2758596896855305,2.0829787010982996,2.4423558973221677,2.1900629154200586,2.4767847722854324,2.200889498014523,2.176305855401227,2.0365305410379606,2.4085502528089004,1.9613300353349679,1.9610572753197832,2.063741185888926,2.2171062447511,1.9225624411819482,1.9192221944782468,1.937498437810841,1.930214293941244,1.7239238759228237,1.7074058249753465,1.7402439791747721,1.6283183033447812,1.7530772464904352,1.6707025592296934,1.64543765791232,1.5669060921399762,1.6751443772158947,1.6563402467680923,1.6387304212604181,1.5754006549965545,1.596621972874794,1.585160212821649,1.590468903768549,1.5938973067760411,1.7614430756183603,1.5939898795112086,1.5537064380458976,1.4899490720576614,1.5086654382720697,1.696620057857074,1.5158724351710158,1.554667789082055,1.494354881089016,1.5232271203809562,1.5982097021725952,1.4092822875143634,1.4023234628779528,1.4934205941518508,1.5792898196577465,1.423508973596989,1.4756286677296726,1.435830013307898,1.4510909011583037,1.4468112311366403,1.4190548511102452,1.4258074332048556,1.4238665130292647,1.3870115302770236,1.4836994717323078,1.4512658933042129,1.3413270575665952,1.3760711208491971,1.3854005080822445,1.3377209512938233,1.3280966915162955,1.3163583069865572,1.2863795382955738,1.4179027318258852,1.3704554324968374,1.374392484173922,1.333044444376327,1.3459448978174497,1.4177711756221352,1.3572057362632268,1.3412830708344499]

    # for i in 1:nvalues
    #     avgtau2[i] = mean(tau5[(i-1)*5+1:i*5])
    #     stderrvtau[i] = std(tau5[(i-1)*5+1:i*5]) / sqrt(5)
    # end


    # find the linear fits of the persistence times in loglog
    criticalval1 = 5
    criticalval2 = 5
    criticalindex1 = findmin(abs.(alllegendvalues .- criticalval1))[2]
    criticalindex2 = findmin(abs.(alllegendvalues .- criticalval2))[2]
    # params, stderrs, linearmodel = logloglinearfit(alllegendvalues[1:end-5], avgtaus[1:end-5])
    
    params, stderrs, linearmodel = logloglinearfit(alllegendvalues[1:criticalindex2], avgtaus[1:criticalindex2])
    params2, stderrs2, linearmodel2 = logloglinearfit(alllegendvalues[criticalindex1:end], avgtaus[criticalindex1:end])

    println(params)
    println(stderrs)

    println(params2)
    println(stderrs2)

    # println(avgtau_stderrs)
    # println(stderrvtau)

    # plot line
    linex = collect(minimum(alllegendvalues):1:maximum(alllegendvalues))
    liney = linearmodel(linex)
    liney2 = 1 ./ (2 .+ (linex ./ 2))
    
    # plot fit lines
    plt.clf()
    # plt.plot(linex, liney, "--", c="k", linewidth=0.5)
    # plt.plot(linex, linearmodel2(linex), "--", c="k", linewidth=0.5)
    # plt.plot(linex, liney2, "--", c="r", linewidth=0.5)
    # plt.plot(linex, linearmodel2(linex), "--", c="k", linewidth=0.5)
    
    # plot potential theoretical fit
    lowlinex = collect(.5:.001:2)
    noflockyline = 1 ./ (2 .+ (lowlinex ./ 2))
    plt.plot(lowlinex, noflockyline, "--", c="r", linewidth=0.5)


    # plt.ylim([.3, 8])


    # determine marker
    local marker
    local markercolor
    if sweeptype == densitysweep
        marker = "^"
        markercolor = "firebrick"
    elseif sweeptype == interactionstrengthsweep
        marker = "^"
        markercolor = "forestgreen"
    end
    # plot data on top
    plt.plot(alllegendvalues, avgtaus, marker, markeredgecolor=markercolor,  markerfacecolor="none", markersize=5)
    plt.errorbar(alllegendvalues, avgtau2, yerr=interval95p.*stderrvtau, fmt="none", elinewidth=0.25, ecolor="k", barsabove=false, capsize=1.5)
    
    plt.xlabel(legendtitle, fontsize=24, fontname=FONT)
    plt.ylabel(L"$\tau$", fontsize=28, fontname=FONT)

    plt.tick_params(axis="x", top=true, bottom=true, direction="in", which = "both")
    plt.tick_params(axis="y", left=true, right=true, direction="in", which = "both")

    # plt.ylim([.9,11])
    axes = plt.gca()

    plt.xscale("log")
    plt.yscale("log")
    # axes.log()
    
    if saveplot
        # get dir path to save plot
        datadirname = dirname(txtfilename)
        outputfilename = joinpath(datadirname, "$(persistence_time_plot_name_prefix)_from_txt_fitlength-$(fit_length).pdf")
        outputfilename = checkfilename(outputfilename)
        plt.savefig(outputfilename, bbox_inches = "tight", pad_inches=0.1)
    end


    if show
        plt.show()
    end


    plt.plot(alllegendvalues, avgtaus)
    plt.show()

end


"""
Produces a line plot of the orientation self-correlation as a function of the time interval `dt`
of all the files in a given directory, and plots the lines on the same axes.

Set `savetxt` to save a .txt file with the orientation self-correlation function results
Set `show` to display the plot as well as saving it.
"""
function oc_plot_dir_persistencetimes(dirname::String, sweepname::String, sweeptype::SweepType, settletime::Real=-1, maxdt::Real=-1, saveplot::Bool=true, savetxt::Bool=true, show::Bool=false, title::Bool=true)
    
    ORIENTATION_SELFCORRELATION_DIRNAME = "Orienation Self-Correlation"

    FILE_NAME_PREFIX = "orientation_selfcorrelation_plot"
    persistence_time_prefix = "persistence_time_plot"
    
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

        # load spins and Simulation Params
        spins, simparams = loadsimspins(inputfilename)
        # add simparams to array
        push!(sp_array, simparams)
        
        # run correlation function and add to plot
        ocmat_i::Matrix{Float64} = orientationcorrelation(spins, simparams, settletime, maxdt)


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
        interactionstr = "$(sp_array[1].interaction)"
        legendtitle = "Interaction Fliprate"
        sortedorder = sortperm([sp.interactionfliprate for sp in sp_array])
        legendlabels = (sp -> "$(sp.interactionfliprate)").(sp_array[sortedorder])
        paramsstr = "t=$(Int64(sp_array[1].totaltime))  N=$(sp_array[1].numparticles)  Boxwidth=$(sp_array[1].boxwidth). $(interactionstr)"
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
    plt.xlabel("$(raw"Time Interval $t$")")
    plt.ylabel("Orientation Self-Correlation\n$(raw"$C_p(t)$")")

    if title
        plt.title("Orientation Self-Correlation of Active Particle Simulation\n$(paramsstr)")
    end
    
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


    # now fit the data to an exponential to find the persistence times
    if savetxt
        fit_length = 10
        fit_samples = Int64(round(fit_length / sp_array[1].snapshot_dt))
        xdata = full_ocmat[1, :]

        paramslist = []
        stderrlist = []
        for i in sortedorder
            ydata = full_ocmat[1+i, :]
            tau, tau_stderr = fitdecay(xdata[1:fit_samples], ydata[1:fit_samples])

            push!(paramslist, tau)
            push!(stderrlist, tau_stderr)
        end


        tau_txtfilepath = joinpath(getanalysissweepdir(ORIENTATION_SELFCORRELATION_DIRNAME, sp_array, sweepname), "$(persistence_time_prefix)-settletime$(settletime).txt")
        tau_realpath = checkfilename(tau_txtfilepath)

        open(tau_realpath, "w") do io
            println(io, "Fits of Orientation Self-Correlation Data")
            println(io, "Sweep Type: $(sweeptype)")
            println(io, "Params String: $(paramsstr)")
            println(io, "Data Matrix")
            writedlm(io, permutedims(legendlabels), ",")    # xdata: sweep values
            writedlm(io, permutedims(paramslist), ",")    # ydata: Tau values
            writedlm(io, permutedims(stderrlist), ",")    # stderr: Tau Standard Error values
        end
            
    end
end





# sweepdir = fixpath("work/data/sweeps/alignsimple/densitysweep/Aug13-density-sweep-B50")
# sweepdir_temp = fixpath("work/data/sweeps/alignsimple/densitysweep/tempstore")
# datatxt = fixpath("work/analysis/Orienation Self-Correlation/Align Simple/Aug15_InteractionSweep/Aug15_InteractionSweep-full/orientation_selfcorrelation_plot-settletime100.0.txt")
# datatxt = fixpath("work/analysis/Orienation Self-Correlation/Align Simple/Aug13_DensitySweep/total/orientation_selfcorrelation_plot-settletime100.0.txt")
datatxtnoint = fixpath("work/analysis/Orienation Self-Correlation/No Interaction/test1/orientation_selfcorrelation_plot-settletime100.0.txt")



# sweepdir_int = fixpath("work/data/sweeps/alignsimple/interactionsweep/Aug15-B50-Isweep/interactionsweep")
# sweepdir_int_t1 = fixpath("work/data/sweeps/alignsimple/interactionsweep/Aug15-B50-Isweep/t1")
# oc_plot_dir_persistencetimes(sweepdir_int, "Aug15_InteractionSweep-middle2", interactionstrengthsweep, 100.0, 20, true, true, false, false)
oc_average_reps_txt(datatxtnoint, 8, true, true)
