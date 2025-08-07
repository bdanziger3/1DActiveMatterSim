########################
# analysis_functions.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing functions that compute scientific functions on analysis data from correlation functions.
########################

using PyPlot
using PyCall
using ColorSchemes
using FFTW
using LsqFit

include("../simulation/sim_structs.jl")


"""
Plots the FFT of a given data matrix.

Assumes `datamatrix` is an array with 2 rows. x-data in row 0 and y-data in row 1
"""
function plotfft(datamatrix::Matrix{<:Real}, simparams::SimulationParameters, outfilepath::String="", show::Bool=false)

    settleddata = datamatrix[2,10000:end]
    transformeddata = rfft(settleddata)
    # transformeddata = fft(settleddata)

    nsamples = length(settleddata)
    samplingfreq = 1 / simparams.snapshot_dt
    
    freqs = rfftfreq(nsamples, samplingfreq)
    # freqs = fftfreq(nsamples, samplingfreq)
    nfreqs = length(freqs)

    maxabs = maximum(abs.(transformeddata))

    m = findall(abs.(transformeddata) .>= 0.05 * maxabs)

    println(m[end])


    plt.clf()
    plt.grid(true, zorder=0)
    plt.plot(freqs, transformeddata ./ nsamples)
    plt.xlim([0, freqs[m[end]]])
    plt.ylabel("Normalised Amplitude")
    plt.xlabel("Frequency (Hz)")
    plt.title("FFT of Mean Spin\n N=$(simparams.numparticles)  Boxwidth=$(simparams.boxwidth)  t=$(Int64(simparams.totaltime))  $(simparams.interaction) I=$(simparams.interactionfliprate)")
    # plt.plot(freqs[1:20000], transformeddata[1:20000])
    # plt.plot(freqs, transformeddata)
    # plt.plot(datamatrix[1,:], reconstructed)
    # plt.plot(datamatrix[1,:], datamatrix[2,:])
    # plt.plot(datamatrix[1,:], reconstructed.-datamatrix[2,:])

    if outfilepath != ""
        plt.savefig(outfilepath, bbox_inches = "tight", pad_inches=0.1)
    end
    
    
    if show
        plt.show()
    end


    plt.clf()
    plt.grid(true, zorder=0)
    plt.plot(1 ./ (2 .* freqs[2:100]), transformeddata[2:100] ./ nsamples)
    plt.ylabel("Normalised Amplitude")
    plt.xlabel("Reversal Time (s)")
    plt.xlim([0, 50])
    
    if outfilepath != ""
        plt.savefig("$(outfilepath[1:end-4])_revtime.pdf", bbox_inches = "tight", pad_inches=0.1)
    end




end

function getfftreversaltime(podata::Matrix{<:Real}, snapshot_dt::Float64=0.01)::Float64
    settledindex = 10000
    settleddata = podata[:, settledindex:end]

    # Run FFT on data
    transformeddata_i = rfft(settleddata[2, :])

    nsamples = size(settleddata)[2]
    samplingfreq = 1 / snapshot_dt
    freqs = rfftfreq(nsamples, samplingfreq)

    max_val, max_i = findmax(abs.(transformeddata_i))
    maxfreq = freqs[max_i]
    
    f_avg = sum(freqs .* abs.(transformeddata_i)) / sum(abs.(transformeddata_i)) # Spectral centroid
    reversaltime_avg = sum((1 ./ (2 .* freqs[2:end])) .* abs.(transformeddata_i[2:end])) / sum(abs.(transformeddata_i[2:end])) # Spectral centroid

    println("maxfreq: $(maxfreq), avgfreq: $(f_avg), average reversal time: $(reversaltime_avg).")

    return reversaltime_avg
end


function plotfftsweep(datamatrix::Matrix{<:Real}, sweepvalues::Array{<:Real}, outputplotdir::String="", sweepcolormap=ColorSchemes.balance, snapshot_dt::Float64=0.01, show::Bool=false)

    nshots = size(datamatrix)[1] - 1
    settledindex = 10000

    settleddata = datamatrix[:, settledindex:end]

    plt.clf()
    plt.grid(true, zorder=0)

    maxfreqs = zeros(nshots)
    avgfreqs = zeros(nshots)
    avgrevtimes = zeros(nshots)

    for i in 1:nshots
        transformeddata_i = rfft(settleddata[i+1, :])

        nsamples = size(settleddata)[2]
        samplingfreq = 1 / snapshot_dt
        freqs = rfftfreq(nsamples, samplingfreq)
        nfreqs = length(freqs)

        linecolor = get(sweepcolormap, (i / nshots))
        plt.plot(freqs[1:100], transformeddata_i[1:100] ./ nsamples, color=(linecolor.r, linecolor.g, linecolor.b))

        max_val, max_i = findmax(abs.(transformeddata_i))
        maxfreqs[i] = freqs[max_i]
        
        f_avg = sum(freqs .* abs.(transformeddata_i)) / sum(abs.(transformeddata_i)) # Spectral centroid
        reversaltime_avg = sum((1 ./ (2 .* freqs[2:end])) .* abs.(transformeddata_i[2:end])) / sum(abs.(transformeddata_i[2:end])) # Spectral centroid
        avgfreqs[i] = f_avg
        avgrevtimes[i] = reversaltime_avg


    end

    legendlabels = (val -> "$val").(sweepvalues)
    plt.legend(legendlabels, title="Interaction Strength", loc="upper right")
    
    if !isempty(outputplotdir)
        filesavepath = joinpath(outputplotdir, "FFT_plot.pdf")
        plt.savefig(filesavepath, bbox_inches = "tight", pad_inches=0.1)
    end
    
    if show
        plt.show()
    end

    plt.clf()


    plt.grid(true, "both", zorder=0)
    plt.plot(sweepvalues, avgrevtimes, "o", markerfacecolor="none", markeredgecolor="seagreen", label="Data")
    plt.xlabel(raw"Interaction Fliprate")
    plt.ylabel(raw"Average Reversal Time $\langle T \rangle$")

    if !isempty(outputplotdir)
        filesavepath = joinpath(outputplotdir, "FFT_Reversal_Time_plot.pdf")
        plt.savefig(filesavepath, bbox_inches = "tight", pad_inches=0.1)
    end

    if show
        plt.show()
    end
    plt.clf()



    ## LOG LOG PLOT
    # compute best fit line of log

    # remove 0 if in data
    if sweepvalues[1] == 0
        sweepvalues = sweepvalues[2:end]
        avgrevtimes = avgrevtimes[2:end]
    end


    model(x, p) = (p[1] .* x) .+ p[2]

    # Generate noisy data
    xdata = log10.(sweepvalues)
    ydata = log10.(avgrevtimes)

    # Initial parameter guess: [a, b, c]
    p0 = [0.3, 0.0]

    # fit only data after ~I=10 since that is when we see the clusters form visually
    min_sweep = 2.5
    fit = curve_fit(model, xdata[sweepvalues .>= min_sweep], ydata[sweepvalues .>= min_sweep], p0)
    p_est = fit.param
    println("Estimated parameters: ", p_est)

    bestfit_x = collect(1:200:1001)
    bestfit_logx = log10.(bestfit_x)
    bestfit_logy = model(bestfit_logx, p_est)
    bestfit_y = exp10.(bestfit_logy)



    plt.grid(true, "both", zorder=0)
    # plt.plot(sweepvalues, avgrevtimes, "-o")
    plt.plot(sweepvalues, avgrevtimes, "o", markerfacecolor="none", markeredgecolor="seagreen", label="Data")
    plt.plot(bestfit_x, bestfit_y, "--", label="Best Fit", color="black")
    
    plt.text("m = $(p[1])")

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel(raw"Interaction Fliprate ")
    plt.ylabel(raw"Average Reversal Time $\langle T \rangle$")

    
    if !isempty(outputplotdir)
        filesavepath = joinpath(outputplotdir, "FFT_Reversal_Time_loglog_plot.pdf")
        plt.savefig(filesavepath, bbox_inches = "tight", pad_inches=0.1)
    end

    if show
        plt.show()
    end



end