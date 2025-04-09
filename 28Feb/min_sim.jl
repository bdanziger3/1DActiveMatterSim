import Base.==
using Random
using PyPlot
using PyCall
using DelimitedFiles
using DataStructures
using DataFrames


# Minumum approach to the same code I had before

"""
Simulation parameters on which to run particles
"""
struct SimulationParameters
    numparticles::Int
    totaltime::Real
    dt::Real
    v0::Real
    temp::Real
    boxwidth::Real
    starttime::Real
    # startinglocations::Array{Real}
    # randomstarts::Bool
    # randomspins::Bool

    function # Inner constructor with some default values
        SimulationParameters(numparticles, totaltime, dt, v0, temp, boxwidth=1, starttime=0)
        return new(numparticles, totaltime, dt, v0, temp, boxwidth, starttime)
    end
    function # Inner constructor for null object with all default values
        SimulationParameters(paramsdict::Dict)
        return new(paramsdict["numparticles"], paramsdict["totaltime"], paramsdict["dt"], paramsdict["v0"], paramsdict["temp"], paramsdict["boxwidth"], paramsdict["starttime"])
    end
    function # Inner constructor for dictionary input
        SimulationParameters()
        return new(0, 0, 0, 0, 0, 0, 0)
    end
end

function asdict(simparams::SimulationParameters)
    paramsdict::Dict = Dict("numparticles" => simparams.numparticles, "totaltime" => simparams.totaltime, "dt" => simparams.dt, "v0" => simparams.v0, "temp" => simparams.temp, "boxwidth" => simparams.boxwidth, "starttime" => simparams.starttime)
    return paramsdict
end

function csv_serialize(simparams::SimulationParameters)
    return "$(simparams.numparticles),$(simparams.totaltime),$(simparams.dt),$(simparams.v0),$(simparams.temp),$(simparams.boxwidth),$(simparams.starttime)"
end

function json_serialize(simparams::SimulationParameters)::String
    # paramsdict::Dict = Dict("numparticles" => simparams.numparticles, "totaltime" => simparams.totaltime, "dt" => simparams.dt, "v0" => simparams.v0, "temp" => simparams.temp, "boxwidth" => simparams.boxwidth, "starttime" => simparams.starttime)
    # json_str = JSON.json(asdict(simparams))
    # return json_str
    return JSON.json(asdict(simparams), 2)
end

function json_deserialize_simparams(json_str::String)::SimulationParameters
    paramsdict = JSON.parse(json_str)
    simparams = SimulationParameters(paramsdict)
    return simparams
end

function gettimes(simparams::SimulationParameters)::Array{Real}
    endtime = simparams.starttime + simparams.totaltime
    return collect(simparams.starttime:simparams.dt:endtime)
end

function getntimes(simparams::SimulationParameters)::Int64
    return Int64(floor(simparams.totaltime / simparams.dt) + 1)
end

function asarray(simparams::SimulationParameters)::Array{Any}
    ar = [simparams.numparticles, simparams.totaltime, simparams.dt, simparams.v0, simparams.temp, simparams.boxwidth, simparams.starttime]
    return ar
end

function assertparams(simparams::SimulationParameters)::Bool
    @assert getntimes(simparams) == size(gettimes(simparams))[1] "error in construction of times in `gettimes()`."
    @assert dt <= totaltime "timestep `dt` must be less than `totaltime`."
    
    return true
end


struct SimulationData
    simparams::SimulationParameters
    times::Array{Real}
    positions::Matrix{Real}
    wrappedpositions::Matrix{Real}
    spins::Matrix{Int8}

end

function assertdim(simdata::SimulationData)
    # assert simparams is consistent
    assertparams(simdata.simparams)
    # assert times are correct
    ntimes::Int64 = getntimes(simdata.simparams)
    nparticles::Int64 = simparams.numparticles
    possize = size(simdata.positions)
    wrappedpossize = size(simdata.wrappedpositions)
    spinssize = size(simdata.spins)
    # assert all data arrays have the correct number of times
    @assert ntimes == possize[1] "Positions data has incorrect number of times. Should have $(ntimes), found $(possize[1])"
    @assert ntimes == wrappedpossize[1] "Wrapped positions data has incorrect number of times. Should have $(ntimes), found $(wrappedpossize[1])"
    @assert ntimes == spinssize[1] "Spins data has incorrect number of times. Should have $(ntimes), found $(spinssize[1])"
    
    # assert all data arrays have the correct number of particles
    @assert nparticles == possize[2] "Positions data has incorrect number of particles. Should have $(nparticles), found $(possize[2])"
    @assert nparticles == wrappedpossize[2] "Wrapped positions data has incorrect number of particles. Should have $(nparticles), found $(wrappedpossize[2])"
    @assert nparticles == spinssize[2] "Spins data has incorrect number of particles. Should have $(nparticles), found $(spinssize[2])"
end

function json_serialize(simdata::SimulationData)
    datadict::Dict = Dict(Pair{String, Dict}("simparams", asdict(simdata.simparams)), Pair{String, Array{Real}}("times", simdata.times), Pair{String, Matrix{Real}}("positions", simdata.positions), Pair{String, Matrix{Real}}("wrappedpositions", simdata.wrappedpositions), Pair{String, Matrix{Int8}}("spins", simdata.spins))
    # println(datadict)
    json_str = JSON.json(datadict, 2)
    # println(json_str)
    # return json_str
end

function Base.:(==)(simdata1::SimulationData, simdata2::SimulationData)::Bool
    if simdata1.simparams != simdata2.simparams
        return false
    end

    ## TODO finish implementing checking data

    return true
end

# function json_deserialize_simdata(json_str::String)
#     simdata_1 = JSON.parse(json_str)
#     println(simdata_1)
#     simparams = SimulationParameters(simdata_1["simparams"])
#     # println(simparams)
#     simdata = SimulationData(simparams, simdata_1["times"], simdata_1["positions"], simdata_1["wrappedpositions"], simdata_1["spins"])
#     # println(simdata)
#     return simdata
# end

# function is

function simdata2df(simdata::SimulationData)
    # builds a data frame from the simulation data

    df = DataFrame(particlelabel=Int[], time=Real[], position=Real[], wrappedposition=Real[], spin=Int[])

    # rows consist of data for 1 particle at 1 time
    times = gettimes(simdata.simparams)
    for t_i in 1:length(simdata.times)
        for i in 1:simparams.numparticles
            push!(df, (i, times[t_i], simdata.positions[t_i, i], simdata.wrappedpositions[t_i, i], simdata.spins[t_i, i]))
        end
    end

    return df
end

function df2simdata(df::DataFrame, simparams::SimulationParameters=SimulationParameters())
    sorteddf = sort(df, [:time, :particlelabel])

    times::Array{Real} = unique(sorteddf.time)
    ntimes = length(times)
    nparticles = maximum(sorteddf.particlelabel)
    positions::Matrix{Real} = zeros(ntimes, nparticles)
    wrappedpositions::Matrix{Real} = zeros(ntimes, nparticles)
    spins::Matrix{Int8} = ones(ntimes, nparticles)

    for time_i in 1:length(times)
        timestamp = times[time_i]
        subdf = filter(:time => t -> t == timestamp, sorteddf)

        if size(subdf, 1) != nparticles
            println("Data for $(size(subdf, 1)) particles found at time $(timestamp). Expected $(nparticles)")
        end
        
        positions[time_i, :] = subdf.position
        wrappedpositions[time_i, :] = subdf.wrappedposition
        spins[time_i, :] = subdf.spin
    end

    simdata = SimulationData(simparams, times, positions, wrappedpositions, spins)

    return simdata

end





function randlinflip(dt::Real, temp::Real)
    return rand() < (dt * temp)
end

function wrap(positions::Array{Real}, boxwidth::Real)
    cwwrap = x -> x - (boxwidth * round(x / boxwidth))
    return cwwrap.(positions)
end

function runstep(currpositions::Array{Real}, currwrappedpositions::Array{Real}, currspins::Array{Int8}, simparams::SimulationParameters)
    # update positions in place
    currpositions .+= (currspins .* simparams.v0 * simparams.dt)

    # handle wrapping of positions (periodic boundary conditions)
    currwrappedpositions = wrap(currpositions, simparams.boxwidth)

    # update spins in place
    flips::Array{Bool} = (x -> randlinflip(simparams.dt, simparams.temp)).(currspins)
    currspins[flips] .*= -1

    # return currpositions, currwrappedpositions, currspins
end


function runsim(simparams::SimulationParameters)

    # set initial particle states
    currpositions::Array{Real} = zeros(1, simparams.numparticles)
    currwrappedpositions::Array{Real} = zeros(1, simparams.numparticles)
    currspins::Array{Int8} = fill!(Array{Int8}(undef, 1, simparams.numparticles), 1)
    
    times = collect(simparams.starttime:simparams.dt:simparams.totaltime)
    nsteps = length(times) - 1

    mx::Matrix{Real} = zeros(nsteps+1, simparams.numparticles)
    mwrappedx::Matrix{Real} = zeros(nsteps+1, simparams.numparticles)
    ms::Matrix{Int8} = zeros(nsteps+1, simparams.numparticles)

    mx[1,:] = copy(currpositions)
    mwrappedx[1,:] = copy(currwrappedpositions)
    ms[1,:] = copy(currspins)

    # run simulation `nsims` times
    for i in 1:nsteps

        runstep(currpositions, currwrappedpositions, currspins, simparams)

        mx[i+1,:] = copy(currpositions)
        mwrappedx[i+1,:] = copy(currwrappedpositions)
        ms[i+1,:] = copy(currspins) 
    end

    simdata = SimulationData(simparams, times, mx, mwrappedx, ms)
    return simdata
end


# get particle densities
function particledensity(simdata::SimulationData, nbins::Real=100)
    n = simdata.simparams.numparticles
    times = simdata.times
    wrappedx = simdata.wrappedpositions

    binwidth = simdata.simparams.boxwidth / nbins
    xgrid::Array{Real} = collect(-simdata.simparams.boxwidth/2:binwidth:simdata.simparams.boxwidth/2)

    densities::Matrix{Real} = zeros(wrappedx.size[1], length(xgrid))

    bin = (x -> round(x / binwidth) + round(length(xgrid) / 2))

    for i in 1:wrappedx.size[1]

        densitymap = counter(bin.(wrappedx[i,:]))

        bins = collect(1:length(xgrid))
        density = (b -> densitymap[b]).(bins)

        densities[i,:] = density ./ n

    end


    return xgrid, densities

end


function set_visual(xlabel = "", ylabel = "",size= (3,2)) # size = (3,2) if figure fills half a page width
    plt.clf() # clears old plots
    fig, ax = plt.subplots(figsize = (size[1],size[2])) # creates plot
    fig = plt.gcf()
    fig.patch.set_alpha(0) # makes area outside of plot transparent

    ax.set_xlabel(xlabel) # sets axis labels
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="y", direction="in", which = "both") # makes ticks face inward
    ax.tick_params(axis="x", direction="in", which = "both")
    return fig,ax # returns figure and axis such that it can be manipulated outside of the function
end

function simpleplotvertxy(simdata::SimulationData, filepath::String = "", title::String = "")
    fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6,10))

    
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(title)   

    for p in 1:simdata.simparams.numparticles
        plt.plot(simdata.positions[:,p], simdata.times)
    end

    plt.xlim(-simdata.simparams.boxwidth / 2, simdata.simparams.boxwidth / 2)
    ax.invert_yaxis()

    # if filepath == ""
    #     filepath = makedefaultfilename()
    # end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 



function densityplot(simdata::SimulationData, filepath::String = "", title::String = "", wrapped::Bool=true, spacebins::Int=100)
    # make particle density plot

    wrappedx = simdata.wrappedpositions
    xlims = [-simdata.simparams.boxwidth / 2, simdata.simparams.boxwidth / 2]
    if !wrapped
        wrappedx = simdata.positions

        m = maximum([abs(minimum(wrappedx)), abs(maximum(wrappedx))])
        xlims = [-m, m]
    end
    xpoints = wrappedx.size[2]
    times = simdata.times
    ntimes = length(times)
    singlelistt = zeros(length(wrappedx))
    singlelistx = zeros(length(wrappedx))

    for i in 1:length(times)
        singlelistt[(i-1)*xpoints+1:(i)*xpoints] = times[i] * ones(xpoints)
        singlelistx[(i-1)*xpoints+1:(i)*xpoints] = wrappedx[i,:]
    end

    # println(singlelistt)
    # println(singlelistx)
    # println(wrappedx)


    fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6,10))

    
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(title)   

    plt.hist2d(singlelistx, singlelistt, bins=[spacebins, length(times)], cmap="Blues")

    plt.xlim(xlims)
    plt.clim([0, simdata.simparams.numparticles / 10])
    ax.invert_yaxis()
    plt.colorbar(label="Density")


    if filepath != ""
        plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

        println("plot saved in current directory as $(filepath)")
    end
end 


# function periodicplotvertxy(simdata::ParticleSimData, filepath::String = "", title::String = "")
#     boundary_crosses = periodicboundary_paths(simdata)
#     println(boundary_crosses)

#     fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6, 10))

    
#     ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    
#     plt.title(title)   

#     for i in 1:length(simdata.simparams.numparticles)
#         path = simdata.positions[:,i]
#         plt.plot(path.positions, path.times)
#     end
        
#     for i in keys(boundary_crosses)
#         path = simdata.particlepaths[i]
#         # Access the color of the last added line
#         line_color = ax.get_lines()[i].get_color()

#         # Add line with the same color to the left and right
#         for j in 1:boundary_crosses[i][1]
#             plot(path.positions .+ (j * simdata.boxwidth), path.times, color=line_color)
#         end
#         for j in 1:boundary_crosses[i][2]
#             plot(path.positions .- (j * simdata.boxwidth), path.times, color=line_color)
#         end
#     end
        
#     plt.xlim(-simdata.boxwidth / 2, simdata.boxwidth / 2)
#     ax.invert_yaxis()

#     if filepath == ""
#         filepath = makedefaultfilename()
#     end

#     plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

#     println("plot saved in current directory as $(filepath)")
# end 



# periodicplotvertxy(simdata, "./plots/spnv$(T)-$(nsims)_per.pdf", plottitle)

# for j in 1:nsims
#     xarray::Array{Real} = []

#     for i in 1:1000
#         timestep!(p, .01)
#         push!(xarray, p.x)
#     end

#     plt.plot(1:length(xarray), xarray)

# end


# plt.savefig("single_particle.pdf", bbox_inches = "tight",pad_inches=0.01)

