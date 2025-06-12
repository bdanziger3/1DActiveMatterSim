using DataFrames

@enum InteractionType nointeraction alignsimple antialignsimple

"""
Simulation parameters on which to run particles
"""
struct SimulationParameters
    numparticles::Int64
    totaltime::Float64
    dt::Float64
    v0::Float64
    fliprate::Float64
    boxwidth::Float64
    interaction::InteractionType
    interactionfliprate::Float64    # probability of a flip for each timestep if in interaction
    starttime::Float64
    randomstarts::Bool


    function # Inner constructor with some default values
        SimulationParameters(numparticles, totaltime, dt, v0, fliprate, boxwidth=1,  interaction=nointeraction, interactionfliprate=Inf64, starttime=0, randomstarts=false)
        return new(numparticles, totaltime, dt, v0, fliprate, boxwidth, interaction, interactionfliprate, starttime, randomstarts)
    end
    function # Inner constructor for dictionary input
        SimulationParameters(paramsdict::Dict)
        return new(paramsdict["numparticles"], paramsdict["totaltime"], paramsdict["dt"], paramsdict["v0"], paramsdict["fliprate"], paramsdict["boxwidth"], paramsdict["interaction"], paramsdict["interactionfliprate"], paramsdict["starttime"], paramsdict["randomstarts"])
    end
    function # Inner constructor for null object with all default values
        SimulationParameters()
        return new(0, 0, 0, 0, 0)
    end
end


function SimulationParameters(simparaminfo::Array)
    return SimulationParameters(simparaminfo...)
end

function SimulationParameters(paramsstr::String)
    # Construct SimulationParameters struct from a comma-separated string
    basic_simparam = SimulationParameters()
    fieldtypes = typeof.(getfield.(Ref(basic_simparam), fieldnames(SimulationParameters)))
    data_list = Array{Any, 1}(undef, length(fieldtypes))
    str_list = split(paramsstr, ",")
    for i = 1:length(fieldtypes)
        if fieldtypes[i] isa Type{<:Enum}
            data_list[i] = eval(Symbol(str_list[i]))
        else
            data_list[i] = parse(fieldtypes[i], str_list[i])
        end
    end

    return SimulationParameters(data_list)
end

# const simparam_fields = ["numparticles", "totaltime", "dt", "v0", "fliprate", "boxwidth", "interaction", "interactionfliprate", "starttime"]
# const simparam_filedtypes = [Int64, Float64, Float64, Float64, Float64, Float64, InteractionType, Float64, Float64]
function asdict(simparams::SimulationParameters)
    paramsdict::Dict = Dict(fieldnames(SimulationParameters) .=> getfield.(Ref(simparams), fieldnames(SimulationParameters)))
    return paramsdict
end

function csv_serialize(simparams::SimulationParameters)
    return join(asarray(simparams), ",")
end

function json_serialize(simparams::SimulationParameters)::String
    # paramsdict::Dict = Dict("numparticles" => simparams.numparticles, "totaltime" => simparams.totaltime, "dt" => simparams.dt, "v0" => simparams.v0, "fliprate" => simparams.fliprate, "boxwidth" => simparams.boxwidth, "starttime" => simparams.starttime)
    # json_str = JSON.json(asdict(simparams))
    # return json_str
    return JSON.json(asdict(simparams), 2)
end

function json_deserialize_simparams(json_str::String)::SimulationParameters
    paramsdict = JSON.parse(json_str)
    simparams = SimulationParameters(paramsdict)
    return simparams
end

function gettimes(simparams::SimulationParameters)::Array{Float64}
    endtime = simparams.starttime + simparams.totaltime
    return collect(simparams.starttime:simparams.dt:endtime)
end

function getntimes(simparams::SimulationParameters)::Int64
    return Int64(floor(simparams.totaltime / simparams.dt) + 1)
end

function asarray(simparams::SimulationParameters)::Array{Any}
    ar::Array{Any} = [simparams.numparticles, simparams.totaltime, simparams.dt, simparams.v0, simparams.fliprate, simparams.boxwidth, simparams.interaction, simparams.interactionfliprate, simparams.starttime, simparams.randomstarts]
    return ar
end

function assertparams(simparams::SimulationParameters)::Bool
    @assert getntimes(simparams) == size(gettimes(simparams))[1] "error in construction of times in `gettimes()`."
    @assert simparams.dt <= simparams.totaltime "timestep `dt` must be less than `totaltime`."

    return true
end

"""
Returns a copy of the `SimulationParameters` struct object but with `randomstarts` set to false and updated `starttime` and `totaltime`
"""
function newstarts(simparams::SimulationParameters, timeextension::Float64)::SimulationParameters
    # get array of original simparams
    sparray = asarray(simparams)

    # update `totaltime` and `starttime`
    sparray[2] = timeextension
    sparray[9] = simparams.totaltime

    # change `randomstarts` to `false`
    sparray[10] = false
    return SimulationParameters(sparray)
end


mutable struct SimulationData
    simparams::SimulationParameters
    times::Array{Float64}
    positions::Matrix{Float64}
    spins::Matrix{Int8}
end

function assertdim(simdata::SimulationData)
    # assert simparams is consistent
    assertparams(simdata.simparams)
    # assert times are correct
    ntimes::Int64 = getntimes(simdata.simparams)
    nparticles::Int64 = simparams.numparticles
    possize = size(simdata.positions)
    spinssize = size(simdata.spins)
    # assert all data arrays have the correct number of times
    @assert ntimes == possize[1] "Positions data has incorrect number of times. Should have $(ntimes), found $(possize[1])"
    @assert ntimes == spinssize[1] "Spins data has incorrect number of times. Should have $(ntimes), found $(spinssize[1])"
    
    # assert all data arrays have the correct number of particles
    @assert nparticles == possize[2] "Positions data has incorrect number of particles. Should have $(nparticles), found $(possize[2])"
    @assert nparticles == spinssize[2] "Spins data has incorrect number of particles. Should have $(nparticles), found $(spinssize[2])"
end

function json_serialize(simdata::SimulationData)
    datadict::Dict = Dict(Pair{String, Dict}("simparams", asdict(simdata.simparams)), Pair{String, Array{Float64}}("times", simdata.times), Pair{String, Matrix{Float64}}("positions", simdata.positions), Pair{String, Matrix{Int8}}("spins", simdata.spins))
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

    df = DataFrame(particlelabel=Int[], time=Real[], position=Real[], spin=Int[])

    # rows consist of data for 1 particle at 1 time
    times = gettimes(simdata.simparams)
    for t_i in 1:length(simdata.times)
        for i in 1:simparams.numparticles
            push!(df, (i, times[t_i], simdata.positions[t_i, i], simdata.spins[t_i, i]))
        end
    end

    return df
end

function df2simdata(df::DataFrame, simparams::SimulationParameters=SimulationParameters())
    sorteddf = sort(df, [:time, :particlelabel])

    times::Array{Float64} = unique(sorteddf.time)
    ntimes = length(times)
    nparticles = maximum(sorteddf.particlelabel)
    positions::Matrix{Float64} = zeros(ntimes, nparticles)
    spins::Matrix{Int8} = ones(ntimes, nparticles)

    for time_i in 1:length(times)
        timestamp = times[time_i]
        subdf = filter(:time => t -> t == timestamp, sorteddf)

        if size(subdf, 1) != nparticles
            println("Data for $(size(subdf, 1)) particles found at time $(timestamp). Expected $(nparticles)")
        end
        
        positions[time_i, :] = subdf.position
        spins[time_i, :] = subdf.spin
    end

    simdata = SimulationData(simparams, times, positions, spins)

    return simdata

end