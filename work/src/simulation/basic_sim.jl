import Base.==
using Random
using PyPlot
using PyCall
using DelimitedFiles
using DataStructures
using DataFrames
using ProgressMeter

include("./sim_structs.jl")
include("./interactions.jl")




# Minumum approach to the 1D Simulations

# Wrap positions as the simulation runs
function wrap(positions::AbstractArray{Float64}, boxwidth::Float64)::AbstractArray{Float64}
    cwwrap = x -> x - (boxwidth * round(x / boxwidth))
    return cwwrap.(positions)
end

function runstep!(currpositions::Array{Float64}, currspins::Array{Int8}, simparams::SimulationParameters)::Nothing
    # update positions in place
    currpositions .+= (currspins .* simparams.v0 * simparams.dt)

    # handle wrapping of positions (periodic boundary conditions)
    currpositions .= wrap(currpositions, Float64(simparams.boxwidth))
    # update spins in place
    flips::Array{Bool} = calcspinflips(simparams, currspins, currpositions)
    currspins[flips] .*= -1

    return nothing
end


function runsim(simparams::SimulationParameters, startpositions::Union{Array{Float64}, Nothing}=nothing, startspins::Union{Array{Int8}, Nothing}=nothing)::SimulationData
    # set initial particle states
    # if no start positions specified, start all at 0
    # if no start spins specified, start all at +1
    currpositions::Array{Float64} = zeros(1, simparams.numparticles)
    currspins::Array{Int8} = fill!(Array{Int8}(undef, 1, simparams.numparticles), 1)
    if simparams.randomstarts
        for i in 1:simparams.numparticles
            currpositions[i] = (rand() - 0.5) * simparams.boxwidth
            currspins[i] = rand([1,-1])
        end
    else
        if !isnothing(startpositions)
            currpositions = copy(startpositions)
        end
        if !isnothing(startspins)
            currspins = copy(startspins)
        end
    end

    
    endtime = simparams.starttime + simparams.totaltime
    times = collect(simparams.starttime:simparams.dt:endtime)
    nsteps = length(times) - 1

    stepsbetweensaves = Int64(round(simparams.snapshot_dt / simparams.dt))
    nsaves = Int64(floor(nsteps / stepsbetweensaves))

    mx::Matrix{Float64} = zeros(nsaves+1, simparams.numparticles)
    ms::Matrix{Int8} = zeros(nsaves+1, simparams.numparticles)

    mx[1,:] = copy(currpositions)
    ms[1,:] = copy(currspins)

    # run simulation `nsims` times
    stepsuntilsave = stepsbetweensaves
    saves = 1
    @showprogress 1 "Running Simulation (N=$(simparams.numparticles), nsteps=$(nsteps))..." for i in 1:nsteps

        runstep!(currpositions, currspins, simparams)

        stepsuntilsave -= 1

        if stepsuntilsave == 0
            mx[saves+1,:] = copy(currpositions)
            ms[saves+1,:] = copy(currspins)

            stepsuntilsave = stepsbetweensaves
            saves += 1
        end
    end

    simdata = SimulationData(simparams, times, mx, ms)
    return simdata
end


function extendsim(inputfilename::String, time::Number)::SimulationData
    finalstate::SimulationData = loadsim_nlines(inputfilename, -1, 1)
    # get final state of existing file

    newsimparams = newstarts(finalstate.simparams, time)

    extended_simdata = runsim(newsimparams, finalstate.positions, finalstate.spins)

    return extended_simdata
    
    # run simulation starting from this state

    # save new data by appending it to existing file

end

