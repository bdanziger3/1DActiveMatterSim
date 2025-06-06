import Base.==
using Random
using PyPlot
using PyCall
using DelimitedFiles
using DataStructures
using DataFrames
using ProgressMeter

include("./interactions.jl")
include("./sim_structs.jl")


# Minumum approach to the 1D Simulations

# Wrap positions as the simulation runs
function wrap(positions::AbstractArray{Float64}, boxwidth::Float64)
    cwwrap = x -> x - (boxwidth * round(x / boxwidth))
    return cwwrap.(positions)
end

function runstep!(currpositions::Array{Float64}, currspins::Array{Int8}, simparams::SimulationParameters)
    # update positions in place
    currpositions .+= (currspins .* simparams.v0 * simparams.dt)


    # handle wrapping of positions (periodic boundary conditions)
    old_p1 = currpositions[1]
    currpositions .= wrap(currpositions, Float64(simparams.boxwidth))
    # update spins in place
    flips::Array{Bool} = applyspininteraction(simparams, currspins, currpositions)
    currspins[flips] .*= -1
end


function runsim(simparams::SimulationParameters)

    # set initial particle states
    # currpositions::Array{Float64} = initialpositions(simparams)
    currpositions::Array{Float64} = zeros(1, simparams.numparticles)
    currspins::Array{Int8} = fill!(Array{Int8}(undef, 1, simparams.numparticles), 1)
    
    times = collect(simparams.starttime:simparams.dt:simparams.totaltime)
    nsteps = length(times) - 1

    mx::Matrix{Float64} = zeros(nsteps+1, simparams.numparticles)
    ms::Matrix{Int8} = zeros(nsteps+1, simparams.numparticles)

    mx[1,:] = copy(currpositions)
    ms[1,:] = copy(currspins)

    # run simulation `nsims` times
    @showprogress 1 "Loading..." for i in 1:nsteps
    # for i in 1:nsteps

        runstep!(currpositions, currspins, simparams)

        mx[i+1,:] = copy(currpositions)
        ms[i+1,:] = copy(currspins) 
    end

    simdata = SimulationData(simparams, times, mx, ms)
    return simdata
end

