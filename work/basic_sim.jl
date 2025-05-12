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



# Minumum approach to the same code I had before



function wrap(positions::AbstractArray{Float64}, boxwidth::Float64)
    cwwrap = x -> x - (boxwidth * round(x / boxwidth))
    return cwwrap.(positions)
end

function runstep!(currpositions::Array{Float64}, currwrappedpositions::Array{Float64}, currspins::Array{Int8}, simparams::SimulationParameters)
    # update positions in place
    currpositions .+= (currspins .* simparams.v0 * simparams.dt)

    # handle wrapping of positions (periodic boundary conditions)
    newcurrwrappedpositions = wrap(currpositions, Float64(simparams.boxwidth))
    currwrappedpositions .= newcurrwrappedpositions

    # update spins in place
    flips::Array{Bool} = applyinteraction(simparams, currspins, currwrappedpositions)
    # flips::Array{Bool} = (x -> randlinflip(simparams.dt, simparams.fliprate)).(currspins)
    currspins[flips] .*= -1

    # return currpositions, currwrappedpositions, currspins
end


function runsim(simparams::SimulationParameters)

    # set initial particle states
    currpositions::Array{Float64} = zeros(1, simparams.numparticles)
    currwrappedpositions::Array{Float64} = zeros(1, simparams.numparticles)
    currspins::Array{Int8} = fill!(Array{Int8}(undef, 1, simparams.numparticles), 1)
    
    times = collect(simparams.starttime:simparams.dt:simparams.totaltime)
    nsteps = length(times) - 1

    mx::Matrix{Float64} = zeros(nsteps+1, simparams.numparticles)
    mwrappedx::Matrix{Float64} = zeros(nsteps+1, simparams.numparticles)
    ms::Matrix{Int8} = zeros(nsteps+1, simparams.numparticles)

    mx[1,:] = copy(currpositions)
    mwrappedx[1,:] = copy(currwrappedpositions)
    ms[1,:] = copy(currspins)

    # run simulation `nsims` times
    @showprogress 1 "Loading..." for i in 1:nsteps
    # for i in 1:nsteps

        runstep!(currpositions, currwrappedpositions, currspins, simparams)

        mx[i+1,:] = copy(currpositions)
        mwrappedx[i+1,:] = copy(currwrappedpositions)
        ms[i+1,:] = copy(currspins) 
    end

    simdata = SimulationData(simparams, times, mx, mwrappedx, ms)
    return simdata
end

