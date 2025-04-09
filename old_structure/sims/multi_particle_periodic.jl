using Random
using PyPlot
using PyCall
using DelimitedFiles

include("../particle_types/particle.jl")
include("../data_types/SimData.jl")
include("../plotting/plot_template.jl")
# include("../plotting/vertical_graph.jl")


"""
Simulation parameters on which to run particles
"""
struct MultiPeriodicSimulationParameters
    totaltime::Real
    dt::Real
    boxwidth::Real
    ensemble::NonInteractingEnsemble
    # startinglocations::Array{Real}
end


function runsim(sim::MultiPeriodicSimulationParameters)

    # save initial particle params 
    initialstate = sim.ensemble.state

    initialx::Array{Real} = initialstate.x
    initialS::Array{Int8} = initialstate.S
    initialv0::Array{Real} = initialstate.v0
    initialT::Real = sim.ensemble.flipconditionparam

    nparticles = length(initialx)
    
    timesteps = 0:sim.dt:sim.totaltime

    nsteps = length(timesteps)

    xmat = Array{Float64}(undef, length(timesteps), nparticles)
    smat = Array{Float64}(undef, length(timesteps), nparticles)

    # run simulation on array of particle states all at once

    ens = NonInteractingEnsemble(initialx, initialS, initialv0, initialT)

    # update array with location data and move time forward
    for i in 1:nsteps
        # println(xmat)
        # println(smat)
        xmat[i, :] = ens.state.x
        smat[i, :] = ens.state.S
        timestep!(ens, sim.dt)
    end
    # println(xmat)
    # println(smat)


    simdata = EnsembleSimData(xmat, smat, sim.totaltime, sim.dt, sim.boxwidth)
    return simdata
end


# setup simulation params and run the simulation
if length(ARGS) == 0
    nparticles = 2
    dt = 0.5
    L = 2
    totalTime = 1

    v0 = 1
    T = 1
elseif length(ARGS) == 4
    nparticles = parse(Int64, ARGS[1])
    dt = parse(Float64, ARGS[2])
    L = parse(Float64,  ARGS[3])
    totalTime = parse(Float64, ARGS[4])

    v0 = 1
    T = 1
elseif length(ARGS) == 6
    nparticles = parse(Int64, ARGS[1])
    dt = parse(Float64, ARGS[2])
    L = parse(Float64,  ARGS[3])
    totalTime = parse(Float64, ARGS[4])

    v0 = parse(Float64, ARGS[5])
    T = parse(Float64, ARGS[6])
else
    println("Unexpected number of arguments. Using default values.\nRun with either 4 or 6 arguments")
    nparticles = 2
    dt = 0.5
    L = 2
    totalTime = 1

    v0 = 1
    T = 1
end


x = zeros(nparticles)
S = ones(nparticles)
v0 = ones(nparticles)

sim = MultiPeriodicSimulationParameters(totalTime, dt, L, NonInteractingEnsemble(x, S, v0, T))
simdata = runsim(sim)

plottitle = "$(nparticles) Particles \n T = $(T), dt = $(dt), v0 = $(v0)"

print(simdata.xmat)
print(simdata.smat)
writedlm("./sims/sim_output_data/multi_sim_$(nparticles)-$(T)-$(v0[1])_x.txt", simdata.xmat, ",")
writedlm("./sims/sim_output_data/multi_sim_$(nparticles)-$(T)-$(v0[1])_S.txt", simdata.smat, ",")


