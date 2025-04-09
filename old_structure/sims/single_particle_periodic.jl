using Random
using PyPlot
using PyCall

include("../particle_types/particle.jl")
include("../data_types/SimData.jl")
include("../plotting/plot_template.jl")
include("../plotting/vertical_graph.jl")


"""
Simulation parameters on which to run particles
"""
struct PeriodicSimulationParameters
    totaltime::Real
    dt::Real
    boxwidth::Real
    particlelist::Array{NonInteractingParticle}
    # startinglocations::Array{Real}
end


function runsim(sim::PeriodicSimulationParameters, nsims::Int64)

    # save initial particle params 
    initialstate = sim.particlelist[1].state

    initialx::Real = initialstate.x
    initialS::Int8 = initialstate.S
    initialv0::Real = initialstate.v0
    initialT::Real = sim.particlelist[1].flipconditionparam
    
    timesteps = 0:sim.dt:sim.totaltime

    simdata = SimData(ParticlePath[], sim.totaltime, sim.dt, sim.boxwidth)

    # run simulation `nsims` times
    for i in 1:nsims
        xarray::Array{Real} = []

        p = NonInteractingParticle(initialx, initialS, initialv0, initialT)

        # update array with location data and move time forward
        for _ in timesteps
            push!(xarray, p.state.x)
            timestep!(p, sim.dt)
        end

        # create and store ParticlePath to output data
        path = ParticlePath("Particle $(i)", timesteps, xarray)
        push!(simdata.particlepaths, path)
    end

    

    matplotlib.pyplot.close()

    return simdata
end

# setup simulation params and run the simulation
boxwidth = 2
T = 1
dt = 0.1
totalTime = 5
sim = PeriodicSimulationParameters(totalTime, dt, boxwidth, [NonInteractingParticle(0, 1, 1, T)])
nsims = 15
simdata = runsim(sim, nsims)

plottitle = "One Particle \n T = $(T), dt = $(dt)"
# handle plotting
simpleplotvertxy(simdata, "./plots/spnv$(T)-$(nsims).pdf", "Hello")
periodicplotvertxy(simdata, "./plots/spnv$(T)-$(nsims)_per.pdf", plottitle)

# for j in 1:nsims
#     xarray::Array{Real} = []

#     for i in 1:1000
#         timestep!(p, .01)
#         push!(xarray, p.x)
#     end

#     plt.plot(1:length(xarray), xarray)

# end


# plt.savefig("single_particle.pdf", bbox_inches = "tight",pad_inches=0.01)

