include("../src/simulation/basic_sim.jl")
include("../src/file_management/simdata_files.jl")
# include("./plot_lib.jl")
# include("./activesim1d.jl")

# import ActiveSim1D
# using ActiveSim1D

data_dir = "./work/data"


# println(length(simdata.positions[:,1]))
# println(length(simdata.times))
# setup simulation params and run the simulation
N::Int = 3
boxwidth::Real = 1
fliprate::Real = 1
v0::Real = 1 #.5 * boxwidth / 3
dt::Real = 0.0001
totaltime::Real = .001
interaction = alignsimple
interactionfliprate = 300
randomstarts = true

simparams_small = SimulationParameters(N, totaltime, dt, v0, fliprate, boxwidth, interaction, interactionfliprate, 0, randomstarts)


# longer tests for saving and loading data filestring
N::Int = 10
boxwidth::Real = 1
T::Real = 1
v0::Real = 1
dt::Real = 0.1
totaltime::Real = 10
simparams_10 = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)

snaptshot_dt::Real = 1
simparams_10_snapshot = SimulationParameters(N, totaltime, dt, v0, T, boxwidth, nointeraction, 1, 0, false, snaptshot_dt)
