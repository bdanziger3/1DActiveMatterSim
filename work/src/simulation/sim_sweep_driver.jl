using CSV
using Dates

include("./basic_sim.jl")
include("../file_management/simdata_files.jl")
include("../utils/paths.jl")
# include("./plot_lib.jl")
# include("./activesim1d.jl")

# import ActiveSim1D
# using ActiveSim1D


### Full Sweeps
flipratesweep = [10, 50, 100, 150, 200, 250, 300]
# boxwidthsweep = [25, 50, 75, 100, 125, 150, 175, 200]
# densitysweep = [.1, .5, 1, 5, 10, 15, 20, 50, 100]
# densitysweep = [5, 10, 15, 20, 50, 100]

### Single Value
# flipratesweep = [100]
boxwidthsweep = [100]
densitysweep = [100]

serialized = true

fliprate::Real = 1
v0::Real = 1
dt::Real = 1e-4
totaltime::Real = 100
interaction = turnaway
randomstarts = true
snapshot_dt = 1e-2

segment_splits = 1


datasweeps_dir = "$(getrootabspath())/work/data/sweeps/$(interaction)/interactionsweep"


for i in flipratesweep
    for b in boxwidthsweep
        for d in densitysweep

            # calculate needed number of particles to get density
            n = d * b

            segment_time = Int64(round(totaltime / segment_splits))
            simparams = SimulationParameters(n, segment_time, dt, v0, fliprate, b, interaction, i, 0, randomstarts, snapshot_dt)
            
            data_file = joinpath(datasweeps_dir, "N$(simparams.numparticles)-B$(simparams.boxwidth)-$(simparams.interaction)-$(Int64(round(simparams.interactionfliprate))).txt")
            
            
            # run first sim to create the file
            simdata = runsim(simparams)
            savesim(simdata, data_file, serialized ? rwtserialized : rowwisetxt)

            println("Saved first $(segment_time)s at $(now())")
            for s in 1:segment_splits-1
                #### EXTENDING SIM
                extended_sd = extendsim(data_file, segment_time, serialized)
                appendsim(extended_sd, data_file, serialized)
                println("Saved extended sim $(s) at $(now())")
            end

        end
    end
end
