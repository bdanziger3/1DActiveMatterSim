using CSV
using Dates

include("./basic_sim.jl")
include("../file_management/simdata_files.jl")
include("../utils/paths.jl")
# include("./plot_lib.jl")
# include("./activesim1d.jl")

# import ActiveSim1D
# using ActiveSim1D

NREPS = 5


### Full Sweeps
# flipratesweep = [1, 2, 5, 10, 50, 100, 150, 200, 250, 300]
# flipratesweep = [0.5, 1.5, 2.5, 3, 3.5, 4, 4.5]
# flipratesweep = [0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
# boxwidthsweep = [1, 2, 5, 10, 20, 50, 100]
densitysweep = [.1, .2, .5, 1, 2, 3, 4, 5, 6, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
# densitysweep = [5, 10, 15, 20, 50, 100]

### Single Value
flipratesweep = [100]
boxwidthsweep = [50]
# densitysweep = [10]

serialized = true

fliprate::Real = 1
v0::Real = 1
dt::Real = 1e-4
totaltime::Real = 500
interaction = alignsimple
randomstarts = true
snapshot_dt = 1e-2

segment_splits = 2

local sweeptype = ""
if length(boxwidthsweep) == 1 && length(densitysweep) == 1
    sweeptype = "interactionsweep"
elseif length(boxwidthsweep) == 1 && length(flipratesweep) == 1
    sweeptype = "densitysweep"
elseif length(flipratesweep) == 1 && length(densitysweep) == 1
    sweeptype = "boxwidthsweep"
end

datasweeps_dir = "$(getrootabspath())/work/data/sweeps/$(interaction)/$(sweeptype)"

if !isdir("$(getrootabspath())/work/data/sweeps/$(interaction)")
    mkdir("$(getrootabspath())/work/data/sweeps/$(interaction)")
end
if !isdir(datasweeps_dir)
    mkdir(datasweeps_dir)
end

for i in flipratesweep
    for b in boxwidthsweep
        for d in densitysweep
            for r in 1:NREPS

                # calculate needed number of particles to get density
                n = d * b

                segment_time = Int64(round(totaltime / segment_splits))
                simparams = SimulationParameters(n, segment_time, dt, v0, fliprate, b, interaction, i, 0, randomstarts, snapshot_dt)
                
                data_file = joinpath(datasweeps_dir, "N$(simparams.numparticles)-B$(simparams.boxwidth)-$(simparams.interaction)-$(simparams.interactionfliprate)_rep$r.txt")
                
                
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

                # now collapse the segments into 1 segment
                # collapsesegments(data_file, "$(data_file[1:end-4])_t$totaltime.txt", true)
                collapsesegments(data_file, data_file, true)
                # rm(data_file)
            end
        end
    end
end
