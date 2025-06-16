using CSV

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
ctime::Real = 100
fliprate::Real = 1
v0::Real = 1 #.5 * boxwidth / 3
dt::Real = 0.0001
totaltime::Real = .001
interaction = alignsimple
interactionfliprate = 300
randomstarts = true
simparams = SimulationParameters(N, totaltime, dt, v0, fliprate, boxwidth, interaction, interactionfliprate, 0, randomstarts)
# nsims = 15


"""
1. Run a simulation and load the data.
2. Reload the saved data and extend it
3. Combine the two data structs and append the files together
"""
function test_extendsim()
    timeextension = .0005
    sd1 = runsim(simparams)

    filestring = "./test_data/test_extendsim_sd1"
    file1 = "$(filestring).txt"
    file2 = "$(filestring)_ext.txt"
    savesim(sd1, file1, rowwisetxt)
    
    sd2 = extendsim(file1, timeextension)
    savesim(sd2, file2, rowwisetxt)

    # check params
    @assert sd2.simparams.starttime == sd1.simparams.totaltime
    @assert sd2.simparams.totaltime == timeextension
    @assert sd2.simparams.randomstarts == false

    sp1 = asarray(sd1.simparams)
    sp2 = asarray(sd2.simparams)

    @assert sp1[1] == sp2[1]
    for i in 3:8
        @assert sp1[i] == sp2[i]
    end


    appendsim(sd2, file1)
end


test_extendsim()