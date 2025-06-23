include("../src/file_management/simdata_files.jl")
include("../src/file_management/data_serialisation.jl")


# tests for saving and loading data filestring
N::Int = 100
boxwidth::Real = 1
fliprate::Real = 1
v0::Real = 1
dt::Real = 0.1
totaltime::Real = 1
simparams = SimulationParameters(N, totaltime, dt, v0, fliprate, boxwidth)

"""
Serialises and deserialises the spins and positions of a simulation and confirms that the values are the same
"""
function test_serialisation()
    # # run simple simulation if data file does not exist
    # datafile_path = "test_data/serialisation_test_data.txt"
    # if !isfile(datafile_path)
    #     simdata = runsim(simparams)
    #     savesim(simdata, datafile_path, rowwisetxt, true)
    # end

    # # load data and try to serialise the spins
    # loadedsim::SimulationData = loadsim(datafile_path, rowwisetxt)


    spins = rand([1,-1], N)
    # spins = [1, 1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1]

    println(spins)

    t0 = time()
    serialisedstr_ascii = serialisespins(spins, ascii128)
    t1 = time()

    println(serialisedstr_ascii)

    td0 = time()
    newspins_ascii = deserialisespins(serialisedstr_ascii, N, ascii128)
    td1 = time()

    @assert newspins_ascii == spins "ASCII deserialised spins are not the same as the initial spins"

    println("Timing results:\nASCII: encoding: $(t1 - t0) seconds, decoding: $(td1 - td0) seconds")
    println(td1)
end

    



test_serialisation()