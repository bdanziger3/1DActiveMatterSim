include("../src/file_management/simdata_files.jl")



# tests for saving and loading data filestring
N::Int = 10
boxwidth::Real = 1
T::Real = 5
v0::Real = 1
dt::Real = 0.1
totaltime::Real = 10
simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)

"""
Runs a short simulation, saves the data, then loads and resaves the data.
Checks that both simdata objects are the same
Checks that both saved files are identical
"""
function test_save_load()
    # run simple simulation
    simdata = runsim(simparams)

    # save data and reload it to a new `SimulationData` struct
    datafile_path_1 = "test_data/save_load_rowwise_test.txt"
    datafile_path_2 = "test_data/save_load_rowwise_test_copy.txt"
    savesim(simdata, datafile_path_1, rowwisetxt, true)
    loadedsim = loadsim(datafile_path_1, rowwisetxt)

    # save again to a different file
    savesim(loadedsim, datafile_path_2, rowwisetxt, true)

    # test simdata objects are identical
    # @assert simdata.simparams == loadedsim.simparams "`SimulationParameters` of the two `SimulationData` objects are not identical"
    println(simdata.simparams)
    println(loadedsim.simparams)
    @assert simdata == loadedsim

    # test files are identical
    linenumber = 0
    open(datafile_path_1, "r") do io1
        open(datafile_path_2, "r") do io2
            while !eof(io1) && !eof(io2)
                line1 = readline(io1)
                line2 = readline(io2)
                linenumber += 1
                @assert line1 == line2 "Files differ on line $(linenumber)."
            end
            @assert eof(io1) && eof(io2)  "One file ended and the other did not after $(linenumber) lines."
        end
    end

    return true

end

function test_save_load_seq()
    # run simple simulation
    simdata = runsim(simparams)

    # save data and reload it to a new `SimulationData` struct
    datafile_path_1 = "test_data/save_load_sequentialtxt_test.txt"
    datafile_path_1_saved_as = "test_data/save_load_sequentialtxt_test_seq.txt"
    datafile_path_2 = "test_data/save_load_sequentialtxt_test_copy.txt"
    datafile_path_2_saved_as = "test_data/save_load_sequentialtxt_test_copy_seq.txt"
    savesim(simdata, datafile_path_1, sequentialtxt, true)
    loadedsim = loadsim(datafile_path_1_saved_as, sequentialtxt)

    # save again to a different file
    savesim(loadedsim, datafile_path_2, sequentialtxt, true)

    # test simdata objects are identical
    # @assert simdata.simparams == loadedsim.simparams "`SimulationParameters` of the two `SimulationData` objects are not identical"
    println(simdata.simparams)
    println(loadedsim.simparams)
    @assert simdata == loadedsim

    # test files are identical
    linenumber = 0
    open(datafile_path_1_saved_as, "r") do io1
        open(datafile_path_2_saved_as, "r") do io2
            while !eof(io1) && !eof(io2)
                line1 = readline(io1)
                line2 = readline(io2)
                linenumber += 1
                @assert line1 == line2 "Files differ on line $(linenumber)."
            end
            @assert eof(io1) && eof(io2)  "One file ended and the other did not after $(linenumber) lines."
        end
    end

    return true

end

test_save_load()