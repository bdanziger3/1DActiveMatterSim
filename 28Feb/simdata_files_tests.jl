include("./simdata_files.jl")



# tests for saving and loading data filestring

"""
Runs a short simulation, saves the data, then loads and resaves the data.
Checks that both simdata objects are the same
Checks that both saved files are identical
"""
function test_save_load()
    # TODO: FAILS BECAUSE SOME INITIAL SIMPARAMS ARE LOADED AS INTS.
    # POTENTIAL FIX: FORCE THE SIMPARAMS TO BE FLOATS

    # run simple simulation
    N::Int = 10
    boxwidth::Real = 1
    T::Real = 5
    v0::Real = 1
    dt::Real = 0.1
    totaltime::Real = 10
    simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)
    nsims = 15

    simdata = runsim(simparams)

    # save data and reload it to a new `SimulationData` struct
    datafile_path_1 = "save_load_sequentialtxt_test"
    datafile_path_2 = "save_load_sequentialtxt_test_copy"
    savesim(simdata, datafile_path_1, sequentialtxt, true)
    loadedsim = loadsim("$(datafile_path_1)_simdata.txt", sequentialtxt)

    # save again to a different file
    savesim(loadedsim, datafile_path_2, sequentialtxt, true)

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

test_save_load()