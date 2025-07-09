include("./test_fixtures.jl")
include("../src/utils/.paths.jl")



# # tests for saving and loading data filestring
# N::Int = 10
# boxwidth::Real = 1
# T::Real = 5
# v0::Real = 1
# dt::Real = 0.1
# totaltime::Real = 10
# simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)

# snaptshot_dt::Real = 1
# simparams_snapshot = SimulationParameters(N, totaltime, dt, v0, T, boxwidth, nointeraction, 1, 0, false, snaptshot_dt)

"""
Runs a short simulation, saves the data, then loads and resaves the data.
Checks that both simdata objects are the same
Checks that both saved files are identical
"""
function test_save_load()
    # run simple simulation
    simdata = runsim(simparams_10)

    # save data and reload it to a new `SimulationData` struct
    datafile_path_1 = "test_data/save_load_rowwise_test.txt"
    datafile_path_2 = "test_data/save_load_rowwise_test_copy.txt"
    savesim(simdata, datafile_path_1, rowwisetxt, true)
    loadedsim = loadsim(datafile_path_1, rowwisetxt)

    # save again to a different file
    savesim(loadedsim, datafile_path_2, rowwisetxt, true)

    # test simdata objects are identical
    @assert simdata.simparams == loadedsim.simparams "`SimulationParameters` of the two `SimulationData` objects are not identical"
    @assert simdata == loadedsim

    # test files are identical
    linenumber = 0
    open(datafile_path_1, "r") do io1
        open(datafile_path_2, "r") do io2
            while !eof(io1) && !eof(io2)
                line1 = readline(io1)
                line2 = readline(io2)
                linenumber += 1
                if linenumber != 2  # since line 2 has a timestamp and will differ
                    @assert line1 == line2 "Files differ on line $(linenumber)."
                end
            end
            @assert eof(io1) && eof(io2)  "One file ended and the other did not after $(linenumber) lines."
        end
    end

    return true

end

"""
Run an simple simulation (that has `snapshot_dt != dt` and try to load it
    
Make sure the final and original data are the same.
"""
function test_save_load_snaptshot()
    # run simple simulation and try to 
    simdata = runsim(simparams_10_snapshot)

    # save data and reload it to a new `SimulationData` struct
    datafile_path_1 = "test_data/save_load_rowwise_snapshot_test.txt"
    datafile_path_2 = "test_data/save_load_rowwise_snapshot_test_copy.txt"
    savesim(simdata, datafile_path_1, rowwisetxt, true)
    loadedsim = loadsim(datafile_path_1, rowwisetxt)

    # save again to a different file
    savesim(loadedsim, datafile_path_2, rowwisetxt, true)

    # test simdata objects are identical
    @assert simdata.simparams == loadedsim.simparams "`SimulationParameters` of the two `SimulationData` objects are not identical"
    @assert simdata == loadedsim

    # test files are identical
    linenumber = 0
    open(datafile_path_1, "r") do io1
        open(datafile_path_2, "r") do io2
            while !eof(io1) && !eof(io2)
                line1 = readline(io1)
                line2 = readline(io2)
                linenumber += 1
                if linenumber != 2  # since line 2 has a timestamp and will differ
                    @assert line1 == line2 "Files differ on line $(linenumber)."
                end
            end
            @assert eof(io1) && eof(io2)  "One file ended and the other did not after $(linenumber) lines."
        end
    end

    return true

end

function test_save_load_seq()
    # run simple simulation
    simdata = runsim(simparams_10)

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


"""
Load an extended file with 2 or more segments
"""
function test_load_extended_file()
    # make a file with 3 sim segments

    # extended sim file
    datafile = fixpath("/work/tests/test_data/test_extendsim_sd1.txt")
    datafile_long = fixpath("/work/data/22-6/N10000-nointeraction-t100-sn0.01.txt")
    # run simple simulation and try to 
    t0 = time()
    simdata = getsimsegments(datafile)
    t1 = time()

    println(t1 - t0)
    t0 = time()
    simdata = getsimsegments(datafile_long)
    t1 = time()
    println(t1 - t0)


    println(simdata)

    return false

end

"""
1. Run a simulation and load the data.
2. Reload the saved data and extend it
3. Combine the two data structs and append the files together
"""
function test_extendsim_save_load(nsegments::Int=2)
    @assert nsegments >= 2 "Test only works for 2 or more segments."

    timeextension = .0005
    sd1 = runsim(simparams_small)
    sdlist::Array{SimulationData} = Array{SimulationData}[]
    push!(sdlist, sd1)

    filestring = "./test_data/test_extendsim_saveload_sd"
    file1 = "$(filestring).txt"
    savesim(sd1, file1, rowwisetxt, true)

    totaltimesum = sd1.simparams.totaltime
    
    # extend segments to get to `nsegments` total segments
    for s in 2:nsegments
        file_s = "$(filestring)_segment_$(s).txt"
        sd2 = extendsim(file1, timeextension)
        savesim(sd2, file_s, rowwisetxt, true)
        appendsim(sd2, file1)
        push!(sdlist, sd2)

        # check params
        @assert sdlist[s].simparams.starttime == totaltimesum "Segment $(s)'s `startime` is not the sum of the previous segemnts' `totaltime`s. `startime` = $(sdlist[s].simparams.starttime), total time sum = $(totaltimesum)."
        @assert sdlist[s].simparams.totaltime == timeextension "Segment $(s) does not have a `totaltime` equal to the expected extension time."
        @assert sdlist[s].simparams.randomstarts == false "Segment $(s) does not have `randomstarts` set to `false` as expected."
        
        sp1 = asarray(sd1.simparams)
        sp2 = asarray(sd2.simparams)
        
        @assert sp1[1] == sp2[1] "Sim Params of segment $(s) don't match the original for `$(fieldnames(SimulationParameters)[1])`."
        for i in 3:8
            @assert sp1[i] == sp2[i] "Sim Params of segment $(s) don't match the original for `$(fieldnames(SimulationParameters)[i])`."
        end
        
        totaltimesum += sd2.simparams.totaltime
    end

    # now check that the correct data is loaded
    for i in 1:nsegments
        loadedsegment = loadnthsimsegment(file1, i)

        @assert loadedsegment == sdlist[i] "Segment $(i) does not match."
    end

end





# Run the actual tests

# @assert test_save_load()
# @assert test_save_load_snaptshot()
# @assert test_save_load_seq()

# test_load_extended_file()
test_extendsim_save_load(2)
test_extendsim_save_load(3)

println("PASSED")