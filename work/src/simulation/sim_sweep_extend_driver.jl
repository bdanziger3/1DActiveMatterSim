using CSV
using Dates

include("./basic_sim.jl")
include("../file_management/simdata_files.jl")
include("../utils/paths.jl")

# Extension time per new segment
extensiontime = 900
nsegments = 1
serialized = true

# Directory with simdata files
dir_to_extend = fixpath("work/data/sweeps/alignsimple/interactionsweep")
# dir_to_extend = fixpath("work/data/sweeps/alignsimple/interactionsweep/test1")


# Go through directory
filepahts = joinpath.(dir_to_extend, readdir(dir_to_extend))
filenames = readdir(dir_to_extend)

# #### EXTENDING SIM
for filepath in filepahts

    # extend `nsegments` times if it is a .txt file
    if isfile(filepath)

        for i in 1:nsegments
            extended_sd = extendsim(filepath, extensiontime, serialized)
            appendsim(extended_sd, filepath, serialized)
        end

        # collapse file to one segment
        collapsesegments(filepath, "$(filepath[1:end-4])_t1000.txt")

        # remove original file
        rm(filepath)
    end
end


