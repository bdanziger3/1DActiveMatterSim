########################
# convertsimdata.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# Script file that can be executed from the command line to collapse, serialize, and deserialize data files.

# Run from the terminal as follows
# `julia convertsimdata.jl -function sourcefilepath outputfilepath`

# where `-function` is one of `-c`, `-s`, or `-d`, for compress, serialize, or deserialize
########################

include("../simdata_files.jl")

function main(args::Vector{String})
    # check for correct args
    @assert startswith(args[1], "-") "Must specify a function to run: `-c`: 'collapse', `-s`: 'serialize', or `-d`: 'deserialize'."
    @assert length(args) >= 2 "Must provide at least 2 arguments: function to run and source file path. Optional 3rd argument for output file path."
    
    local outputfilepath
    if length(args) == 2
        outputfilepath = findfreefilename((args[2])[1:end-4])
    elseif isfile(args[3])
        outputfilepath = findfreefilename((args[3])[1:end-4])
    else
        outputfilepath = args[3]
    end


    if startswith(args[1], "-c")
        # collapse segments
        collapsesegments(args[2], outputfilepath, true)
    
    elseif startswith(args[1], "-s")
        # serialize / compress file
        compresssimfile(args[2], outputfilepath, true)
    
    elseif startswith(args[1], "-d")
        # deserialize / uncompress file
        uncompresssimfile(args[2], outputfilepath, true)
    end
    
    
    # print out the output file path to stdout
    println(stdout, outputfilepath)
    return outputfilepath
end


main(ARGS)