########################
# serialization_driver.jl
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing driver code to condense data file segments and serialize/deserialize data files
########################

include("./simdata_files.jl")
include("../utils/paths.jl")


"""
Collapses all data files in a given directory so that each only has one segment of data
Set `deleteoriginals=true` 
"""
function collapsefilesindir(dirname::String, deleteoriginals::Bool=false)

    # open directory and collapse all files
    filelist = readdir(dirname)

    for filename in filelist
        inputfilename = joinpath(dirname, filename)
        outputfilename = "$(inputfilename[1:end-4])_collapsed.txt"
        collapsesegments(inputfilename, outputfilename)

        # delete original if toggled
        if deleteoriginals
            rm(inputfilename)
        end
    end
end


"""
Collapses all data files in a given directory so that each only has one segment of data
Set `deleteoriginals=true` 
"""
function uncompressfiledir(dirname::String, deleteoriginals::Bool=false)

    # open directory and collapse all files
    filelist = readdir(dirname)

    for filename in filelist
        inputfilename = joinpath(dirname, filename)
        outputfilename = "$(inputfilename[1:end-4])_deserialized.txt"
        
        try
            uncompresssimfile(inputfilename, outputfilename)
        catch e
            println("Could not uncompress $(inputfilename)")
        else
            # delete original if toggled
            if deleteoriginals
                rm(inputfilename)
            end
        end

    end
end


datadir = fixpath("work/data/sweeps/densitysweep")
datadir_collapsed = fixpath("work/data/sweeps/densitysweep/uncompressed")



collapsefilesindir(datadir)

# uncompressfiledir(datadir_collapsed, true)
