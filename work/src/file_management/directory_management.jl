########################
# directory_management.jl
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# File containing functions to compress, uncompress, and collapse all data files in a directory.
########################

include("./simdata_files.jl")

function compresssimdir(dirname::String, outputdir::String="", deletefiles::Bool=false)

     # load every file in the directory
    datafiles_all = readdir(dirname)
    datafiles = []
    outputfiles = []

    if isempty(outputdir)
        outputdir = dirname
    end

    # only look at .txt files
    for file in datafiles_all
        if isfile(joinpath(dirname, file)) && endswith(file, ".txt")
            push!(datafiles, joinpath(dirname, file))
            push!(outputfiles, joinpath(outputdir, file))
        end
    end

    for (i, filename) in enumerate(datafiles)
        # compressedfilename = "$((outputfiles[i])[1:end-4])_compressed.txt"
        compressedfilename = "$((outputfiles[i])[1:end-17]).txt"
        compresssimfile(filename, compressedfilename)
        
        # delete original file if flag is on
        if deletefiles
            rm(filename)
        end
    end
end


function uncompresssimdir(dirname::String, outputdir::String="", deletefiles::Bool=false)

     # load every file in the directory
    datafiles_all = readdir(dirname)
    datafiles = []
    outputfiles = []

    if isempty(outputdir)
        outputdir = dirname
    end

    # only look at .txt files
    for file in datafiles_all
        if isfile(joinpath(dirname, file)) && endswith(file, ".txt")
            push!(datafiles, joinpath(dirname, file))
            push!(outputfiles, joinpath(outputdir, file))
        end
    end

    for (i, filename) in enumerate(datafiles)
        uncompressedfilename = "$((outputfiles[i])[1:end-4])_uncompressed.txt"
        collapsedfilename = "$((outputfiles[i])[1:end-4])_collapsed.txt"

        # collapse and then compress
        collapsesegments(filename, collapsedfilename)
        uncompresssimfile(collapsedfilename, uncompressedfilename)

        # rm collapsed file
        rm(collapsedfilename)

        # delete original file if flag is on
        if deletefiles
            rm(filename)
        end
    end
end
