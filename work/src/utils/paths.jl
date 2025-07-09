########################
# paths.jl
# 
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)
# 
# File containing helper methods to get absolute paths from relative paths on generic device
########################



const ROOT_DIR_ABS_PATH1::String = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev"
const ROOT_DIR_ABS_PATH2::String = "/home/s2696227/Documents/Dissertation/1DActiveSolids"

const ROOT_DIRNAME_BD = "Dev"
const ROOT_DIRNAME = "1DActiveSolids"

"""
Gets the absolute path of the root of the 1DActiveSolids Repository on your machine.
"""
function getrootabspath()::String
    currdir::String = @__DIR__

    # look at parent directory until at the root
    while !(endswith(currdir, ROOT_DIRNAME) || endswith(currdir, ROOT_DIRNAME_BD))
        currdir = dirname(currdir)
    end
    return currdir
end


"""
Gets the absolute path on your machine to a file given as a relative path to the root directory `1DActiveSolids`
"""
function relpath2abspath(path::String)::String
    return joinpath(getrootabspath(), path)
end

"""
Takes in an arbitrary path (absolute or relative) and returns the absolute version
"""
function fixpath(path::String)::String
    # check if its an absolute or relative path
    abspathstart = (@__DIR__)[1:7]

    # if it's an absolute path, all good
    if startswith(path, abspathstart)
        return path
    else
        # if not, append it to the absolute of the root path
        return relpath2abspath(path)
    end
end

"""
gets the relative path from the root of an absolute path
"""
function abspath2relpath(path::String)
    return relpath(path, getrootabspath())
end




