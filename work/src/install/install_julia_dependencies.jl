########################
# install_julia_dependencies.jl
# 
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)
# 
# Finds all Julia packages listed in `dependeencies.md` and installs them.
########################

# ENV["PYTHON"] = "/home/s2696227/python310/bin/python"
using Pkg

DEPENDENCIES_FILENAME = "dependencies.md"
full_dependencies_filepath = joinpath(dirname(dirname(dirname(dirname(@__FILE__)))), DEPENDENCIES_FILENAME)

println(full_dependencies_filepath)
function main()
    """Opens the dependencies file and installs each of the Python Packages listed"""
    # open file and get all lines
    dependencies_file = open(full_dependencies_filepath)
    lines = readlines(dependencies_file)
    close(dependencies_file)

    reading_pkgs = false
    for (i, line) in enumerate(lines)
        if contains(lowercase(line), "julia packages")
            # turn on flags to start reading
            reading_pkgs = true

        elseif reading_pkgs

            # turn off flag at end of list
            if !startswith(line, "-")
                reading_pkgs = false
            else
                # get package name and install
                pkginfo = split(strip(line[2:end]), " ")
                pkgname = pkginfo[1]
                pkgversion = pkginfo[2]

                # build PyCall
                Pkg.add(pkgname)
            end
        end
    end
end



main()