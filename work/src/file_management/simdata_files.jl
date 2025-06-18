include("../simulation/basic_sim.jl")
# include("./sim_structs.jl")

using JSON
using Dates

@enum DataFileType seperatefiles sequentialtxt rowwisetxt json

function getfirstsaveddate(filename::String)
    file = open(filename, "r")
    firstline = readline(file)
    dateline = readline(file)
    saveddate::DateTime = DateTime(0)
        
    # if the current line has a date, read it
    datestartindex = length("Saved at ") + 1
    try 
        saveddate = DateTime(dateline[datestartindex:end])
    catch e
        saveddate = DateTime(0)
    end
    close(file)

    return saveddate
end

function checkfilename(outputfilename::String)
    if isfile(outputfilename)
        saveddate::DateTime = getfirstsaveddate(outputfilename)
        println("Overwrite $(outputfilename) with simulation data saved at $(saveddate)? Y/N ")

        response = readline()
        if lowercase(response) != "y"
            # find new file_name
            fileprefix = outputfilename[1:end-4]
            i = 0
            while isfile("$(fileprefix)_$(i).txt")
                i += 1
            end
            newname = "$(fileprefix)_$(i).txt"
            println("Saving instead as $(newname)")
            return newname
        end
    end

    # Use original name if file doesn't exist or user confirms overwrite
    return outputfilename
end

function loadsim(inputfilename::String, filetype::DataFileType)::SimulationData
    # initialize parameters and matrices        
    ntimes::Int64 = 0
    posmatrix::Matrix{Float64} = zeros(0, 0)
    spinsmatrix::Matrix{Int8} = zeros(0, 0)
    simdates::Array{DateTime} = Array{DateTime}[]

    if filetype == rowwisetxt
        # read file
        file = open(inputfilename, "r")

        line = readline(file)
        # Check that data file header is correct
        @assert (contains(lowercase(line), "row wise txt") || contains(lowercase(line), "rowwise txt")) "Rowwisetxt data file does not have correct first line. File may be corrupted."
        
        line = readline(file)

        # if the current line is a date, read it. Else skip to next line
        datestartindex = length("Saved at ") + 1
        try 
            push!(simdates, DateTime(line[datestartindex:end]))
        catch e
            push!(simdates, DateTime(0))
        else
            line = readline(file)
        end

        @assert contains(lowercase(line), "simulation parameters") "Rowwisetxt data file does not have header. File may be corrupted."
       
        # read next line to get number of position data points
        simparaminfo_str = readline(file)
        simparams = SimulationParameters(simparaminfo_str)
        ntimes = getnsaves(simparams)
            
        # now read through the data in the order: [[positions], spins]]
        line = readline(file)
        @assert contains(lowercase(line), "particle states") "Rowwisetxt data file does not have correct `Particle States` header. File may be corrupted."
        
        # positions and spins
        # read in data from file
        posmatrix = zeros(ntimes, simparams.numparticles)
        spinsmatrix = zeros(ntimes, simparams.numparticles)
        for i in 1:ntimes                    
            dataline = readline(file)
            particle_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
            posmatrix[i,:] = copy(particle_i_data[1:simparams.numparticles])
            spinsmatrix[i,:] = Int8.(copy(particle_i_data[simparams.numparticles+1:end]))
        end

        # close file
        close(file)
            
        # store matrix data as `SimulationData` object and return
        simdata = SimulationData(simparams, gettimes(simparams), posmatrix, spinsmatrix)
        return simdata

    elseif filetype == sequentialtxt
        # read file
        file = open(inputfilename, "r")

        line = readline(file)
        # skip first line if contains header
        if contains(lowercase(line), "sequential txt")
            line = readline(file)
        end
        
        # Check that data file header is correct
        @assert contains(lowercase(line), "simulation parameters") "Sequentialtxt data file does not have correct header. File may be corrupted."
       
        # read next line to get number of position data points
        simparaminfo_str = readline(file)
        simparams = SimulationParameters(simparaminfo_str)
        ntimes = getnsaves(simparams)
            
        # now read through the data in the order: [positions, wrappedpositions, spins]
        line = readline(file)
        @assert contains(lowercase(line), "positions") "Sequentialtxt data file does not have correct `Positions` header. File may be corrupted."
        
        # read in position data from file
        posmatrix = zeros(ntimes, simparams.numparticles)
        for i in 1:ntimes                    
            dataline = readline(file)
            pos_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
            posmatrix[i,:] = copy(pos_i_data)
        end
        
        # read in spins data from file
        line = readline(file)
        @assert contains(lowercase(line), "spins") "Sequentialtxt data file does not have correct `Spins` header line. File may be corrupted."
        
        spinsmatrix = zeros(ntimes, simparams.numparticles)
        for i in 1:ntimes                    
            dataline = readline(file)
            spins_i_data::Array{Int8} = parse.(Int8, split(dataline, ","))
            spinsmatrix[i,:] = copy(spins_i_data)
        end

        # close file
        close(file)
            
        # store matrix data as `SimulationData` object and return
        simdata = SimulationData(simparams, gettimes(simparams), posmatrix, spinsmatrix)
        return simdata
    end
end


"""
Reads only the `nlines` lines of data in a data file starting at line `start_line`.
`startline` is indexed at 1

use `startline=-1` to read the last `nlines` lines of the file

"""
function loadsim_nlines(inputfilename::String, startline::Int, nlines::Int, filetype::DataFileType=rowwisetxt)::SimulationData
    # initialize parameters and matrices        
    ntimes::Int64 = nlines
    posmatrix::Matrix{Float64} = zeros(nlines, 0)
    spinsmatrix::Matrix{Int8} = zeros(nlines, 0)
    simdates::Array{DateTime} = Array{DateTime}[]

    # read file
    file = open(inputfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert (contains(lowercase(line), "row wise txt") || contains(lowercase(line), "rowwise txt")) "Rowwisetxt data file does not have correct first line. File may be corrupted."
    
    line = readline(file)
    # if the current line is a date, read it. Else skip to next line
    datestartindex = length("Saved at ") + 1
    try 
        push!(simdates, DateTime(line[datestartindex:end]))
    catch e
        push!(simdates, DateTime(0))
    else
        line = readline(file)
    end
    @assert contains(lowercase(line), "simulation parameters") "Rowwisetxt data file does not have header. File may be corrupted."
    
    # read next line to get number of position data points
    simparaminfo_str = readline(file)
    simparams = SimulationParameters(simparaminfo_str)
    ntimes = getnsaves(simparams)
        
    # now read through the data in the order: [[positions], spins]]
    line = readline(file)
    @assert contains(lowercase(line), "particle states") "Rowwisetxt data file does not have correct `Particle States` header. File may be corrupted."
    
    # positions and spins

    # skip to start line
    if startline == -1
        linestoskip = ntimes - nlines
    else
        linestoskip = startline - 1
    end

    for _ in 1:linestoskip
        readline(file)
    end

    # read in `nlines` lines of data from file
    posmatrix = zeros(nlines, simparams.numparticles)
    spinsmatrix = zeros(nlines, simparams.numparticles)
    for i in 1:nlines                    
        dataline = readline(file)
        particle_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
        posmatrix[i,:] = copy(particle_i_data[1:simparams.numparticles])
        spinsmatrix[i,:] = Int8.(copy(particle_i_data[simparams.numparticles+1:end]))
    end

    # close file
    close(file)
        
    # store matrix data as `SimulationData` object and return
    simdata = SimulationData(simparams, gettimes(simparams), posmatrix, spinsmatrix)
    return simdata
end

"""
Saves simulation data from a `SimulationData` object to a data file.

Can specify the output file type with `filetype`.

By default, asks for confirmation if overwriting an existing file.
Set `forceclear` to `true` to clear the file without warning.

"""
function savesim(simdata::SimulationData, outputfilename::String, filetype::DataFileType=rowwisetxt, forceclear::Bool=false)         
    # Ask to confirm if overwriting existing file
    if !forceclear
        filenametouse = checkfilename(outputfilename)
    else
        filenametouse = outputfilename
    end
    
    if filetype == sequentialtxt
        filestring = "$(filenametouse[1:end-4])_seq.txt"
        open(filestring, "w") do io
            println(io, "Sequential txt")   # Sequential .txt
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simdata.simparams))
            println(io, "Positions")
            writedlm(io, simdata.positions, ",")
            println(io, "Spins")
            writedlm(io, simdata.spins, ",")
        end
    elseif filetype == rowwisetxt
        open(filenametouse, "w") do io
            println(io, "Row Wise txt")  # Row Wise .txt
            println(io, "Saved at $(round(now(), Dates.Second(1)))")  # Timestamp
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simdata.simparams))
            println(io, "Particle States ([positions], [spins])")
            # construct particle state matrix by concatenating states of position and spin
            particlestates::Matrix{Any} = zeros(getnsaves(simdata.simparams), 2 * simdata.simparams.numparticles)
            particlestates[:,1:simdata.simparams.numparticles] = simdata.positions
            particlestates[:,simdata.simparams.numparticles+1:end] = simdata.spins
            writedlm(io, particlestates, ",")
        end
    end

end


"""
Saves simulation data from a `SimulationData` object to an existing data file.

Appends to the bottom of an existing file with a new header for the timestamp and Simulation Parameters
"""
function appendsim(simdata::SimulationData, outputfilename::String)         
    # Only works with Row Wise txt files
    # Read file to make sure it's the right type
    file = open(outputfilename, "r")
    firstline = readline(file)
    @assert contains(lowercase(firstline), "row wise txt") "Incorrect file header on first line. File not identified as Row Wise txt"
    close(file)

    # now append new data to the end
    open(outputfilename, "a") do io
        println(io, "Saved at $(round(now(), Dates.Second(1)))")  # Timestamp
        println(io, "Simulation Parameters")
        println(io, csv_serialize(simdata.simparams))
        println(io, "Particle States ([positions], [spins])")
        # construct particle state matrix by concatenating states of position and spin
        particlestates::Matrix{Any} = zeros(getnsaves(simdata.simparams, simdata.simparams.dt), 2 * simdata.simparams.numparticles)
        particlestates[:,1:simdata.simparams.numparticles] = simdata.positions
        particlestates[:,simdata.simparams.numparticles+1:end] = simdata.spins
        writedlm(io, particlestates, ",")
    end
end

