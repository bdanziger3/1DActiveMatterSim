include("../simulation/basic_sim.jl")
# include("./sim_structs.jl")

using JSON
using Dates

const DATE_START_INDEX = length("Saved at ") + 1    # 9
const SEGMENT_HEADER_LINES = 4

@enum DataFileType seperatefiles sequentialtxt rowwisetxt json

function getfirstsaveddate(filename::String)
    file = open(filename, "r")
    firstline = readline(file)
    dateline = readline(file)
    saveddate::DateTime = DateTime(0)
        
    # if the current line has a date, read it
    try 
        saveddate = DateTime(dateline[DATE_START_INDEX:end])
    catch e
        saveddate = DateTime(0)
    end
    close(file)

    return saveddate
end

function checkfilename(outputfilename::String)::String
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

"""
Checks if data directory exists for the specified date. If not, makes the directory. Returns the dir name

By default checks for the current date.
"""
function getdatedir(date::Date=today(), currdir::String=pwd())
    while basename(currdir) != "work"
        currdir = dirname(currdir)
    end

    newdirname = "$(currdir)/data/$(Dates.day(date))-$(Dates.month(date))"
    if !isdir(newdirname)
        # if dir doesn't exist, make it and return the name
        mkdir(newdirname)
    end

    return newdirname

end

function getsimsegments(inputfilename::String)
    # initialize parameters and matrices        
    simdates::Array{DateTime} = Array{DateTime}[]
    simtimes::Array{Float64} = Array{Float64}[]
    simsaves::Array{Int64} = Array{Int64}[]
    simstarts::Array{Int64} = Array{Int64}[]
    linecounter::Int64 = 0

    # read file
    file = open(inputfilename, "r")

    line = readline(file)
    linecounter += 1
    # Check that data file header is correct
    @assert (contains(lowercase(line), "row wise txt") || contains(lowercase(line), "rowwise txt")) "Rowwisetxt data file does not have correct first line. File may be corrupted."
    
    line = readline(file)
    linecounter += 1

    # each instance of the while loop begins on the timestamp
    while !eof(file)
        # if the current line is a date, read it. Else skip to next line
        try 
            push!(simdates, DateTime(line[DATE_START_INDEX:end]))
        catch e
            push!(simdates, DateTime(0))
        else
            line = readline(file)
            linecounter += 1
        end

        @assert contains(lowercase(line), "simulation parameters") "Rowwisetxt data file does not have header. File may be corrupted."
    
        # read next line to get number of position data points
        simparaminfo_str = readline(file)
        linecounter += 1
        simparams = SimulationParameters(simparaminfo_str)
        
        # assert start times are correct
        @assert simparams.starttime == sum(simtimes)
        
        push!(simsaves, getnsaves(simparams))
        push!(simtimes, simparams.totaltime)
        
        # now skip the particle state data
        line = readline(file)
        linecounter += 1
        @assert contains(lowercase(line), "particle states") "Rowwisetxt data file does not have correct `Particle States` header. File may be corrupted."
        
        push!(simstarts, linecounter+1)

        linestoskip = getnsaves(simparams)
        for _ in 1:linestoskip
            readline(file)
            linecounter += 1
        end

        line = readline(file)
        linecounter += 1
    end    
    # close file
    close(file)

    simdata = SimulationData(simparams, getsavetimes(simparams), posmatrix, spinsmatrix)
    return simdata
end


"""
Gets the simdata of the nth segment of a data file.

`segment_n` is indexed at 1
"""
function loadnthsimsegment(inputfilename::String, segment_n::Int, startlines::Array{Int}=nothing)
    # initialize parameters and matrices        
    if isnothing(startlines)
        startlines = getsimsegments(inputfilename)[1]
    end

    startline::Int64 = startlines[segment_n]


    # read file
    file = open(inputfilename, "r")

    # skip to first particle state line in the desired segment
    linestoskip = startline - 3
        for _ in 1:linestoskip
            readline(file)
        end

    # read next line to get simulation params
    simparaminfo_str = readline(file)
    simparams = SimulationParameters(simparaminfo_str)
    ntimes = getnsaves(simparams)

    line = readline(file)
    # Check that data file header is correct
    @assert contains(lowercase(line), "particle states") "Rowwisetxt data file does not have correct `Particle States` header. File may be corrupted."
    
    # positions and spins
    # read in data from file
    posmatrix = zeros(ntimes, simparams.numparticles)
    spinsmatrix = zeros(ntimes, simparams.numparticles)
    @showprogress 1 "Loading Simulation Segment $(segment_n) (N=$(simparams.numparticles), nsteps=$(ntimes))..." for i in 1:ntimes
        dataline = readline(file)
        particle_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
        posmatrix[i,:] = copy(particle_i_data[1:simparams.numparticles])
        spinsmatrix[i,:] = Int8.(copy(particle_i_data[simparams.numparticles+1:end]))
    end
    
    # close file
    close(file)
    return simstarts, simsaves, simtimes, simdates, simstarts
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
        try 
            push!(simdates, DateTime(line[DATE_START_INDEX:end]))
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
        @showprogress 1 "Loading Simulation (N=$(simparams.numparticles), nsteps=$(ntimes))..." for i in 1:ntimes
        # for i in 1:ntimes                    
            dataline = readline(file)
            particle_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
            posmatrix[i,:] = copy(particle_i_data[1:simparams.numparticles])
            spinsmatrix[i,:] = Int8.(copy(particle_i_data[simparams.numparticles+1:end]))
        end

        # close file
        close(file)
            
        # store matrix data as `SimulationData` object and return
        simdata = SimulationData(simparams, getsavetimes(simparams), posmatrix, spinsmatrix)
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
    try 
        push!(simdates, DateTime(line[DATE_START_INDEX:end]))
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
    simdata = SimulationData(simparams, getsavetimes(simparams), posmatrix, spinsmatrix)
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
        particlestates::Matrix{Any} = zeros(getnsaves(simdata.simparams), 2 * simdata.simparams.numparticles)
        particlestates[:,1:simdata.simparams.numparticles] = simdata.positions
        particlestates[:,simdata.simparams.numparticles+1:end] = simdata.spins
        writedlm(io, particlestates, ",")
    end
end



"""
Saves simulation data from a `SimulationData` object to a data file.

Can specify the output file type with `filetype`.

By default, asks for confirmation if overwriting an existing file.
Set `forceclear` to `true` to clear the file without warning.

"""
function savesim_compressed(simdata::SimulationData, outputfilename::String, filetype::DataFileType=rowwisetxt, forceclear::Bool=false)         
    # Ask to confirm if overwriting existing file
    if !forceclear
        filenametouse = checkfilename(outputfilename)
    else
        filenametouse = outputfilename
    end
    
    @assert filetype == rowwisetxt "Can only compress `roweisetxt` files."
    
    open(filenametouse, "w") do io
        println(io, "Row Wise txt")  # Row Wise .txt
        println(io, "Saved at $(round(now(), Dates.Second(1)))")  # Timestamp
        println(io, "Simulation Parameters")
        println(io, csv_serialize(simdata.simparams))
        println(io, "Particle States ([positions], [spins])")
        # construct relative particle state matrix by concatenating states of position and spin
        particlestates::Matrix{Any} = zeros(getnsaves(simdata.simparams), 2 * simdata.simparams.numparticles)
        particlestates[:,1:simdata.simparams.numparticles] = simdata.positions
        particlestates[:,simdata.simparams.numparticles+1:end] = simdata.spins
        writedlm(io, particlestates, ",")
    end

end
