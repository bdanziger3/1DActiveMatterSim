include("../simulation/basic_sim.jl")
# include("./sim_structs.jl")

using JSON

@enum DataFileType seperatefiles sequentialtxt rowwisetxt json


dft::DataFileType = seperatefiles

function loadsim(inputfilename::String, filetype::DataFileType)::SimulationData
    # initialize parameters and matrices        
    ntimes::Int64 = 0
    posmatrix::Matrix{Float64} = zeros(0, 0)
    spinsmatrix::Matrix{Int8} = zeros(0, 0)

    if filetype == rowwisetxt
        # read file
        file = open(inputfilename, "r")

        line = readline(file)
        # Check that data file header is correct
        @assert (contains(lowercase(line), "row wise txt") || contains(lowercase(line), "rowwise txt")) "Rowwisetxt data file does not have correct first line. File may be corrupted."
        
        line = readline(file)
        @assert contains(lowercase(line), "simulation parameters") "Rowwisetxt data file does not have header. File may be corrupted."
       
        # read next line to get number of position data points
        simparaminfo_str = readline(file)
        simparams = SimulationParameters(simparaminfo_str)
        ntimes = getntimes(simparams)
            
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
        ntimes = getntimes(simparams)
            
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

    # read file
    file = open(inputfilename, "r")

    line = readline(file)
    # Check that data file header is correct
    @assert (contains(lowercase(line), "row wise txt") || contains(lowercase(line), "rowwise txt")) "Rowwisetxt data file does not have correct first line. File may be corrupted."
    
    line = readline(file)
    @assert contains(lowercase(line), "simulation parameters") "Rowwisetxt data file does not have header. File may be corrupted."
    
    # read next line to get number of position data points
    simparaminfo_str = readline(file)
    simparams = SimulationParameters(simparaminfo_str)
    ntimes = getntimes(simparams)
        
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

By default appends the data to the file at `outputfilename` following the structure of the `rowwisetxt` DataFileType.
Set `clearfile` to `true` to clear the file before writing instead.

"""
function savesim(simdata::SimulationData, outputfilename::String, filetype::DataFileType=rowwisetxt, clearfile::Bool=false)         
    if filetype == sequentialtxt
        filestring = "$(outputfilename)_simdata.txt"
        open(filestring, clearfile ? "w" : "a") do io
            println(io, "Sequential txt")   # Sequential .txt
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simdata.simparams))
            println(io, "Positions")
            writedlm(io, simdata.positions, ",")
            println(io, "Spins")
            writedlm(io, simdata.spins, ",")
        end
    elseif filetype == rowwisetxt
        filestring = "$(outputfilename)_simdata.txt"
        open(filestring, clearfile ? "w" : "a") do io
            println(io, "Row Wise txt")  # Row Wise .txt
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simdata.simparams))
            println(io, "Particle States ([positions], [spins])")
            # construct particle state matrix by concatenating states of position and spin
            particlestates::Matrix{Any} = zeros(getntimes(simdata.simparams), 2 * simdata.simparams.numparticles)
            particlestates[:,1:simdata.simparams.numparticles] = simdata.positions
            particlestates[:,simdata.simparams.numparticles+1:end] = simdata.spins
            writedlm(io, particlestates, ",")
        end
    elseif filetype == seperatefiles
        simdatafile_prefix = outputfilename
        writedlm("$(simdatafile_prefix)_x.txt", simdata.positions, ",")
        writedlm("$(simdatafile_prefix)_xwrap.txt", simdata.wrappedpositions, ",")
        writedlm("$(simdatafile_prefix)_S.txt", simdata.spins, ",")
        writedlm("$(simdatafile_prefix)_t.txt", simdata.times, ",")
    end

end





# sd = loadsim("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N3-strong-alignsimple_simdata.txt", sequentialtxt)

# savesim(sd, "saveas_rwt.txt", rowwisetxt, true)

# sdr = loadsim("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/saveas_rwt.txt_simdata.txt", rowwisetxt)

# println(sd == sdr)

# savesim(sdr, "saveas_rwt2.txt", rowwisetxt, true)
