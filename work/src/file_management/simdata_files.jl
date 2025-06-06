include("../simulation/basic_sim.jl")
# include("./sim_structs.jl")

using JSON

@enum DataFileType seperatefiles sequentialtxt json


dft::DataFileType = seperatefiles

function loadsim(inputfilename::String, filetype::DataFileType)::SimulationData
    if filetype == sequentialtxt
        # read file
        file = open(inputfilename, "r")

        line = readline(file)
        # Check that data file header is correct
        @assert contains(lowercase(line), "simulation parameters") "Sequentialtxt data file does not have correct first line. File may be corrupted."
       
        # read next line to get number of position data points
        simparaminfo_str = readline(file)
        simparams = SimulationParameters(simparaminfo_str)
        ntimes::Int64 = getntimes(simparams)
            
        # now read through the data in the order: [positions, wrappedpositions, spins]
        line = readline(file)
        @assert contains(lowercase(line), "positions") "Sequentialtxt data file does not have correct `Positions` header. File may be corrupted."
        
        # positions
        # initialize matrix and read in data from file
        posmatrix::Matrix{Float64} = zeros(ntimes, simparams.numparticles)
        for i in 1:ntimes                    
            dataline = readline(file)
            pos_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
            posmatrix[i,:] = copy(pos_i_data)
        end
        
        # wrappedpositions
        # initialize matrix and read in data from file
        line = readline(file)
        @assert contains(lowercase(line), "wrapped positions") "Sequentialtxt data file does not have correct `Wrapped Positions` header line. File may be corrupted."
        
        wrappedposmatrix::Matrix{Float64} = zeros(ntimes, simparams.numparticles)
        for i in 1:ntimes                    
            dataline = readline(file)
            wrappedpos_i_data::Array{Float64} = parse.(Float64, split(dataline, ","))
            wrappedposmatrix[i,:] = copy(wrappedpos_i_data)
        end
        
        # spins
        # initialize matrix and read in data from file
        line = readline(file)
        @assert contains(lowercase(line), "spins") "Sequentialtxt data file does not have correct `Spins` header line. File may be corrupted."
        
        spinsmatrix::Matrix{Int8} = zeros(ntimes, simparams.numparticles)
        for i in 1:ntimes                    
            dataline = readline(file)
            spins_i_data::Array{Int8} = parse.(Int8, split(dataline, ","))
            spinsmatrix[i,:] = copy(spins_i_data)
        end

        # close file
        close(file)
            
        # store matrix data as `SimulationData` object and return
        simdata = SimulationData(simparams, gettimes(simparams), posmatrix, wrappedposmatrix, spinsmatrix)
        return simdata
    end
end


"""
Saves simulation data from a `SimulationData` object to a data file.

Can specify the output file type with `filetype`.

By default appends the data to the file at `outputfilename`.
Set `clearfile` to `true` to clear the file before writing instead.

"""
function savesim(simdata::SimulationData, outputfilename::String, filetype::DataFileType, clearfile::Bool=false)         
    if filetype == sequentialtxt
        filestring = "$(outputfilename)_simdata.txt"
        open(filestring, clearfile ? "w" : "a") do io
            println(io, "Simulation Parameters")
            println(io, csv_serialize(simdata.simparams))
            println(io, "Positions")
            writedlm(io, simdata.positions, ",")
            println(io, "Spins")
            writedlm(io, simdata.spins, ",")
        end
    elseif filetype == seperatefiles
        simdatafile_prefix = outputfilename
        writedlm("$(simdatafile_prefix)_x.txt", simdata.positions, ",")
        writedlm("$(simdatafile_prefix)_xwrap.txt", simdata.wrappedpositions, ",")
        writedlm("$(simdatafile_prefix)_S.txt", simdata.spins, ",")
        writedlm("$(simdatafile_prefix)_t.txt", simdata.times, ",")
    end

end





# loadsim("testfunc_simdata.txt", sequentialtxt)

