include("../src/file_management/simdata_files.jl")
include("../src/file_management/data_serialisation.jl")
include("test_fixtures.jl")


"""
Serializes and deserializes randomly generated spins and confirms that the values are the same before and after.
"""
function test_spinserialisation()
    N::Int64 = 1e4
    # Generate random spins and serialize and deserialize
    spins = rand([1,-1], N)
    # println(spins)
    
    t0 = time()
    serializedstr_ascii = serializespins(spins, ascii128)
    t1 = time()
    # println(serializedstr_ascii)

    td0 = time()
    newspins_ascii = deserializespins(serializedstr_ascii, N, ascii128)
    td1 = time()

    @assert newspins_ascii == spins "ASCII deserialized spins are not the same as the initial spins"

    println("Spin Serialization Timing results:\nASCII: encoding: $(t1 - t0) seconds, decoding: $(td1 - td0) seconds")
end

"""
Serializes and deserializes simulated positions and confirms that the values are the same before and after.
"""
function test_positionserialisation()
    simparams = simparams_10


    # run simple simulation if data file does not exist
    datafile_path = "test_data/serialisation_test_data.txt"
    if !isfile(datafile_path)
        simdata = runsim(simparams)
        savesim(simdata, datafile_path, rowwisetxt, true)
    end

    # load data and try to serialize the spins
    loadedsim::SimulationData = loadsim(datafile_path, rowwisetxt)

    # println(loadedsim.positions)
    nsaves = getnsaves(simparams)

    str_array::Array{String} = Array{String}(undef, 1, nsaves - 1)

    encodingtimes::Array{Float64} = zeros(1, nsaves-1)
    for i in 2:nsaves
        t0 = time()
        str_array[i-1] = serializepositions(loadedsim.positions[i, :], loadedsim.positions[1, :], loadedsim.simparams)
        t1 = time()

        encodingtimes[i-1] = t1 - t0
    end

    # now decode the string
    decodingtimes::Array{Float64} = zeros(1, nsaves-1)
    
    positions_deserialized::Matrix{Float64} = zeros(nsaves, simparams.numparticles)
    positions_deserialized[1, :] = loadedsim.positions[1, :]
    for i in 2:nsaves
        t0 = time()
        positions_deserialized[i, :] = deserializepositions(str_array[i-1], loadedsim.positions[1, :], loadedsim.simparams)
        t1 = time()

        decodingtimes[i-1] = t1 - t0
    end

    # check that deserialized positions are the same as the originals
    threshold = simparams.v0 * simparams.dt * .001  # within 0.1% of dx 
    @assert all(abs.(positions_deserialized .- loadedsim.positions) .<= threshold) "Deserialized positions do not equal the original positions."
    
    println("Position Serialization Timing results:\nEncoding: mean: $(sum(encodingtimes) / length(encodingtimes)) seconds, max: $(maximum(encodingtimes)) seconds.\nDecoding: mean: $(sum(decodingtimes) / length(decodingtimes)) seconds, max: $(maximum(decodingtimes))")
end

    



test_spinserialisation()
test_positionserialisation()

println("PASSED")