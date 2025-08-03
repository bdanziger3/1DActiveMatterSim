include("../simulation/basic_sim.jl")
# include("../sim÷ulation/sim_structs.jl")

using JSON
using Dates

const ASCIIBITS::Int64 = 6

const ASCIIMININT = 40

const POS_SPINS_SEPARATOR = " "

"""
Calculate the total number of encoding characters each of width `encodingbits` needed to encode the data that is of size `databits`.

For spin data, `databits` should be the number of particles
    The returned value is the total number of characters needed to encode ALL spins

For position data, `databits` should be the number of bits per particle position calculated with `maxposbits`.
    The returned value is the number of characters needed PER position

For ascii128 serialisation, `encodingbits` should be 6
"""
function neededbits(databits::Int, encodingbits::Int64=ASCIIBITS)::Int
    nbits = div(databits, encodingbits)

    if databits % encodingbits != 0
        nbits += 1
    end

    return nbits
end

"""
Calculates the maximum number of bits needed to encode 1 particle's relative position data
"""
function bitsperpos(simparams::SimulationParameters)::Int64
    # max relative position is maximum number of steps that fit in the box
    dx::Float64 = simparams.v0 * simparams.dt
    maxrelativepos::Int64 = Int64(ceil(simparams.boxwidth / dx))

    # calculate how many bits needed to encode the data
    # maxrelativepos + 1 possible relative positions (+ or - direction plus 0)
    maxposbits::Int64 = Int64(ceil(log2(maxrelativepos + 1)))

    return maxposbits
end

"""
Gets the number of characters per position value based on the Simulation Parameters
"""
function charsperpos(simparams::SimulationParameters)::Int64
    return neededbits(bitsperpos(simparams), ASCIIBITS)
end

"""
Serializes a list of spins into a compressed string.

Serialisation schemes:

`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit ascii characters.
However, uses only characters in the range of 40-126 to prevent control, newline, whitespace, and special characters
"""
function serializespins(spins::Array{<:Number})::String
    # convert spins to bits (1 -> 0 ; -1 -> 1)
    bits::BitArray = (s -> s == -1).(spins)
    numparticles = length(spins)

    # there are 128 ASCII characters, but only 94 of them [33, 126]
    # are writable characters that are not controls or whitespace.
    # Skip 33-39 because they might be special characters: use 6-bit starting at 40: [40, 103]
    encodingints::Array{UInt8} = zeros(1, neededbits(numparticles, ASCIIBITS))
    for (i, b) in pairs(bits)
        slot = div(i-1, ASCIIBITS) + 1
        # add 2^(n) to the appropriate index of the array
        encodingints[slot] += (b * 2^((mod(i-1, ASCIIBITS))))
    end
    
    encodingints .+= ASCIIMININT        
    encoding = join(Char.(encodingints))
    
    return encoding
end


"""
Deserializes a string into a list of spin values.

Deserializes into a UInt and then converts (0 -> 1, 1 -> -1)

Serialisation schemes:
`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1)
            and converts this to a string of 6-bit integers represented by ascii characters
            
            However, uses only characters in the range of 40-126 to prevent control, newline, whitespace, and special characters
"""
function deserializespins(datastr::String, numparticles::Int)::Array{Int8}
    spins::Array{Int8} = zeros(1, numparticles)

    # convert ASCII string into array of Int8s
    ints = Int8.(collect(datastr)) .- ASCIIMININT
    # convert into one binary row array
    binspins = transpose(vcat(digits.(ints, base=Int8(2), pad=ASCIIBITS)...))
    # convert from binary to spins (0 -> 1, 1 -> -1)
    spins = (s -> s==1 ? Int8(-1) : Int8(1)).(binspins)
    # remove padding on the end by only returning `numparticles` spins
    return spins[1:numparticles]
end


"""
Serializes a list of positions into a compressed string encoding distance relative to initial positions.

Serialisation schemes:
`ascii128`: Converts the set of positions (in units of initeger steps of size `dx` [`dx = dt * v_0`] from the particle's initial position) values into a 6-bit signed integer
            then converts this to a string of 6-bit integers represented ascii characters.

            However, uses only characters in the range of 40-126 to prevent control, newline, whitespace, and special characters
"""
function serializepositions(positions::Array{<:Number}, initialpositions::Array{<:Number}, simparams::SimulationParameters)::String
    # represent positions as integer number of steps from original Positions
    dx::Float64 = simparams.v0 * simparams.dt
    unsignedrelativepositions::Array{<:Int} = Int64.(round.(mod.(positions - initialpositions, simparams.boxwidth) / dx))

    # calculate how many 6-bit ASCII characters are needed per position
    charspp::Int64 = charsperpos(simparams)
    
    # We are encoding as unsigned ints, so use `base=2^(bits)`
    unsignedasciiintmatrix = mapreduce(permutedims, vcat, digits.(unsignedrelativepositions, base=UInt8(2^(ASCIIBITS)), pad=charspp))
    # convert to unsigned ints, shift by `ASCIIMININT`, and convert to chars
    charmatrix = (x -> Char(x + ASCIIMININT)).(unsignedasciiintmatrix)
    stringencoding::String = join(permutedims(charmatrix))

    return stringencoding
end

"""
Deserializes a string into a list of position values.

Serialisation schemes:
`ascii128`: Converts the string back into 6-bit unsigned integers, then converts to signed integers,
            and lastly uses the initial particle positions and the size of `dx` to get the original positions.
"""
function deserializepositions(datastr::String, initialpositions::Array{<:Number}, simparams::SimulationParameters)::Array{Float64}
    # calculate encoding parameters
    dx::Float64 = simparams.v0 * simparams.dt
    charspp::Int64 = charsperpos(simparams) # how many 6-bit ASCII characters are needed per position
    
    # initialize results array
    positions::Array{Float64} = zeros(1, simparams.numparticles)
    # reshape and convert characters to unsigned ints
    unsignedintmatrix = permutedims(reshape((Int8.(collect(datastr)) .- ASCIIMININT), charspp, :))

    # sum up signed ints (scaled by digit place for a 6-bit signed int) to get relative positions
    sumcol::Array{Int64} = zeros(simparams.numparticles, 1)
    
    for i in 1:charspp
        sumcol += Int64.(unsignedintmatrix[:,i]) * (2^((ASCIIBITS) * (i-1)))
    end
    
    # convert relative positions back to absolute positions
    positions = (sumcol * dx) .+ initialpositions

    return wrap(positions, simparams.boxwidth)
end

"""
Calculate the number of decimal places needed that preserves the position information.

Keep `EXCESS_DECIMAL_PLACES=2` decimal places smaller than `dx`
"""
function roundpositions(positions::Array{Float64}, simparams::SimulationParameters, returnarray::Bool=false)::Any
    EXCESS_DECIMAL_PLACES::Int64 = 2
    dx::Float64 = simparams.dt * simparams.v0
    ndecimalplaces::Int64 = Int64(ceil(-log10(dx))) + EXCESS_DECIMAL_PLACES
    roundedpositions = round.(positions, digits=ndecimalplaces)

    if returnarray
        return roundedpositions
    else
        return join(roundedpositions, ",")
    end
end


"""
Take a string of regular data from a file and return the serialized version
"""
function serializedatafileline(dataline::String, initialpositions::Array{Float64}, simparams::SimulationParameters, isfirstline::Bool=false)::String
    # load in data, compress, save outstring
    particle_data = parse.(Float64, split(dataline, ","))
    positions = copy(particle_data[1:simparams.numparticles])
    spins = Int8.(copy(particle_data[simparams.numparticles+1:end]))

    positionstring = serializepositions(positions, initialpositions, simparams)
    spinsstring = serializespins(spins)

    compressedline::String = ""

    if isfirstline
        roundedpositionstring = roundpositions(positions, simparams)
        compressedline = "$(roundedpositionstring)$(POS_SPINS_SEPARATOR)$(spinsstring)"
    else
        compressedline = "$(positionstring)$(POS_SPINS_SEPARATOR)$(spinsstring)"
    end

    return compressedline
end

"""
Take a string of serialized data from a file and return the deserialized version as a string.
"""
function deserializedatafileline(dataline::String, initialpositions::Array{Float64}, simparams::SimulationParameters, isfirstline::Bool=false)::String
    # load in string, separate positions and spins
    datalinelist = split(dataline, POS_SPINS_SEPARATOR)
    posstring::String = datalinelist[1]
    encodedspins::String = datalinelist[2]

    spindata = deserializespins(encodedspins, simparams.numparticles)

    if isfirstline
        # For the first state, keep the position data and use the deserialized spins
        return "$(posstring),$(join(spindata, ','))"
    else
        # decode the positions as well based on the initial positions
        positiondata::Array{Float64} = deserializepositions(posstring, initialpositions, simparams)
        roundedpositionstring = roundpositions(positiondata, simparams)
        return "$(roundedpositionstring),$(join(spindata, ','))"
    end
end


"""
Take a string of serialized data from a file and return either the spins or position data.
"""
function deserializelineonedatatype(dataline::String, datatype::String, simparams::SimulationParameters, initialpositions::Array{Float64}=zeros(), isfirstline::Bool=false)::Array{<:Real}
    # load in string, separate positions and spins
    datalinelist = split(dataline, POS_SPINS_SEPARATOR)

    if contains(lowercase(datatype), "spin")
        # spins
        encodedspins::String = datalinelist[2]
        return deserializespins(encodedspins, simparams.numparticles)

    else
        # positions
        posstring::String = datalinelist[1]

        if isfirstline
            # For the first state, keep the position data
            return parse.(Float64, split(posstring, ","))
        else
            # decode the positions based on the initial positions
            positiondata::Array{Float64} = deserializepositions(posstring, initialpositions, simparams)
            return roundpositions(positiondata, simparams, true)
        end
    end
end

"""
Take a string of serialized data from a file and return the positions and spins as arrays.
"""
function deserializelinegetarrays(dataline::String, initialpositions, simparams::SimulationParameters, isfirstline::Bool=false)::Tuple{Array{<:Real}, Array{<:Real}}
    # load in string, separate positions and spins
    datalinelist = split(dataline, POS_SPINS_SEPARATOR)

    posstring::String = datalinelist[1]
    encodedspins::String = datalinelist[2]
    spins = deserializespins(encodedspins, simparams.numparticles)

    if isfirstline
        # For the first state, keep the position data
        return parse.(Float64, split(posstring, ",")), spins
    else
        # decode the positions based on the initial positions
        positiondata::Array{Float64} = deserializepositions(posstring, initialpositions, simparams)
        return roundpositions(positiondata, simparams, true), spins
    end
end

