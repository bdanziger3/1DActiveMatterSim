include("../simulation/basic_sim.jl")
# include("../sim÷ulation/sim_structs.jl")

using JSON
using Dates

@enum Serialisation utf8 ascii128 none

const UTF8BITS::Int8 = 20
const ASCIIBITS::Int8 = 6

const ASCIIMININT = 40

const BITSDICT::Dict{Serialisation, Int8} = Dict(utf8 => UTF8BITS, ascii128 => ASCIIBITS)

const POS_SPINS_SEPARATOR = " "

"""
Calculate the total number of encoding characters each of width `encodingbits` needed to encode the data that is of size `databits`.

For spin data, `databits` should be the number of particles
    The returned value is the total number of characters needed to encode ALL spins

For position data, `databits` should be the number of bits per particle position calculated with `maxposbits`.
    The returned value is the number of characters needed PER position

For utf8 serialisation, `encodingbits` should be 7
For ascii128 serialisation, `encodingbits` should be 20
"""
function neededbits(databits::Int, encodingbits::Union{Int, Serialisation})::Int
    if typeof(encodingbits) == Serialisation
        encodingbits = BITSDICT[encodingbits]
    end

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
    # (2 * maxrelativepos) + 1 possible relative positions (+ or - direction plus 0)
    maxposbits::Int64 = Int64(ceil(log2((2 * maxrelativepos) + 1)))

    return maxposbits
end

"""
Gets the number of characters per position value based on the Simulation Parameters
"""
function charsperpos(simparams::SimulationParameters, serialisation::Serialisation=ascii128)::Int64
    return neededbits(bitsperpos(simparams), serialisation)
end

"""
Serializes a list of spins into a compressed string.

Serialisation schemes:
`utf8`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of utf-8 characters

`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit ascii characters.
However, uses only characters in the range of 40-126 to prevent control, newline, whitespace, and special characters
"""
function serializespins(spins::Array{<:Number}, serialisation::Serialisation=utf8)::String
    # convert spins to bits (1 -> 0 ; -1 -> 1)
    bits::BitArray = (s -> s == -1).(spins)
    numparticles = length(spins)
    # builds an array of appropriately-sized integers based on the bits
    # 20-bit integers for `utf8` and 7-bit integers for `ascii128`
    # intencoding::Array{UInt32} = zeros(1, neededbits(length(spins), serialisation))
    # wordwidth = BITSDICT[serialisation]
    # for (i, b) in pairs(bits)
    #     # add 2^(n) to the appropriate index of the array
    #     arrayindex = div(i-1, wordwidth) + 1
    #     bitexp = (i-1) % wordwidth
    #     intencoding[arrayindex] += b * 2^(bitexp)
    # end
    if serialisation == ascii128
        # there are 128 ASCII characters, but only 94 of them [33, 126]
        # are writable characters that are not controls or whitespace.
        # Skip 33-39 because they might be special characters: use 6-bit starting at 40: [40, 103]
        encodingints::Array{UInt8} = zeros(1, neededbits(numparticles, ascii128))
        for (i, b) in pairs(bits)
            slot = div(i-1, ASCIIBITS) + 1
            # add 2^(n) to the appropriate index of the array
            encodingints[slot] += (b * 2^((mod(i-1, ASCIIBITS))))
        end
        
        encodingints .+= ASCIIMININT        
        encoding = join(Char.(encodingints))
        
        return encoding
    else
        return join(spins, ",")
    end
end


"""
Deserializes a string into a list of spin values.

Deserializes into a UInt and then converts (0 -> 1, 1 -> -1)

Serialisation schemes:
`utf8`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of utf-8 characters
`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit-8 characters
"""
function deserializespins(datastr::String, numparticles::Int, serialisation::Serialisation=utf8)::Array{Int8}
    spins::Array{Int8} = zeros(1, numparticles)
    
    if serialisation == ascii128
        # convert ASCII string into array of Int8s
        ints = Int8.(collect(datastr)) .- ASCIIMININT
        # convert into one binary row array
        binspins = transpose(vcat(digits.(ints, base=Int8(2), pad=ASCIIBITS)...))
        # convert from binary to spins (0 -> 1, 1 -> -1)
        spins = (s -> s==1 ? Int8(-1) : Int8(1)).(binspins)
        # remove padding on the end by only returning `numparticles` spins
        return spins[1:numparticles]
    else
        # interpret is as a normal string with 1s and -1s separated by commas
        return transpose(parse.(Int8, split(datastr, ",")))
    end
end



"""
Serializes a list of positions into a compressed string encoding distance relative to initial positions.

Serialisation schemes:
`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit ascii characters.
However, uses only characters in the range of 40-126 to prevent control, newline, whitespace, and special characters
"""
function serializepositions(positions::Array{<:Number}, initialpositions::Array{<:Number}, simparams::SimulationParameters)::String
    # represent positions as integer number of steps from original Positions
    dx::Float64 = simparams.v0 * simparams.dt
    relativepositions::Array{<:Int} = Int64.(round.((positions - initialpositions) / dx))

    # calculate how many 6-bit ASCII characters are needed per position
    charspp::Int64 = charsperpos(simparams)
    
    # We are encoding as signed ints, so use `base=2^(bits-1)`
    signedasciiintmatrix = mapreduce(permutedims, vcat, digits.(relativepositions, base=Int8(2^(ASCIIBITS-1)), pad=charspp))

    # convert to unsigned ints, shift by `ASCIIMININT`, and convert to chars
    charmatrix = (x -> x < 0 ? Char(x + 2^ASCIIBITS + ASCIIMININT) : Char(x + ASCIIMININT)).(signedasciiintmatrix)
    stringencoding::String = join(permutedims(charmatrix))

    return stringencoding
end

"""
Deserializes a string into a list of position values.

Deserializes into a UInt and then converts (0 -> 1, 1 -> -1)

Serialisation schemes:
`utf8`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of utf-8 characters
`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit-8 characters
"""
function deserializepositions(datastr::String, initialpositions::Array{<:Number}, simparams::SimulationParameters)::Array{Float64}
    # calculate encoding parameters
    dx::Float64 = simparams.v0 * simparams.dt
    charspp::Int64 = charsperpos(simparams) # how many 6-bit ASCII characters are needed per position
    
    # initialize results array
    positions::Array{Float64} = zeros(1, simparams.numparticles)
    # reshape and convert characters to unsigned ints
    unsignedintmatrix = permutedims(reshape((Int8.(collect(datastr)) .- ASCIIMININT), charspp, :))
    # convert back into signed ints
    signedintmatrix = (x -> x > 2^(ASCIIBITS-1) ? x - 2^ASCIIBITS : x).(unsignedintmatrix)

    # sum up signed ints (scaled by digit place for a 6-bit signed int) to get relative positions
    sumcol::Array{Int64} = zeros(simparams.numparticles, 1)
    
    for i in 1:charspp
        sumcol += Int64.(signedintmatrix[:,i]) * (2^((ASCIIBITS - 1) * (i-1)))
    end
    
    # convert relative positions back to absolute positions
    positions = (sumcol * dx) .+ initialpositions

    return positions
end