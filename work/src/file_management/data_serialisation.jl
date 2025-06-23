include("../simulation/basic_sim.jl")
# include("../sim÷ulation/sim_structs.jl")

using StringEncodings

using JSON
using Dates

@enum Serialisation utf8 ascii128 none

UTF8BITS::Int8 = 20
ASCIIBITS::Int8 = 6

ASCIICHARS = UInt8(94)
ASCIIMININT = 40

BITSDICT::Dict{Serialisation, Int8} = Dict(utf8 => UTF8BITS, ascii128 => ASCIIBITS)

"""
Calculate the total number of encoding characters each of width `encodingbits` needed to encode the data that is of size `databits`.

For spin data, `databits` should be the number of particles

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
Serialises a list of spins into a compressed string.

Serialisation schemes:
`utf8`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of utf-8 characters

`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit ascii characters.
However, uses only characters in the range of 40-126 to prevent control, newline, whitespace, and special characters
"""
function serialisespins(spins::Array{<:Number}, serialisation::Serialisation=utf8)::String
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
    if serialisation == utf8
        # encoding = Char.(digits(bitsum, base=UInt32(2^UTF8BITS)))
        # encoding = Char.(intencoding)
        # return join(encoding)
        return 1
    elseif serialisation == ascii128
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
Deserialises a string into a list of spin values.

Deserialises into a UInt and then converts (0 -> 1, 1 -> -1)

Serialisation schemes:
`utf8`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of utf-8 characters
`ascii128`: Converts the set of spin values into a binary number (1 -> 0, -1 -> 1) and converts this to a string of 7-bit-8 characters
"""
function deserialisespins(datastr::String, numparticles::Int, serialisation::Serialisation=utf8)::Array{Int8}
    spins::Array{Int8} = zeros(1, numparticles)
    
    if serialisation == ascii128
        t1 = time()
        # convert ASCII string into array of Int8s
        ints = Int8.(collect(datastr)) .- ASCIIMININT
        t2 = time()
        # convert into one binary row array
        binspins = transpose(vcat(digits.(ints, base=Int8(2), pad=ASCIIBITS)...))
        t3 = time()
        # convert from binary to spins (0 -> 1, 1 -> -1)
        spins = (s -> s==1 ? Int8(-1) : Int8(1)).(binspins)
        t4 = time()
        # remove padding on the end by only returning `numparticles` spins
        println("$(t2-t1),$(t3-t2),$(t4-t3),$(t4)")
        return spins[1:numparticles]
    else
        # interpret is as a normal string with 1s and -1s separated by commas
        return transpose(parse.(Int8, split(datastr, ",")))
    end
end

