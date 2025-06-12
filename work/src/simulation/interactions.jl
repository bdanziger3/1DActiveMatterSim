# Library of interaction functions between particles to use in the simulation
# include("./sim_structs.jl")




function randlinflip(dt::Real, fliprate::Real)
    return rand() < (dt * fliprate)
end

function findnn(currpositions::Matrix{<:Real}, particlei::Int)
    # return indices of left and right neighbors for particle i
    
    nparticles = length(currpositions)
    x_i = currpositions[particlei]
    left = -1
    leftdx = Inf
    right = -1
    rightdx = Inf
    ontop = []
    for p_i in 1:nparticles
        if p_i != particlei

            d_p = currpositions[p_i] - x_i
            if d_p < 0
                if abs(d_p) < leftdx
                    leftdx = abs(d_p)
                    left = p_i
                end
            elseif  d_p > 0
                if abs(d_p) < rightdx
                    rightdx = abs(d_p)
                    right = p_i
                end
            else
                push!(ontop, p_i)
            end
        end
    end

    # if no overlaps
    if length(ontop) == 0
        return [left, right]
    elseif length(ontop) == 1
        if right == -1 || leftdx <= rightdx
            return [left, ontop[1]]
        else
            return [ontop[1], right]
        end
    else
        # if at least 2 overlapping particles, return them all 
        return ontop
    end
end


function alignOn2(simparams::SimulationParameters, currpositions::Matrix{<:Real}, currspins::Matrix{Int8}, antialign::Bool = false)
    nparticles =  length(currpositions)
    canflip::Array{Bool} = Array{Bool}(undef, 1, nparticles)
    canflip .= false

    # naive approach: check all particles to find nearest neighbors    
    for p_i in 1:nparticles
        # find nearest neighbors of each particle
        neighbors = findnn(currpositions, p_i)
        
        # calculate total spin of all neighboring particles
        spinsum = 0
        for p in neighbors
            if p != -1
                spinsum += currspins[p]
            end
        end
        polarity = sign(spinsum)

        # set p_i as elligible for an interaction flip if the spinsum is nonzero and opposite from the current spin of p_i
        if currspins[p_i] == -polarity
            canflip[p_i] = true
        end
    end

    # now for all particles that can flip, check if they do flip
    # flips = (x -> randomflips(simparams.dt, simparams.interactionfliprate)).(canflip)

    flips_int::Array{Bool} = Array{Bool}(undef, 1, simparams.numparticles)

    for p in 1:nparticles
        if canflip[p]
            flips_int[p] = randlinflip(simparams.dt, simparams.interactionfliprate)
        else
            flips_int[p] = false
        end
    end

    return flips_int
end 

function alignOn2_quick(simparams::SimulationParameters, currpositions::Matrix{<:Real}, currspins::Matrix{Int8}, antialign::Bool = false)
    nparticles =  length(currpositions)
    canflip::Array{Bool} = Array{Bool}(undef, 1, nparticles)
    canflip .= false

    # naive approach
    # check all particles to find nearest neighbors
    mindx = simparams.v0 * simparams.dt

    # can remember pairs so only have to check half of all pairs
    leftdict = Dict{Int, Int}()
    rightdict = Dict{Int, Int}()
    overlapdict = Dict{Int, Set{Int}}()
    overlap3set = Set{Int}()
    overlap2neighbors = Dict{Int, Tuple{Int}}()
    for p_i in 1:nparticles
        leftneighbor = -1
        leftdist = -Inf
        rightneighbor = -1
        rightdist = Inf

        # get discretized position
        positionint = Int64(round(currpositions[p_i] / mindx))

        if p_i in keys(leftdict) && p_i in keys(rightdict)

        end

        # find nieghbors if not found already
        if !((p_i in overlap3set) || (p_i in keys(leftdict)) && (p_i in keys(rightdict)))
            for p_j in 1:nparticles
                if p_i != p_j
                    d_ji = currpositions[p_j] - currpositions[p_i]
                    
                    if abs(d_ji) < (mindx / 2)
                        # if particles on top of each other, count how many at each position

                        if !(positionint in keys(overlapdict))
                            overlapdict[positionint] = Set([p_i, p_j])
                        else
                            push!(overlapdict[positionint], p_i)
                            push!(overlapdict[positionint], p_j)
                        end

                    elseif d_ji < 0 && d_ji > leftdist
                        leftdist = d_ji
                        leftneighbor = p_j
        
                    elseif d_ji > 0 && d_ji < rightdist
                        rightdist = d_ji
                        rightneighbor = p_j
                    end
                end
            end
        end

        # if overlapping with 2 or more particles, use average spins of all overlap
        if positionint in keys(overlapdict) && length(overlapdict[positionint]) >= 3
            for p in overlapdict[positionint]
                push!(overlap3set, p)
            end
            # find total of all spins of particles overlapping with p_i
            spinsum = 0
            for p in overlapdict[positionint]
                if p != p_i
                    spinsum += currspins[p]
                end
            end

            # set p_i as elligible for an interaction flip if the spinsum is nonzero and opposite from the current spin of p_i
            if abs(spinsum) > 0.5 && spinsum * currspins[p_i] < 0
                canflip[p_i] = true
            end

            continue

        elseif positionint in keys(overlapdict) && length(overlapdict[positionint]) == 2
            # if first time here, add to dictionary
            if !(positionint in keys(overlap2neighbors))
                overlap2neighbors[positionint] = Tuple{Int}([leftneighbor, rightneighbor])
            end
            
            # find other particle
            other_p = -1
            for p_k in overlapdict[positionint]
                if p_k != p_i
                    other_p = p_k
                end
            end
            # if only overlapping with 1 particle, pick a left and a right
            if (p_i in keys(leftdict)) && !(p_i in keys(rightdict))
                rightdict[p_i] = other_p
                leftdict[other_p] = p_i
            elseif !(p_i in keys(leftdict)) && (p_i in keys(rightdict))
                leftdict[p_i] = other_p
                rightdict[other_p] = p_i
            elseif !(p_i in keys(leftdict)) && !(p_i in keys(rightdict))
                # if neither assigned yet, make p_i the left ones
                rightdict[p_i] = other_p
                leftdict[p_i] = leftneighbor

                leftdict[other_p] = p_i
            else
                # if both already assigned, reassign just these
                rightdict[p_i] = other_p
                leftdict[other_p] = p_i
            end

        else
            # not in overlap, just use neighbors
            leftdict[p_i] = leftneighbor
            rightdict[leftneighbor] = p_i
            rightdict[p_i] = rightneighbor
            leftdict[rightneighbor] = p_i
        end

        # if not assigned yet, check if can flip
        if leftdict[p_i] == -1 && rightdict[p_i] == -1
            spinsum = 0
        elseif leftdict[p_i] == -1
            spinsum = currspins[rightdict[p_i]]
        elseif rightdict[p_i] == -1
            spinsum = currspins[leftdict[p_i]]
        else
            spinsum = currspins[leftdict[p_i]] + currspins[rightdict[p_i]]
        end
        # set p_i as elligible for an interaction flip if the spinsum is nonzero and opposite from the current spin of p_i
        if abs(spinsum) > 0.5 && spinsum * currspins[p_i] < 0
            canflip[p_i] = true
        end
    end

    # now for all particles that can flip, check if they do flip
    # flips = (x -> randomflips(simparams.dt, simparams.interactionfliprate)).(canflip)

    flips_int::Array{Bool} = Array{Bool}(undef, 1, simparams.numparticles)
    for p in 1:nparticles
        if flips_int[p]
            flips_int[p] = randlinflip(simparams.dt, simparams.interactionfliprate)
        end
    end

    return flips_int
end 

function calcspinflips(simparams::SimulationParameters, currspins::Matrix{Int8}, currpositions::Matrix{Float64}, randomflips=nothing)
    
    if isnothing(randomflips)
        randomflips = randlinflip
    end

    flips_int::Array{Bool} = zeros(1, simparams.numparticles)
    if simparams.interaction == alignsimple
        # update spins in place from interactions
        flips_int = alignOn2(simparams, currpositions, currspins, false)
    else
        # println(flips_int)
    end
    # println(" ")
    # println(flips_int)
        
    # update spins in place randomly from stochastic noise
    # flips_noise::Array{Bool} = (x -> randomflips(simparams.dt, simparams.fliprate)).(currspins)

    flips_noise::Array{Bool} = Array{Bool}(undef, 1, simparams.numparticles)
    for p in 1:simparams.numparticles
        flips_noise[p] = randomflips(simparams.dt, simparams.fliprate)
    end

    # println(flips_noise)
    # println(flips_noise)
    totalflips = flips_int .|| flips_noise

    return totalflips
end


