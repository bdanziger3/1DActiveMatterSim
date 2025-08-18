using CSV
using Statistics

include("../simulation/basic_sim.jl")
include("./histogram_helper.jl")


@enum SweepType densitysweep boxwidthsweep interactionstrengthsweep nosweep

"""
Mean Squared Displacement

Calculates the squared distance a particle travels over a timestep dt
averaged over all particles and over all intervals of length dt.

Returns a distribution that gives the MSD (in units of length^2) as a function of dt (in units of time)

Calculated from the (unrwrapped) positions
"""
function meansqdisp(positions::Matrix{<:Float64}, simparams::SimulationParameters, settletime::Real=-1, maxdt::Real=-1)::Matrix{<:Real}
    # returns the Mean Squared Displacement data from a SimulationData object

    # if no `settletime` specified, use half the `totaltime`
    if settletime < 0
        settletime = simparams.totaltime / 2
    end

    nsaves = getnsaves(simparams)
    positions = unwrappositions(positions, simparams)
    
    mintimestep::Float64 = simparams.snapshot_dt
    
    # if no `maxdt` provided, defaults to maximum distance between `settletime` and `totaltime`
    local maxtimestep::Float64
    if maxdt < 0
        maxtimestep = simparams.totaltime - settletime
    else
        maxtimestep = maxdt
    end

    ndts::Int64 = Int64(floor(maxtimestep / mintimestep))
    
    # calculate snapshot index of the settletime
    settletime_i = nsaves - ndts

    msdmat::Matrix{Float64} = zeros(2, ndts+1)

    # first column of `msdmat` is 0s

    # evaluate for each interval dt, incremented by number of `mintimestep`s 
    @showprogress 1 "Calculating MSD..." for ndt in 1:ndts
        
        mint0_i = settletime_i
        maxt0_i = nsaves - ndt
        nt0s = maxt0_i - mint0_i + 1

        # integrate over all starting points t0
        t0dispsum = 0
        for t0_i in mint0_i:maxt0_i
            # average over all particles
            particledispsum = msd(positions[t0_i+ndt, :], positions[t0_i, :])
            t0dispsum += particledispsum
        end
        t0dispsum /= nt0s

        # fill in top row as timestep `dt` and bottom row as Mean Squared Displacement
        msdmat[1, ndt+1] = ndt * mintimestep
        msdmat[2, ndt+1] = t0dispsum
    end

    return msdmat
end


"""
Mean Squared Displacement

Calculates the squared distance a particle travels over a timestep dt
averaged over all particles and over all intervals of length dt.

Returns a distribution that gives the MSD (in units of length^2) as a function of dt (in units of time)

Calculated from the (unrwrapped) positions
"""
function meansqdisp(simdata::SimulationData, settletime::Real=-1.0, maxdt::Real=-1.0)::Matrix{<:Real}
    return meansqdisp(simdata.positions, simdata.simparams, settletime, maxdt)
end




"""
Mean Spin

Calculates the average spin over all particles as a function of time
"""
function meanspin(simdata::SimulationData)::Matrix{Float64}
    # returns the Mean Spin as a function of time throughout one simulation

    times = getsavetimes(simdata.simparams)
    nsaves = getnsaves(simdata.simparams)
    spins = copy(simdata.spins)
    
    meanspinmat::Matrix{Float64} = zeros(2, nsaves)

    # for t in 1:nsaves
    @showprogress 1 "Calculating Mean Spin..." for t in 1:nsaves
        meanspinmat[:, t] = [times[t]; sum(spins[t, :]) / simdata.simparams.numparticles]
    end

    return meanspinmat
end

"""
Mean Spin but with different parameters

Calculates the average spin over all particles as a function of time
"""
function meanspin(spins::Matrix{Int8}, simparams::SimulationParameters)::Matrix{Float64}
    # returns the Mean Spin as a function of time throughout one simulation
    times = getsavetimes(simparams)
    nsaves = getnsaves(simparams)
    
    meanspinmat::Matrix{Float64} = zeros(2, nsaves)

    # for t in 1:nsaves
    @showprogress 1 "Calculating Mean Spin..." for t in 1:nsaves
        meanspinmat[:, t] = [times[t]; sum(spins[t, :]) / simparams.numparticles]
    end

    return meanspinmat
end

"""
Polar Order Window

Calculates the average spin within a radius around each point in space as a function of space and time
"""
function polarorderwindow(simdata::SimulationData, radius::Float64)::Matrix{Float64}
    # returns the Polar Order as a function of time throughout one simulation

    nsaves = getnsaves(simdata.simparams)
    wrappedpos = copy(simdata.positions)
    spins = copy(simdata.spins)
    dx = abs(simdata.simparams.v0 * simdata.simparams.dt)

    minx = (-simdata.simparams.boxwidth / 2) + (dx / 2)
    maxx = (simdata.simparams.boxwidth / 2) - (dx / 2)

    xs = collect(minx:dx:maxx)
    nx = size(xs)[1]
    
    polarordermat::Matrix{Float64} = zeros(nsaves, nx)

    for t in 1:nsaves
        for x in 1:nx
            nparticles = 0
            spinsum = 0
            # loop over all particles and find polar order of particles within r of xarray

            for n in 1:simdata.simparams.numparticles
                if abs(wrappedpos[t, n] - xs[x]) <= radius
                    nparticles += 1
                    spinsum += spins[t, n]
                end
            end
            if nparticles > 0
                spinsum /= nparticles
            else
                spinsum = 0
            end
            
            polarordermat[t, x] = spinsum
        end
    end

    return polarordermat
end

"""
Orientation Correlation

C_p(t)

Calculates the average Correlation between a particle's orientation and its orientation at a a time t later

Returns a matrix givng (t,C_p(t)) data
"""
function orientationcorrelation(spins::Matrix{Int8}, simparams::SimulationParameters, settletime::Real=-1, maxdt::Real=-1)::Matrix{Float64}
    # returns the Orientation Corelation as a function of time throughout one simulation

    nsaves = getnsaves(simparams)

    # if no `settletime` provided, defaults to half of the total time
    if settletime == -1
        settletime = simparams.totaltime / 2
    end

    
    # find step number corresponding to settletime
    settleindex = getindexoftime(simparams, settletime)
    
    mintimestep::Float64 = simparams.snapshot_dt
    
    # if no `maxdt` provided, defaults to maximum distance between `settletime` and `totaltime`
    local maxtimestep::Float64
    if maxdt == -1
        maxtimestep = simparams.totaltime - settletime
    else
        maxtimestep = maxdt
    end
    
    dtarray::Array{Float64} = collect(mintimestep:mintimestep:maxtimestep)
    ndts::Int64 = size(dtarray)[1]

    ocmat::Matrix{Float64} = zeros(2, ndts+1)

    dt = mintimestep
    t0 = settleindex


    # evaluate for each interval dt, incremented by number of `mintimestep`s
    # this is the x-axis data for the Orientation Correlation Function
    @showprogress 1 "Calculating Orientation Self-Correlation Densities..." for ndt in 0:ndts        
        # reset `t0` back to first index we're measuring
        t0 = settleindex

        # this is the y-axis data for the Orientation Correlation Function
        t0ocsum = 0

        # for each time interval dt, sweep and average over all possible start times t0
        # i.e. sliding window integration average
        # calculate number of start times
        nt0s = (nsaves - settleindex) - (ndt - 1)
        for t0_i in 0:(nt0s - 1)
            t0ocsum += sum(spins[settleindex + t0_i, :] .* spins[settleindex + t0_i + ndt, :])
        end

        # now divide by total number of products taken to normalize to +-1.
        t0ocsum /= (simparams.numparticles * nt0s)

        ocmat[1, ndt+1] = ndt * mintimestep
        ocmat[2, ndt+1] = t0ocsum
    end

    return ocmat
end


"""
Orientation Correlation for different parameters
"""
function orientationcorrelation(simdata::SimulationData, settletime::Real=-1, maxdt::Real=-1)::Matrix{Float64}
    return orientationcorrelation(simdata.spins, simdata.simparams, settletime, maxdt)
    
end


"""
Takes a list of positions and bins them based on the `binsettings`,
then counts the number of particles in each bin and divides by the width of the bin.
Returns an array of particle densities for each bin with the first bin starting at -boxwidth/2 and the final bin ending at boxwidth/2.

By default, runs the analysis on the entire simulation data.
set `timeindex` to the index of the time in the simulation to get the data for a single frame.
"""
function posdensity_fixedbins(simdata::SimulationData, timeindex::Int64=-1)::Array{Float64}
    # binwidth with respect to the minimum space step `dx`
    dx_per_bin::Float64 = 10
    posbinsettings::BinSettings = positionbinsettings(simdata.simparams.boxwidth, dx_per_bin * snapshot_dx(simdata.simparams))

    # use entire position matrix by default
    local positions::Array{Float64}
    local normalization::Float64
    if timeindex == -1
        positions = simdata.positions
        normalization = getnsaves(simdata.simparams) * posbinsettings.binwidth
    elseif timeindex > 0
        positions = simdata.positions[timeindex, :]
        normalization = posbinsettings.binwidth
    else
        @assert false "Unexpected value for `timeindex`. Must be either -1 or a positive integer."
    end

    # bin the particle positions
    posbincounts::Array{Int64} = bincounts(positions, posbinsettings)

    # convert to particle density (particles / unit of length)
    particledensities = posbincounts ./ normalization

    return particledensities
end

"""
Gets particle densities and then bins that data to get histogram data bin counts
"""
function posdensityhistcounts(simdata::SimulationData, timeindex::Int64=-1)::Array{Float64}
    n_density_bins::Int64 = 100
    
    particledensities = posdensity(simdata, timeindex)::Array{Float64}

    # now bin the particle densities
    densitiesbinsettings = BinSettings(maximum(particledensities), n_density_bins)
    posdensitycounts::Array{Int64} = bincounts(particledensities, densitiesbinsettings)

    return posdensitycounts    
end

"""
Takes a simulation set of positions and, at each timestep, bins them based on the `binsettings`.
Then counts the number of particles in each bin and divides by the width of the bin to get the particle densities of each bin at each timestep.
Remembers the occurence of particle densities regardless of the position of the bin, and counts over all timesteps.
Returns an array of particle densities (of size nsaves x nbins)

"""
function particledensity_alltimes(simdata::SimulationData, binwidth::Float64=1.0)::Array{Float64}
    posbinsettings::BinSettings = positionbinsettings(simdata.simparams.boxwidth, binwidth)

    # bin the particle positions at each timestep
    nsaves = getnsaves(simdata.simparams)
    posbindensities::Matrix{Float64} = zeros(nsaves, posbinsettings.nbins)
    @showprogress 1 "Calculating Particle Densities..." for i in 1:nsaves
        # count particles in each bin and divide by binwidth to get particle densities
        posbindensities[i,:] = bincounts(simdata.positions[i,:], posbinsettings) ./ posbinsettings.binwidth
    end

    return posbindensities
end





"""
Mean Cluster Width

Calculates the average length of clusters of particles that are aligned
Measured in units of number of particles
"""
function meanclusterlength(positions::Matrix{Float64}, spins::Matrix{Int8}, simparams::SimulationParameters)::Matrix{Float64}
    # returns the Mean Spin as a function of time throughout one simulation
    times = getsavetimes(simparams)
    nsaves = getnsaves(simparams)
    
    meanclusterlengthmat::Matrix{Float64} = zeros(2, nsaves)
    meanclusterlengthmat[1, :] = times

    # for t in 1:nsaves
    @showprogress 1 "Calculating Mean Cluster Length..." for t in 1:nsaves
        # sort particles
        sortedindices = sortperm(positions[t, :])
        pos_sorted_spins = spins[t, sortedindices]

        runlengths = []
        currspin = pos_sorted_spins[1]
        currrunlength = 1
        for i in 2:simparams.numparticles
            if currspin == pos_sorted_spins[i]
                # increment if same spin
                currrunlength += 1
            else
                # add run to array
                push!(runlengths, currrunlength)
                currrunlength = 0
                currspin = pos_sorted_spins[i]
            end
        end

        # now add final run
        if currspin != pos_sorted_spins[1]
            push!(runlengths, currrunlength)
        else
            # if last run has same spin as first run, combine them
            runlengths[1] += currrunlength
        end

        meanclusterlengthmat[2, t] = mean(runlengths)
    end

    return meanclusterlengthmat
end

"""
Same function but with `SimulationData` object as an input parameter.
"""
function meanclusterlength(simdata::SimulationData)::Matrix{Float64}
    return meanclusterlength(simdata.positions, simdata.spins, simdata.simparams)
end




    
