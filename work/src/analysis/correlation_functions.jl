using CSV

include("../simulation/basic_sim.jl")



"""
Mean Squared Displacement

Calculates the squared distance a particle travels over a timestep dt
averaged over all particles and over all intervals of length dt.

Returns a distribution that gives the MSD (in units of length^2) as a function of dt (in units of time)

Calculated from the (unrwrapped) positions
"""
function meansqdisp(simdata::SimulationData)::Matrix{<:Real}
    # returns the Mean Squared Displacement data from a SimulationData object

    times = gettimes(simdata.simparams)
    ntimes = getntimes(simdata.simparams)
    positions = copy(simdata.positions)

    mintimestep::Float64 = simdata.simparams.dt
    maxtimestep::Float64 = simdata.simparams.totaltime / 2
    dtarray::Array{Float64} = collect(mintimestep:mintimestep:maxtimestep)
    ndts::Int64 = size(dtarray)[1]

    msdmat::Matrix{Float64} = zeros(2, ndts+1)

    msdmat[1, 1] = 0
    msdmat[2, 1] = 0


    dt = mintimestep

    # evaluate for each interval dt, incremented by number of `mintimestep`s 
    for ndt in 1:ndts
        
        maxt0 = ntimes - ndt

        # integrate over all starting points t0
        t0dispsum = 0
        for t0 in 1:maxt0


            # average over all particles
            particledispsum = 0
            for n in 1:simdata.simparams.numparticles
                particledispsum += (positions[t0+ndt, n] - positions[t0, n])^2
            end
            particledispsum /= simdata.simparams.numparticles

            t0dispsum += particledispsum

        end
        t0dispsum /= maxt0


        msdmat[1, ndt+1] = ndt * mintimestep
        msdmat[2, ndt+1] = t0dispsum
    end




    # for each winow of length dt

    # for each particle


    return msdmat
end




"""
Polar Order

Calculates the average spin over all particles as a function of time
"""
function polarorder(simdata::SimulationData)::Matrix{Float64}
    # returns the Polar Order as a function of time throughout one simulation

    times = gettimes(simdata.simparams)
    ntimes = getntimes(simdata.simparams)
    spins = copy(simdata.spins)
    
    polarordermat::Matrix{Float64} = zeros(ntimes, 2)

    for t in 1:ntimes
        polarordermat[t, :] = [times[t], sum(spins[t, :]) / simdata.simparams.numparticles]
    end

    return polarordermat
end

"""
Polar Order Window

Calculates the average spin within a radius around each point in space as a function of space and time
"""
function polarorderwindow(simdata::SimulationData, radius::Float64)::Matrix{Float64}
    # returns the Polar Order as a function of time throughout one simulation

    times = gettimes(simdata.simparams)
    ntimes = getntimes(simdata.simparams)
    wrappedpos = copy(simdata.wrappedpositions)
    spins = copy(simdata.spins)
    dx = abs(simdata.simparams.v0 * simdata.simparams.dt)

    minx = (-simdata.simparams.boxwidth / 2) + (dx / 2)
    maxx = (simdata.simparams.boxwidth / 2) - (dx / 2)

    xs = collect(minx:dx:maxx)
    nx = size(xs)[1]
    
    polarordermat::Matrix{Float64} = zeros(ntimes, nx)

    for t in 1:ntimes
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
function ocf(simdata::SimulationData)::Matrix{Float64}
    # returns the Orientation Corelation as a function of time throughout one simulation


    times = gettimes(simdata.simparams)
    ntimes = getntimes(simdata.simparams)
    spins = copy(simdata.spins)

    mintimestep::Float64 = simdata.simparams.snapshot_dt
    maxtimestep::Float64 = simdata.simparams.totaltime / 2
    dtarray::Array{Float64} = collect(mintimestep:mintimestep:maxtimestep)
    ndts::Int64 = size(dtarray)[1]

    ocmat::Matrix{Float64} = zeros(2, ndts+1)

    dt = mintimestep
    t0 = 1

    # evaluate for each interval dt, incremented by number of `mintimestep`s 
    for ndt in 0:ndts
        
        # for now only calc at one t0
        t0ocsum = 0

        # average over all particles
        particleocsum = 0
        for n in 1:simdata.simparams.numparticles
            particleocsum += (spins[t0+ndt, n] * spins[t0, n])
        end
        particleocsum /= simdata.simparams.numparticles

        t0ocsum += particleocsum

        


        ocmat[1, ndt+1] = ndt * mintimestep
        ocmat[2, ndt+1] = t0ocsum
    end

    return ocmat
end




    
