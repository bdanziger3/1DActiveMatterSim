using CSV

include("./min_sim.jl")




function meansqdisp(simdata::SimulationData)
    # returns the Mean Squared Displacement data from a SimulationData object

    times = copy(simdata.times)
    positions = copy(simdata.positions)

    maxtimestep::Real = simdata.totaltime / 2
    simlength::Real = simdata.simparams.totaltime - simdata.simparams.starttime
    timesteps::Int = collect(0:simdata.simparams.dt:simlength)


end


N::Int = 1000
boxwidth::Real = 1
T::Real = 5
v0::Real = 1
dt::Real = 0.1
totaltime::Real = 100
nsims::Int = 15

simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)


# load in data from file
dataoutfile = "./"

simdatafile_prefix = "$(dataoutfile)/multi_sim_$(N)-$(T)-$(v0)"

mx::Matrix{Real} = readdlm("$(simdatafile_prefix)_x.txt", ',', Float64)
mwrappedx::Matrix{Real} = readdlm("$(simdatafile_prefix)_xwrap.txt", ',', Float64)
ms::Matrix{Real} = readdlm("$(simdatafile_prefix)_S.txt", ',', Int8)
times::Matrix{Real} = readdlm("$(simdatafile_prefix)_t.txt", ',', Float64)

simdata = SimulationData(simparams, times, mx, mwrappedx, ms)


