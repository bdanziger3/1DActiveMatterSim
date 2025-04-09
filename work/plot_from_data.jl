using CSV

include("./min_sim.jl")


N::Int = 1000
boxwidth::Real = 1
T::Real = 5
v0::Real = 1
dt::Real = 0.1
totaltime::Real = 100
nsims::Int = 15

simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)

dataoutfile = "./"

# # println(simdata.wrappedpositions)
simdatafile_prefix = "$(dataoutfile)/multi_sim_$(N)-$(T)-$(v0)"
# df = inou.File("$(simdatafile_prefix)_x.txt")

mx::Matrix{Real} = readdlm("$(simdatafile_prefix)_x.txt", ',', Float64)
mwrappedx::Matrix{Real} = readdlm("$(simdatafile_prefix)_xwrap.txt", ',', Float64)
ms::Matrix{Real} = readdlm("$(simdatafile_prefix)_S.txt", ',', Int8)
times::Matrix{Real} = readdlm("$(simdatafile_prefix)_t.txt", ',', Float64)



simdata = SimulationData(simparams, times, mx, mwrappedx, ms)



# handle plotting
wrappedplotfile = "$(dataoutfile)/loaded-$(T)-$(nsims)_d.pdf"
unwrappedplotfile = "$(dataoutfile)/loaded-$(T)-$(nsims)_d-uw.pdf"
plottitle = "Density of $(N) Particles \n T = $(T), dt = $(dt)"
# simpleplotvertxy(simdata, "./$(T)-$(nsims).pdf", plottitle)
densityplot(simdata, wrappedplotfile, plottitle)
densityplot(simdata, unwrappedplotfile, "$(plottitle)-unwrapped", false)


