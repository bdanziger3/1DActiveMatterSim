using CSV

include("../simulation/basic_sim.jl")
include("../file_management/simdata_files.jl")
include("./plot_lib.jl")


N::Int = 1000
boxwidth::Real = 1
T::Real = 5
v0::Real = 1
dt::Real = 0.1
totaltime::Real = 100
nsims::Int = 15

simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)

dataoutfile = "./plot_outputs/Apr25/"
# # println(simdata.wrappedpositions)
df = "basic_N10_t1000.0_align_T1_sim_simdata.txt"

simdata = loadsim(df, sequentialtxt)


# handle plotting
wrappedplotfile = "$(dataoutfile)/-int-N-$(simdata.simparams.numparticles)-t$(simdata.simparams.totaltime).pdf"
unwrappedplotfile = "$(dataoutfile)/-int-N-$(simdata.simparams.numparticles)-t$(simdata.simparams.totaltime)-uw.pdf"
plottitle = "Density of $(N) Particles with Alignment Interactions\n T = $(T), dt = $(dt)"
# simpleplotvertxy(simdata, "./$(T)-$(nsims).pdf", plottitle)
densityplot(simdata, wrappedplotfile, plottitle)
densityplot(simdata, unwrappedplotfile, "$(plottitle)-unwrapped", false)


