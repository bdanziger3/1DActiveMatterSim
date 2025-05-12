using CSV

include("./basic_sim.jl")
include("./simdata_files.jl")
include("./plot_lib.jl")


# println(length(simdata.positions[:,1]))
# println(length(simdata.times))
# setup simulation params and run the simulation
N::Int = 100
boxwidth::Real = 1
ctime::Real = 3
fliprate::Real = 1/ctime
v0::Real = .5 * boxwidth/ctime
dt::Real = 0.001
totaltime::Real = 100
interaction = alignsimple
simparams = SimulationParameters(N, totaltime, dt, v0, fliprate, boxwidth, interaction)
nsims = 15

# println(simparams == simparams2)
# println(simparams)
# j1 = json_serialize(simparams)
# j2 = json_deserialize_simparams(j1)

# println(simparams)
# println(j1)
# println(j2)

simdata = runsim(simparams)

# savesim(simdata, "1basic_N$(N)_t$(floor(totaltime / dt))_align_T$(fliprate)_sim", sequentialtxt)
savesim(simdata, "basic_N$(N)_t$(floor(totaltime / dt))_interaction_$(interaction)_T$(fliprate)_sim", sequentialtxt, true)
# simdata = loadsim("txt_data_files/basic_N$(N)_t$(Int64(floor(totaltime / dt)))_noint_sim_simdata.txt", sequentialtxt)

# simdata2 = deepcopy(simdata)

# println(simdata == simdata2)
# println(simdata)
# j3 = json_serialize(simdata)
# s4 = json_deserialize_simdata(j3)
# println(s4)
# simdataDf = simdata2df(simdata)



# simdata2 = df2simdata(simdataDf, simparams)
# simdata3 = df2simdata(simdataDf)

# println(simdata2)
# println(simdata3)


# println(simdata.df)
# println(simdataDf)

# dfoutfile = "./df_dataout.csv"
# dataoutfile = "./plot_outputs/Apr9"

# CSV.write(dfoutfile, simdataDf)



# df2 = DataFrame(CSV.File(dfoutfile))

# df2[1,1] = 400


# for i in 1:size(simdataDf)[2]
#     if simdataDf[!,i] != df2[!,i]
#         println("Not equal in column $(i)")
#         break
#     end
#     if i == size(simdataDf)[2]
#         println("Equal")
#     end
# end

# println(simdataDf == df2)
# println(simdataDf === df2)
# println(isequal(simdataDf, df2))

# handle plotting
# wrappedplotfile = "$(dataoutfile)/$(T)-$(nsims).png"
# unwrappedplotfile = "$(dataoutfile)/$(T)-$(nsims)-uw.png"
# plottitle = "Paths of $(N) Particles \n T = $(T), dt = $(dt)"
# # simpleplotvertxy(simdata, "./$(T)-$(nsims).pdf", plottitle)
# # densityplot(simdata, wrappedplotfile, plottitle)
# # densityplot(simdata, unwrappedplotfile, "$(plottitle)-unwrapped", false)
# simpleplotvertxy(simdata, wrappedplotfile, plottitle)
# simpleplotvertxy(simdata, unwrappedplotfile, "$(plottitle)-unwrapped")


# # # println(simdata.wrappedpositions)
# simdatafile_prefix = "$(dataoutfile)/multi_sim_$(N)-$(T)-$(v0)"
# writedlm("$(simdatafile_prefix)_x.txt", simdata.positions, ",")
# writedlm("$(simdatafile_prefix)_xwrap.txt", simdata.wrappedpositions, ",")
# writedlm("$(simdatafile_prefix)_S.txt", simdata.spins, ",")
# writedlm("$(simdatafile_prefix)_t.txt", simdata.times, ",")



# savesim(simdata, "testfunc2", sequentialtxt, true)
# datafile_path = "testfunc2_simdata.txt"
# loadedsim = loadsim(datafile_path, sequentialtxt)

# savesim(loadedsim, "testfunc3", sequentialtxt, true)
