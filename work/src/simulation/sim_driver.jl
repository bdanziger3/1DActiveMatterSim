using CSV
using Dates

include("./basic_sim.jl")
include("../file_management/simdata_files.jl")
# include("./plot_lib.jl")
# include("./activesim1d.jl")

# import ActiveSim1D
# using ActiveSim1D

data_dir = "../../data"


# println(length(simdata.positions[:,1]))
# println(length(simdata.times))
# setup simulation params and run the simulation
N::Int = 5e3
boxwidth::Real = 100
fliprate::Real = 1
v0::Real = 1
dt::Real = 1e-4
totaltime::Real = .03
interaction = alignsimple
interactionfliprate = 100
randomstarts = true
snapshot_dt = 1e-2
simparams = SimulationParameters(N, totaltime, dt, v0, fliprate, boxwidth, interaction, interactionfliprate, 0, randomstarts, snapshot_dt)


# STARTING SIM
# data_file = "$(getdatedir())/N$(N)-B$(boxwidth)-$(interaction)-$(interactionfliprate)-t$(totaltime)-sn$(snapshot_dt).txt"
# simdata = runsim(simparams)

# savesim(simdata, data_file, rwtserialized)

# println("Saved first $(totaltime)s at $(now())")


data_file = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/26-6/N5000-B100-alignsimple-100-t4-sn0.01_comp.txt"

# EXTENDING SIM
extended_sd = extendsim(data_file, totaltime, true)
appendsim(extended_sd, data_file, true)
println("Saved extended sim 1 at $(now())")

# for i in 1:10
#     extended_sd = extendsim(data_file, totaltime)
#     appendsim(extended_sd, data_file)
#     println("Saved extended sim $(i) at $(now())")
# end



# collapse file
# collapsesegments(data_file, data_file)






# savesim(simdata, "1basic_N$(N)_t$(floor(totaltime / dt))_align_T$(fliprate)_sim", sequentialtxt)
# savesim(simdata, "testnewserial", sequentialtxt, true)
# simdata2 = loadsim("testnewserial_simdata.txt", sequentialtxt)

# simdata2 = deepcopy(simdata)

# println(simdata.simparams == simdata2.simparams)
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
