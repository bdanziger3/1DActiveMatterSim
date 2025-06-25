using CSV

include("./test_fixtures.jl")


# println(length(simdata.positions[:,1]))
# println(length(simdata.times))
# setup simulation params and run the simulation
N::Int = 3
boxwidth::Real = 1
ctime::Real = 100
fliprate::Real = 1
v0::Real = 1 #.5 * boxwidth / 3
dt::Real = 0.0001
totaltime::Real = .001
interaction = alignsimple
interactionfliprate = 300
randomstarts = true
simparams = SimulationParameters(N, totaltime, dt, v0, fliprate, boxwidth, interaction, interactionfliprate, 0, randomstarts)
# nsims = 15


"""
1. Run a simulation and load the data.
2. Reload the saved data and extend it
3. Combine the two data structs and append the files together
"""
function test_extendsim(nsegments::Int=2)
    @assert nsegments >= 2 "Test only works for 2 or more segments."

    timeextension = .0005
    sd1 = runsim(simparams_10)
    sdlist::Array{SimulationData} = Array{SimulationData}[]
    push!(sdlist, sd1)

    filestring = "./test_data/test_extendsim_sd1"
    file1 = "$(filestring).txt"
    savesim(sd1, file1, rowwisetxt, true)
    
    # extend segments to get to `nsegments` total segments
    for s in 2:nsegments
        file_s = "$(filestring)_ext_$(s).txt"
        sd2 = extendsim(file1, timeextension)
        savesim(sd2, file_s, rowwisetxt, true)
        appendsim(sd2, file1)
        push!(sdlist, sd2)

        # check params
        @assert sdlist[s].simparams.starttime == sdlist[s-1].simparams.totaltime
        @assert sdlist[s].simparams.totaltime == timeextension
        @assert sdlist[s].simparams.randomstarts == false
    
        sp1 = asarray(sd1.simparams)
        sp2 = asarray(sd2.simparams)

        @assert sp1[1] == sp2[1] "Sim Params of segment $(s) don't match the original for `$(fieldnames(SimulationParameters)[1])`."
        for i in 3:8
            @assert sp1[i] == sp2[i] "Sim Params of segment $(s) don't match the original for `$(fieldnames(SimulationParameters)[i])`."
        end
    
    end
end


test_extendsim(2)


println("PASSED")