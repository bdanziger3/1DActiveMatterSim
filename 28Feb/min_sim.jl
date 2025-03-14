using Random
using PyPlot
using PyCall
using DelimitedFiles
using DataStructures


# Minumum approach to the same code I had before

"""
Simulation parameters on which to run particles
"""
struct SimulationParameters
    numparticles::Int
    totaltime::Real
    dt::Real
    v0::Real
    temp::Real
    boxwidth::Real
    starttime::Real
    # startinglocations::Array{Real}
    # randomstarts::Bool
    # randomspins::Bool

    function # Inner constructor with default values
        SimulationParameters(numparticles, totaltime, dt, v0, temp, boxwidth=1, starttime=0)
        return new(numparticles, totaltime, dt, v0, temp, boxwidth, starttime)
    end
end

struct SimulationData
    simparams::SimulationParameters
    times::Array{Real}
    positions::Matrix{Real}
    wrappedpositions::Matrix{Real}
    spins::Matrix{Int8}
end


function randlinflip(dt::Real, temp::Real)
    return rand() < (dt * temp)
end

function wrap(positions::Array{Real}, boxwidth::Real)
    cwwrap = x -> x - (boxwidth * round(x / boxwidth))
    return cwwrap.(positions)
end


function runsim(simparams::SimulationParameters)

    # set initial particle states
    currpositions::Array{Real} = zeros(1, simparams.numparticles)
    currwrappedpositions::Array{Real} = zeros(1, simparams.numparticles)
    currspins::Array{Int8} = fill!(Array{Int8}(undef, 1, simparams.numparticles), 1)
    
    times = collect(simparams.starttime:simparams.dt:simparams.totaltime)
    nsteps = length(times) - 1

    mx::Array{Real} = zeros(nsteps+1, simparams.numparticles)
    mwrappedx::Array{Real} = zeros(nsteps+1, simparams.numparticles)
    ms::Array{Real} = zeros(nsteps+1, simparams.numparticles)

    mx[1,:] = copy(currpositions)
    mwrappedx[1,:] = copy(currpositions)
    ms[1,:] = copy(currspins)

    # run simulation `nsims` times
    for i in 1:nsteps

        # update positions in place
        currpositions .+= (currspins .* simparams.v0 * simparams.dt)

        # handle wrapping of positions (periodic boundary conditions)
        currwrappedpositions = wrap(currpositions, simparams.boxwidth)

        # update spins in place
        flips::Array{Bool} = (x -> randlinflip(simparams.dt, simparams.temp)).(currspins)
        currspins[flips] .*= -1


        mx[i+1,:] = copy(currpositions)
        mwrappedx[i+1,:] = copy(currpositions)
        ms[i+1,:] = copy(currspins) 
    end

    simdata = SimulationData(simparams, times, mx, mwrappedx, ms)
    return simdata
end


# get particle densities
function particledensity(simdata::SimulationData, nbins::Real=100)
    n = simdata.simparams.numparticles
    times = simdata.times
    wrappedx = simdata.wrappedpositions

    binwidth = simdata.simparams.boxwidth / nbins
    xgrid::Array{Real} = collect(-simdata.simparams.boxwidth/2:binwidth:simdata.simparams.boxwidth/2)

    densities::Matrix{Real} = zeros(wrappedx.size[1], length(xgrid))

    bin = (x -> round(x / binwidth) + round(length(xgrid) / 2))

    for i in 1:wrappedx.size[1]

        densitymap = counter(bin.(wrappedx[i,:]))

        bins = collect(1:length(xgrid))
        density = (b -> densitymap[b]).(bins)

        densities[i,:] = density ./ n

    end


    return xgrid, densities

end


    

    # matplotlib.pyplot.close()


# setup simulation params and run the simulation
N::Int = 1000
boxwidth::Real = 1
T::Real = 5
v0::Real = 1
dt::Real = 0.01
totaltime::Real = 30
simparams = SimulationParameters(N, totaltime, dt, v0, T, boxwidth)
nsims = 15

simdata = runsim(simparams)





# println(simdata.wrappedpositions)

dataoutfile = "./dataout1"
writedlm("./multi_sim_$(N)-$(T)-$(v0)_x.txt", simdata.positions, ",")
writedlm("./multi_sim_$(N)-$(T)-$(v0)_xwrap.txt", simdata.wrappedpositions, ",")
writedlm("./multi_sim_$(N)-$(T)-$(v0)_S.txt", simdata.spins, ",")
writedlm("./multi_sim_$(N)-$(T)-$(v0)_t.txt", simdata.times, ",")



function set_visual(xlabel = "", ylabel = "",size= (3,2)) # size = (3,2) if figure fills half a page width
    plt.clf() # clears old plots
    fig, ax = plt.subplots(figsize = (size[1],size[2])) # creates plot
    fig = plt.gcf()
    fig.patch.set_alpha(0) # makes area outside of plot transparent

    ax.set_xlabel(xlabel) # sets axis labels
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="y", direction="in", which = "both") # makes ticks face inward
    ax.tick_params(axis="x", direction="in", which = "both")
    return fig,ax # returns figure and axis such that it can be manipulated outside of the function
end

function simpleplotvertxy(simdata::SimulationData, filepath::String = "", title::String = "")
    fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6,10))

    
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(title)   

    for p in 1:simdata.simparams.numparticles
        plt.plot(simdata.positions[:,p], simdata.times)
    end

    plt.xlim(-simdata.simparams.boxwidth / 2, simdata.simparams.boxwidth / 2)
    ax.invert_yaxis()

    # if filepath == ""
    #     filepath = makedefaultfilename()
    # end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 



function densityplot(simdata::SimulationData, filepath::String = "", title::String = "", wrapped::Bool=true, spacebins::Int=100)
    # make particle density plot

    wrappedx = simdata.wrappedpositions
    xlims = [-simdata.simparams.boxwidth / 2, simdata.simparams.boxwidth / 2]
    if !wrapped
        wrappedx = simdata.positions

        m = maximum([abs(minimum(wrappedx)), abs(maximum(wrappedx))])
        xlims = [-m, m]
    end
    xpoints = wrappedx.size[2]
    times = simdata.times
    ntimes = length(times)
    singlelistt = zeros(length(wrappedx))
    singlelistx = zeros(length(wrappedx))

    for i in 1:length(times)
        singlelistt[(i-1)*xpoints+1:(i)*xpoints] = times[i] * ones(xpoints)
        singlelistx[(i-1)*xpoints+1:(i)*xpoints] = wrappedx[i,:]
    end

    # println(singlelistt)
    # println(singlelistx)
    # println(wrappedx)


    fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6,10))

    
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(title)   

    plt.hist2d(singlelistx, singlelistt, bins=[spacebins, length(times)], cmap="Blues")

    plt.xlim(xlims)
    plt.clim([0, simdata.simparams.numparticles / 10])
    ax.invert_yaxis()
    plt.colorbar(label="Density")

    # if filepath == ""
    #     filepath = makedefaultfilename()
    # end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 


# handle plotting
plottitle = "Density of $(N) Particles \n T = $(T), dt = $(dt)"
# simpleplotvertxy(simdata, "./$(T)-$(nsims).pdf", plottitle)
densityplot(simdata, "./$(T)-$(nsims)_d.pdf", plottitle)
densityplot(simdata, "./$(T)-$(nsims)_d-uw.pdf", "$(plottitle)-unwrapped", false)


# println(length(simdata.positions[:,1]))
# println(length(simdata.times))

# function periodicplotvertxy(simdata::ParticleSimData, filepath::String = "", title::String = "")
#     boundary_crosses = periodicboundary_paths(simdata)
#     println(boundary_crosses)

#     fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6, 10))

    
#     ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    
#     plt.title(title)   

#     for i in 1:length(simdata.simparams.numparticles)
#         path = simdata.positions[:,i]
#         plt.plot(path.positions, path.times)
#     end
        
#     for i in keys(boundary_crosses)
#         path = simdata.particlepaths[i]
#         # Access the color of the last added line
#         line_color = ax.get_lines()[i].get_color()

#         # Add line with the same color to the left and right
#         for j in 1:boundary_crosses[i][1]
#             plot(path.positions .+ (j * simdata.boxwidth), path.times, color=line_color)
#         end
#         for j in 1:boundary_crosses[i][2]
#             plot(path.positions .- (j * simdata.boxwidth), path.times, color=line_color)
#         end
#     end
        
#     plt.xlim(-simdata.boxwidth / 2, simdata.boxwidth / 2)
#     ax.invert_yaxis()

#     if filepath == ""
#         filepath = makedefaultfilename()
#     end

#     plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

#     println("plot saved in current directory as $(filepath)")
# end 



# periodicplotvertxy(simdata, "./plots/spnv$(T)-$(nsims)_per.pdf", plottitle)

# for j in 1:nsims
#     xarray::Array{Real} = []

#     for i in 1:1000
#         timestep!(p, .01)
#         push!(xarray, p.x)
#     end

#     plt.plot(1:length(xarray), xarray)

# end


# plt.savefig("single_particle.pdf", bbox_inches = "tight",pad_inches=0.01)

