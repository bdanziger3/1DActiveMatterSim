using PyPlot
using PyCall

include("../data_types/simdata.jl")



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

function testplot()
    fig, ax = set_visual(L"Space ($x$)",L"Position over Time")
    x = 0:0.01:20
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(L"Space ($x$)")
    plt.plot(sin.(x), x)
    ax.invert_yaxis()
    fig.patch.set_facecolor("white")
    plt.savefig("verticaltest1.pdf", bbox_inches = "tight",pad_inches=0.01) 
end


function simpleplotvertxy(xdata::Array{Real}, ydata::Array{Real}, filepath::String = "", clear::Bool = true)
    if clear
        set_visual(L"Space ($x$)", L"Time ($s$)")
    end

    fig, ax = ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(L"Space ($x$)")    
    plt.plot(ydata, xdata)
    ax.invert_yaxis()

    if filepath == ""
        filepath = makedefaultfilename()
    end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 

function simpleplotvertxy(simdata::ParticleSimData, filepath::String = "", title::string = "")
    fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6,10))

    
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    plt.title(title)   

    for path::ParticlePath in simdata.particlepaths
        plt.plot(path.positions, path.times)
    end

    plt.xlim(-simdata.boxwidth / 2, simdata.boxwidth / 2)
    ax.invert_yaxis()

    if filepath == ""
        filepath = makedefaultfilename()
    end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 


function periodicplotvertxy(simdata::ParticleSimData, filepath::String = "", title::String = "")
    boundary_crosses = periodicboundary_paths(simdata)
    println(boundary_crosses)

    fig, ax = set_visual(L"Space ($x$)", L"Time ($s$)", (6, 10))

    
    ax.tick_params(top=true, labeltop=true, bottom=false, labelbottom=false)
    
    plt.title(title)   

    for i in 1:length(simdata.particlepaths)
        path = simdata.particlepaths[i]
        plt.plot(path.positions, path.times)
    end
        
    for i in keys(boundary_crosses)
        path = simdata.particlepaths[i]
        # Access the color of the last added line
        line_color = ax.get_lines()[i].get_color()

        # Add line with the same color to the left and right
        for j in 1:boundary_crosses[i][1]
            plot(path.positions .+ (j * simdata.boxwidth), path.times, color=line_color)
        end
        for j in 1:boundary_crosses[i][2]
            plot(path.positions .- (j * simdata.boxwidth), path.times, color=line_color)
        end
    end
        
    plt.xlim(-simdata.boxwidth / 2, simdata.boxwidth / 2)
    ax.invert_yaxis()

    if filepath == ""
        filepath = makedefaultfilename()
    end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 




