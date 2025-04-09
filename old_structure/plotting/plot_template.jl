using PyPlot
using PyCall

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
    fig, ax = set_visual(L"distance $x/\sigma$",L"velocity $v/v_0$") # L"" for latex
    x = 0:0.01:20
    plt.plot(x,1 .+ x.^2)

    plt.savefig("velocity_over_distance.pdf", bbox_inches = "tight",pad_inches=0.01) 
end

function simpleplotxy(xdata, ydata, filepath::String = "", clear::Bool = true)
    if clear
        set_visual("Timestep", L"$x$ location", clear)
    end

    plt.plot(xdata, ydata)
    if filepath == ""
        filepath = makedefaultfilename()
    end

    plt.savefig(filepath, bbox_inches = "tight",pad_inches=0.01)

    println("plot saved in current directory as $(filepath)")
end 

function makedefaultfilename()
    currentfiles = readdir(".")
    defaultname = "myplot.pdf"

    increment::Int32 = 1
    while defaultname in currentfiles
        defaultname = "myplot$(increment).pdf"
    end

    return defaultname
end
