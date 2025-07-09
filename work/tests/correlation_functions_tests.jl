using PyPlot
using PyCall

include("../src/analysis/correlation_functions.jl")
include("../src/file_management/simdata_files.jl")
include("../src/utils/.paths.jl")


# run basic sim


function msd_test()
    # open simdata
    simplesim_datafile = "../txt_data_files/basic_N10_t1000.0_noint_T1_sim_simdata.txt"
    # simplesim_datafile = "../txt_data_files/basic_N100_t1000_noint_sim_simdata.txt"
    # simplesim_datafile = "../txt_data_files/basic_N1000_t10000_noint_sim_simdata.txt"

    simdata = loadsim(simplesim_datafile, sequentialtxt)

    # analyze correlation
    msdmat = meansqdisp(simdata)


    # handle plotting

    # output file path
    plotfilepath = "../plot_outputs/Apr9/MSD_N10_t100_T1_noint.pdf"


    plt.plot(msdmat[1,2:end], msdmat[2,2:end])
    plt.title("Mean Squared Displacement\n10 Particles\ntmax = 10")
    plt.xlabel("dt")
    plt.ylabel("Mean Squared Displacement")
    plt.loglog()


    # plt.savefig(plotfilepath, bbox_inches = "tight",pad_inches=0.01)
    plt.show()


end


function polarordertest()
    # open simdata
    simplesim_datafile = "../txt_data_files/basic_N10_t100_noint_sim.txt"
    # simplesim_datafile = "../txt_data_files/basic_N100_t1000_noint_sim_simdata.txt"
    # simplesim_datafile = "../txt_data_files/basic_N1000_t10000_noint_sim_simdata.txt"

    simdata = loadsim(simplesim_datafile, sequentialtxt)



    # analyze correlation
    polarmat = polarorder(simdata)


    # handle plotting

    # output file path
    plotfilepath = "../plot_outputs/Apr9/PO_N10_t100_noint.pdf"
    # plotfilepath = "../plot_outputs/Apr9/PO_N100_t1000_noint.pdf"
    # plotfilepath = "../plot_outputs/Apr9/PO_N1000_t10000_noint.pdf"


    plt.plot(polarmat[:,1], polarmat[:,2])
    plt.title("Polar Order\n10 Particles\ntmax = 10")
    plt.xlabel("t")
    plt.ylabel("Polar Order")
    plt.ylim([-1,1])


    # plt.savefig(plotfilepath, bbox_inches = "tight",pad_inches=0.01)
    plt.show()
end


function polarorderwindowtest()

    # open simdata
    simplesim_datafile = "../txt_data_files/basic_N10_t100_noint_sim.txt"
    # simplesim_datafile = "../txt_data_files/basic_N100_t1000_noint_sim_simdata.txt"
    # simplesim_datafile = "../txt_data_files/basic_N1000_t10000_noint_sim_simdata.txt"

    simdata = loadsim(simplesim_datafile, sequentialtxt)

    # analyze correlation
    polarmat = polarorderwindow(simdata, 0.1)


    println(polarmat)

end



function ocf_test()

    LOG = true

    nparticles = 10000

    test_interaction_file_300 = fixpath("/work/data/22-6/N10000-alignsimple-t1-sn0.01.txt")
    test_interaction_file_150 = fixpath("/work/data/23-6/N10000-alignsimple-150-t0.5-sn0.01.txt")
    test_interaction_file_60 = fixpath("/work/data/23-6/N10000-alignsimple-t0.5-sn0.01.txt")
    test_interaction_file_30 = fixpath("/work/data/23-6/N10000-alignsimple-t0.5-sn0.01_0.txt")
    test_interaction_file_10 = fixpath("/work/data/23-6/N10000-alignsimple-10-t0.5-sn0.01.txt")
    test_no_interaction_file = fixpath("/work/data/22-6/N10000-nointeraction-t1-sn0.01.txt")

    interactionfliprates = [0, 10, 30, 60, 150, 300]
    legendlabels = string.(interactionfliprates)

    data_files = [test_no_interaction_file, test_interaction_file_10, test_interaction_file_30, test_interaction_file_60, test_interaction_file_150, test_interaction_file_300]

    for file in data_files
        simdata = loadsim(file, rowwisetxt)
    
        # analyze correlation
        ocdmat = orientationcorrelation(simdata, 0.25)
    
        # println(ocdmat)
    
        # handle plotting
        if LOG
            plt.plot(ocdmat[1,1:25], log.(ocdmat[2,1:25]))
        else
            plt.plot(ocdmat[1,1:25], ocdmat[2,1:25])
        end
    end

    plt.title("Orientation Correlatin\n$(nparticles) Particles\ntmax = 0.5")
    plt.xlabel("dt")
    plt.ylabel("Orientation Correlation")
    plt.legend(legendlabels)

    if LOG
        plt.ylim([-1.5,0])
    else
        plt.ylim([0,1])
    end

    plotfilepath::String = ""
    if LOG
        plotfilepath = fixpath("/work/data/correlation_data/N10000-ocf-logplot.pdf")
    else
        plotfilepath = fixpath("/work/data/correlation_data/N10000-ocf.pdf")
    end

    # plt.savefig(plotfilepath, bbox_inches = "tight",pad_inches=0.01)
    # plt.show()


end



ocf_test()

