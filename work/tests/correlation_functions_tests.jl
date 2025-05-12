using PyPlot
using PyCall

include("../corelation_functions.jl")
include("../simdata_files.jl")

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
    # open simdata
    
    # simplesim_datafile = "../txt_data_files/basic_N10_t1000.0_noint_T1_sim_simdata.txt"
    # plotfilepath = "../plot_outputs/Apr25/OCF_N10_t100_T1_noint.pdf"
    
    # simplesim_datafile = "../txt_data_files/basic_N100_t1000_noint_sim_simdata.txt"
    # plotfilepath = "../plot_outputs/Apr25/OCF_N100_t1000_T1_noint.pdf"
    
    # simplesim_datafile = "../txt_data_files/basic_N1000_t10000_noint_sim_simdata.txt"
    # plotfilepath = "../plot_outputs/Apr25/OCF_N1000_t10000_T1_noint.pdf"

    simdata = loadsim(simplesim_datafile, sequentialtxt)

    # analyze correlation
    ocdmat = ocf(simdata)

    # handle plotting

    


    plt.plot(ocdmat[1,2:end], ocdmat[2,2:end])
    plt.title("Orientation Correlatin\n10 Particles\ntmax = 10")
    plt.xlabel("dt")
    plt.ylabel("Orientation Correlation")
    plt.ylim([-1,1])


    plt.savefig(plotfilepath, bbox_inches = "tight",pad_inches=0.01)
    # plt.show()


end



ocf_test()

