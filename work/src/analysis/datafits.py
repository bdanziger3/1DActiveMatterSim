########################
# datafits.py
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# File containing functions to fit analysis data to curves
########################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as scp
from scipy.optimize import curve_fit 

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from file_management.simdata_files import string_list_to_number_list
from utils.paths import fix_path

def fit_oc_decay(datafile:str, plot_data:bool = False) -> list[float]:
    # open file
    f = open(datafile, "r")

    line = f.readline()

    while "sweep values" not in line.lower():
        line = f.readline()
    
    # read sweep values
    line = f.readline()
    sweep_values = string_list_to_number_list(line, ",")

    while "data matrix" not in line.lower():
        line = f.readline()

    decay_params = []
    # read t-data
    line = f.readline()
    tdata = string_list_to_number_list(line, ",")

    # read oc-data and fit to exponential OC = exp(a * t)
    max_index = 1000
    line = f.readline()
    curve_to_fit = lambda t,a: np.exp(- t / a)
    while line != "":
        oc_data = string_list_to_number_list(line, ",")

        popt, pcov = curve_fit(curve_to_fit, tdata[:max_index], oc_data[:max_index])
        decay_params.append(popt[0])

        line = f.readline()

    f.close()

    print(decay_params)

    if plot_data:
        plt.plot(sweep_values, decay_params, "-o")
        plt.show()

    return decay_params


def plot_oc_decay(sweep_values, decay_parameters):


    plt.plot(sweep_values, decay_parameters)
    plt.show()



file = fix_path("work/analysis/Orienation Self-Correlation/Align Simple/Interaction Sweep N1000/orientation_selfcorrelation_plot-settletime-1.0_0.txt")

fit_oc_decay(file, True)
