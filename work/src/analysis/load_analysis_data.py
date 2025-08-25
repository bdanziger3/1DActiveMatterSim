########################
# load_analysis_data.py
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# File containing functions to load output .txt data from analysis functions
########################


import sys
import os
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
import scipy as scp
from scipy.optimize import curve_fit 

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.sim_structs import SimulationParameters
from file_management.simdata_files import string_list_to_number_list
from utils.paths import fix_path

class AnalysisType(Enum):
    MSD = 1
    ORIENTATION_CORRELATION = 2
    POLAR_ORDER = 3
    MEAN_SPIN = 4

class SweepType(Enum):
    NONE = 0
    DENSITY_SWEEP = 1
    BOX_WIDTH_SWEEP = 2
    INTERACTION_STRENGTH_SWEEP = 3


def parse_analysis_type(datafile_header:str) -> AnalysisType:
    if "mean squared displacement data" in datafile_header.lower():
        return AnalysisType.MSD
    elif "orientation self-correlation data" in datafile_header.lower():
        return AnalysisType.ORIENTATION_CORRELATION
    elif "polar order data" in datafile_header.lower():
        return AnalysisType.POLAR_ORDER
    elif "mean spin data" in datafile_header.lower():
        return AnalysisType.MEAN_SPIN
    else:
        assert False, f"Could not parse header {datafile_header}."

def parse_sweep_type(sweeptype_line:str) -> SweepType:
    n_chars_to_skip = len("Sweep Type: ")
    type_str = sweeptype_line[n_chars_to_skip:]
    if "densitysweep" in type_str.lower():
        return SweepType.DENSITY_SWEEP
    elif "boxwidthsweep" in type_str.lower():
        return SweepType.BOX_WIDTH_SWEEP
    elif "interactionstrengthsweep" in type_str.lower():
        return SweepType.INTERACTION_STRENGTH_SWEEP
    else:
        print(f"Could not parse header {sweeptype_line}.")
        return SweepType.NONE



def load_analysis_data(data_path:str):
    """
    Gets the data matrix of a given analysis data file.

    Returns 4 objects:
    `(xy_data_matrix, analysis_type, sweep_type, sim_params_out)`

    `xy_data_matrix`:   The matrix holding the analysis data. Row 0 holds the x-axis data.
                        The remaining lines hold the y-axis data
                        
    `analysis_type`:    The `AnalysisType` of the data

    `sweep_type`:       The `SweepType` of the data. Is `SweepType.NONE` if not a sweep.
    
    `sim_params_out`:   The simulation parameters used. Returns a `SimulationParameters` object if not a sweep.
                        Returns a string to use in the plot title if it is a sweep
    """

    # open file
    f = open(data_path, "r")

    line = f.readline()

    # get the type of data
    analysis_type = parse_analysis_type(line)

    # parse sweep type
    line = f.readline()
    sweep_type = parse_sweep_type(line)

    if sweep_type == SweepType.NONE:
        # if not a sweep, read simulation params and skip to data matrix
        sweep_values = []
        
        # read next line to get number of position data points
        sim_param_info = f.readline()
        sim_param_info_list = string_list_to_number_list(sim_param_info, out_type=SimulationParameters.expected_types)
        sim_params = SimulationParameters.construct_from_list(sim_param_info_list)

        sim_params_out = sim_params

    else:
        # get sweep values if a sweep
        line = f.readline()
        sweep_values = string_list_to_number_list(line)

        # load param string
        line = f.readline()
        params_str = line[len("Params String: "):]

        sim_params_out = params_str

    
    # skip to data matrix
    while "data matrix" not in line.lower() and "plot data" not in line.lower():
        line = f.readline()


    # read x-data
    x_data_line = f.readline()
    x_data_values = string_list_to_number_list(x_data_line)

    # find number of lines
    y_data_lines = len(sweep_values)

    # read in lines
    xy_data_matrix = np.zeros(y_data_lines+1, len(x_data_values))

    xy_data_matrix[0, :] = x_data_values

    for i in y_data_lines:
        y_data_line_i = f.readline()
        xy_data_matrix[i+1, :] = string_list_to_number_list(y_data_line_i)


    return xy_data_matrix, analysis_type, sweep_type, sim_params_out










