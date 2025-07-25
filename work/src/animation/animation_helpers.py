########################
# animation_helpers.py
# Blake Danziger
# 1D Active Solids
# MSc Theoretical Physics Dissertation (2025)

# Helper Functions for plot animation code
########################


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.animation import FuncAnimation
import moviepy

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.sim_structs import SimulationData, SimulationParameters
from file_management.simdata_files import loadsim, loadsim_n_lines, get_simfile_type, prepare_simfile
from utils.paths import fix_path


def print_save_progress(current_frame: int, total_frames: int):
    pct_frames = max(1, total_frames // 100)
    if current_frame % pct_frames == 0:
        print(f"\033[KRendering Progress: {np.round(100 * current_frame / total_frames, 2)} %", end="\r")
    if (total_frames - current_frame) / total_frames < .01 or (total_frames - current_frame) <= 2: 
        print(f"\033[KSaving File...", end="\r")
