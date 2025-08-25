########################
# script_paths.py
# 
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)
# 
# File containing constants of paths of helpful script files.
########################

import os
import sys

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.paths import fix_path

FILE_MANAGEMENT_SCRIPTS_DIR = fix_path("work/src/file_management/scripts")


CONVERT_SIMDATA = os.path.join(FILE_MANAGEMENT_SCRIPTS_DIR, "convertsimdata.jl")


