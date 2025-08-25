########################
# paths.py
# 
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)
# 
# File containing helper methods to get absolute paths from relative paths on generic device
########################

import os
import sys

ROOT_DIR_ABS_PATH1 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev"
ROOT_DIR_ABS_PATH2 = "/home/s2696227/Documents/Dissertation/1DActiveMatterSim"

ROOT_DIRNAME_BD = "Dev"
ROOT_DIRNAME = "1DActiveMatterSim"


def get_root_abspath() -> str:
    """
    Gets the absolute path of the root of the 1DActiveMatterSim Repository on your machine.
    """
    currdir = os.getcwd()

    # look at parent directory until at the root
    while not (currdir.endswith(ROOT_DIRNAME) or currdir.endswith(ROOT_DIRNAME_BD)):
        currdir = os.path.dirname(currdir)
    
    return currdir


def relpath2abspath(path:str) -> str:
    """
    Gets the absolute path on your machine to a file given as a relative path to the root directory `1DActiveMatterSim`
    """
    return os.path.join(get_root_abspath(), path)




def fix_path(path:str) -> str:
    """
    Takes in an arbitrary path (absolute or relative) and returns the absolute version
    """
    # check if its an absolute or relative path
    abspath_start = os.getcwd()[:7]

    # if it's an absolute path, all good
    if path.startswith(abspath_start):
        return path
    else:
        # if not, append it to the absolute of the root path
        return relpath2abspath(path)
    

def abspath2relpath(path:str) -> str:
    """
    gets the relative path from the root of an absolute path
    """
    return os.path.relpath(path, get_root_abspath())




