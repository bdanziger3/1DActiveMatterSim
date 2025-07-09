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
from file_management.simdata_files import loadsim, loadsim_n_lines
from animation.video_file import resize_videofile


clip_filepath = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/8-7/N1000-B100.0-nointeraction-100-T100_0.mp4"
new_clip_filepath = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/8-7/N1000-B100.0-nointeraction-100-T100_0_s.mp4"
clip = moviepy.VideoFileClip(clip_filepath)
resize_videofile(clip, 30, 30, new_clip_filepath)