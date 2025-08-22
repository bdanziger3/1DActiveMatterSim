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
from animation.video_file import resize_videofile
from animation.simdata_redblue_animation import sim_animate, make_mp4s_of_dir, sim_animate_follow
from animation.density_hist_animation import particle_density_animate
from utils.paths import fix_path



# # Concatenate video files
# clip_filepath_1 = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01_TEMP_COLLAPSED_1_follow_first50s.mp4")
# clip_filepath_2 = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01_TEMP_COLLAPSED_1_follow_50-300s.mp4")
# clip_filepath_3 = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01_TEMP_COLLAPSED_1_follow-300-400s.mp4")
# clip_filepath_4 = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01_TEMP_COLLAPSED_1_follow_400-1000s.mp4")
# new_clip_filepath = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t1000-sn0.01_long.mp4")
# newer_clip_filepath = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t1000-sn0.01_faster.mp4")
# clip1 = moviepy.VideoFileClip(clip_filepath_1).subclipped(0, 50)
# clip2 = moviepy.VideoFileClip(clip_filepath_2).subclipped(0, 250)
# clip3 = moviepy.VideoFileClip(clip_filepath_3)
# clip4 = moviepy.VideoFileClip(clip_filepath_4)

# concatenate clips and write the result
# final_clip = moviepy.concatenate_videoclips([clip1, clip2, clip3, clip4])
# final_clip.write_videofile(new_clip_filepath)


for i in [1, 2, 5, 10, 20, 50, 100]:
    file1 = fix_path(f"work/data/sweeps/turnaway/interactionsweep/N500-B50/interactionsweep/N500-B50.0-turnaway-{i}0.0_rep1.txt")
    sim_animate(file1, False, True, 30, False, delete_uncompressed=True, color_track=True)

# file1 = fix_path("work/data/sweeps/turnaway/interactionsweep/N500-B50/interactionsweep/N500-B50.0-turnaway-2.0_rep3.txt")
# sim_animate(file1, False, True, 30, False, delete_uncompressed=True, color_track=True)

# make_mp4s_of_dir(file1, y_offset=True)