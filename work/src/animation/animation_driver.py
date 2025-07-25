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
from file_management.simdata_files import loadsim, loadsim_n_lines, prepare_simfile
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


file_align = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01.txt")
save_filepath = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01.mp4")
save_filepath_faster = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01_faster.mp4")
particle_density_animate(file_align, 100, False, True, fps=40, delete_gif=True)


clip = moviepy.VideoFileClip(save_filepath)
faster_clip = resize_videofile(clip, 60, 30, save_filepath_faster)




## Generate redblue animation
# sim_animate(noint11, SHOW, SAVE, fps=100, y_offset=True, delete_gif=True)
# sim_animate(new_rowwise_short, SHOW, SAVE, fps=600, y_offset=True)
# sim_animate(new_rowwise_sn, SHOW, SAVE, fps=600, y_offset=True)


# data_dir = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/8-7/"
# density_sweep_dir = fix_path("work/data/17-7")

# datafile = fix_path("work/data/12-7/N1000-B100-alignsimple-300-t100-sn0.01_0.txt")


# og = fix_path("work/data/26-6/N1000-B100-alignsimple-300-t100-sn0.01_TEMP_COLLAPSED_1.txt")
# prepared_f = prepare_simfile(og)
# sim_animate_follow(prepared_f, 0, show=False, save=True, fps=5, y_offset=True, delete_gif=True)


# sim_animate(datafile, show=False, save=True, fps=100, y_offset=True, delete_gif=True)
# make_mp4s_of_dir(density_sweep_dir)
