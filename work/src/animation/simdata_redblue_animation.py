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

SPIN_UP_COLOUR = "blue"
SPIN_DOWN_COLOUR = "red"
ALPHA = 0.3
PLOT_YLIM = 0.2
PLOT_DOT_SIZE = 300

DEBUG_MODE_SHORT_NFRAMES = 5

SAVE = True
SHOW = False



def print_save_progress(current_frame: int, total_frames: int):
    pct_frames = max(1, total_frames // 100)
    if current_frame % pct_frames == 0:
        print(f"\033[KRendering Progress: {np.round(100 * current_frame / total_frames, 2)} %", end="\r")
    if (total_frames - current_frame) / total_frames < .3: 
        print(f"\033[KSaving File...", end="\r")

def sim_animate(file_str:str, show:bool = True, save:bool = False, fps:float = 30, y_offset=False, save_filepath=None, debug_mode=None, delete_gif=False):    
    
    # load data from file
    if debug_mode is None:
        sim_data:SimulationData = loadsim(file_str)
    elif debug_mode == "short":
        # load only a few frames for quick debug
        sim_data:SimulationData = loadsim_n_lines(file_str, 0, DEBUG_MODE_SHORT_NFRAMES, change_simparams=True)

    # calculate positions from data

    # Create a figure and axes
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ydata = np.zeros((1,len(xdata)))
    sc, = ax.plot([], [])

    times = sim_data._sim_params.get_save_times()

    # initialise xdata and ydata to al zeros
    offsets = np.zeros([sim_data._sim_params._num_particles, 2])

    init_xdata = offsets[:,0]
    init_ydata = offsets[:,1]

    # If option turned on, spread dots out along y-axis to make the particles easier to differentiate
    if y_offset:
        init_ydata = 2 * PLOT_YLIM * ((np.arange(0, sim_data._sim_params._num_particles) / sim_data._sim_params._num_particles)  - (1/2))
        offsets[:,1] = init_ydata

        # marker_size = int(np.round((2000 / 3) * PLOT_DOT_SIZE / sim_data._sim_params._num_particles))
        marker_size = max(10, np.pow((2 * PLOT_YLIM) / sim_data._sim_params._num_particles, 2))
    else:
        marker_size = PLOT_DOT_SIZE


    up_rgba = mplc.colorConverter.to_rgba(SPIN_UP_COLOUR, alpha=ALPHA)
    down_rgbpa = mplc.colorConverter.to_rgba(SPIN_DOWN_COLOUR, alpha=ALPHA)

    sc = ax.scatter(np.transpose(init_xdata), np.transpose(init_ydata), s=marker_size, c="k", marker='.', animated=True)
    title = ax.text((sim_data._sim_params._box_width/2) - .15, PLOT_YLIM+.01, "t = 0", fontsize=12)
    label = ax.text(-.5, PLOT_YLIM+.01, f"Particles: {sim_data._sim_params._num_particles}\nStochastic Flip Rate: {np.round(sim_data._sim_params._flip_rate, 4)}\nInteraction: {sim_data._sim_params._interaction}\nInteraction Flip Rate: {np.round(sim_data._sim_params._interaction_fliprate, 4)}", fontsize=8)

    def sim_init():
        ax.set_xlim(-sim_data._sim_params._box_width/2, sim_data._sim_params._box_width/2)
        ax.set_ylim(-PLOT_YLIM, PLOT_YLIM)
        ax.get_yaxis().set_visible(False)
        ax.set_xlabel("Particle Position")
        return sc, title

    def sim_update(frame):
        offsets[:,0] = sim_data.positions[frame]
        sc.set_offsets(offsets)
        sc.set_color([up_rgba if x == 1 else down_rgbpa for x in sim_data.spins[frame]])
        title.set_text(f"t = {np.round(times[frame], 5)}")

        return sc, title

    # compress animation and save
    anim_length_s = sim_data._sim_params._total_time
    anim_fps = fps
    max_length = 10
    # max_frames = sim_data._sim_params.get_ntimes()
    # max_frames = max_length * anim_fps
    anim_save_step = max(1, int(np.round(sim_data._sim_params.get_nsaves() / anim_length_s) / anim_fps))
    ani = FuncAnimation(fig, sim_update, frames=np.arange(0, sim_data._sim_params.get_nsaves(), anim_save_step, dtype=int),
                        init_func=sim_init, blit=True, interval=1000/anim_fps)
    
    if save:
        if save_filepath is None:
            in_file_prefix = file_str[:-4]
            save_filepath = f"{in_file_prefix}.mp4"
        if save_filepath.endswith(".mp4"):
            save_filepath_pre = save_filepath[:-4]
            save_filepath_gif = f"{save_filepath_pre}.gif"
        elif save_filepath.endswith(".gif"):
            save_filepath_gif = save_filepath
        
        ani.save(save_filepath_gif, fps=anim_fps, progress_callback=print_save_progress)

        if save_filepath.endswith(".mp4"):
            clip = moviepy.VideoFileClip(save_filepath_gif)
            clip.write_videofile(save_filepath)

            if delete_gif:
                os.remove(save_filepath_gif)



    if show:
        plt.show()



def make_mp4s_of_dir(dir_path:str):
    """
    Saves all txt sim data files in a directory as .mp4 video files
    """

    file_list = os.listdir(dir_path)
    for file_name in file_list:
        if file_name.endswith(".txt"):
            full_file_path = os.path.join(dir_path, file_name)
            print(f"Animating and Saving {full_file_path}...")
            sim_animate(full_file_path, show=False, save=True, fps=100, y_offset=True, delete_gif=True)




noint11 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/22-6/N10000-nointeraction-t1-sn0.01.txt"
# sim_animate(noint11, SHOW, SAVE, fps=100, y_offset=True, delete_gif=True)
# sim_animate(new_rowwise_short, SHOW, SAVE, fps=600, y_offset=True)
# sim_animate(new_rowwise_sn, SHOW, SAVE, fps=600, y_offset=True)
# sim_animate(example_align_intearaction_sim, SHOW, SAVE)


data_dir = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/23-6"

make_mp4s_of_dir(data_dir)
