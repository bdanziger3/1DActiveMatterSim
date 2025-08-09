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

SPIN_UP_COLOR = "red"
SPIN_DOWN_COLOR = "blue"
FOLLOWING_COLOR = "darkorchid"
ALPHA = 0.3
PLOT_YLIM = 0.2
PLOT_DOT_SIZE = 10

DEBUG_MODE_SHORT_NFRAMES = 5

SAVE = True
SHOW = False



def print_save_progress(current_frame: int, total_frames: int):
    pct_frames = max(1, total_frames // 100)
    if current_frame % pct_frames == 0:
        print(f"\033[KRendering Progress: {np.round(100 * current_frame / total_frames, 2)} %", end="\r")
    if (total_frames - current_frame) / total_frames < .01 or (total_frames - current_frame) <= 2: 
        print(f"\033[KSaving File...", end="\r")

def sim_animate(file_str:str, show:bool = True, save:bool = False, fps:float = 30, y_offset=False, save_filepath=None, debug_mode=None, file_suffix=None, delete_uncompressed=False):
    # prepare simfile if needed
    prepared_file_path = prepare_simfile(file_str)
    
    if file_suffix is None:
        file_suffix = ""
    
    # load data from file
    if debug_mode is None:
        sim_data:SimulationData = loadsim(prepared_file_path)
    elif debug_mode == "short":
        # load only a few frames for quick debug
        sim_data:SimulationData = loadsim_n_lines(prepared_file_path, 0, DEBUG_MODE_SHORT_NFRAMES, change_simparams=True)

    # calculate positions from data

    # Create a figure and axes
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ydata = np.zeros((1,len(xdata)))
    sc, = ax.plot([], [])

    times = sim_data._sim_params.get_save_times()

    # initialize xdata and ydata to all zeros
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
        plot_height_fraction = .1
        init_ydata = plot_height_fraction * 2 * PLOT_YLIM * ((np.arange(0, sim_data._sim_params._num_particles) / sim_data._sim_params._num_particles)  - (1/2))
        offsets[:,1] = init_ydata
        marker_size = max(10, np.pow((2 * plot_height_fraction * PLOT_YLIM) / sim_data._sim_params._num_particles, 2))

    up_rgba = mplc.colorConverter.to_rgba(SPIN_UP_COLOR, alpha=ALPHA)
    down_rgbpa = mplc.colorConverter.to_rgba(SPIN_DOWN_COLOR, alpha=ALPHA)

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
        
        ani.save(save_filepath, fps=anim_fps, progress_callback=print_save_progress)


    if show:
        plt.show()

    # delete temporary uncompressed file
    if delete_uncompressed:
        if prepared_file_path != file_str:
            os.remove(prepared_file_path)


def sim_animate_follow(file_str:str, particle_to_follow:int, show:bool = True, save:bool = False, fps:float = 30, y_offset=False, save_filepath=None, debug_mode=None, delete_gif=False):    
    
    # prepare simfile if needed
    prepared_file_path = prepare_simfile(file_str)
    
    # load data from file
    if debug_mode is None:
        sim_data:SimulationData = loadsim(prepared_file_path)
    elif debug_mode == "short":
        # load only a few frames for quick debug
        sim_data:SimulationData = loadsim_n_lines(prepared_file_path, 0, DEBUG_MODE_SHORT_NFRAMES, change_simparams=True)

    # calculate positions from data

    # Create a figure and axes
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ydata = np.zeros((1,len(xdata)))
    sc, = ax.plot([], [])

    times = sim_data._sim_params.get_save_times()

    # initialize xdata and ydata to al zeros
    offsets = np.zeros([sim_data._sim_params._num_particles, 2])

    init_xdata = offsets[:,0]
    init_ydata = offsets[:,1]

    # If option turned on, spread dots out along y-axis to make the particles easier to differentiate
    if y_offset:
        init_ydata = 2 * PLOT_YLIM * ((np.arange(0, sim_data._sim_params._num_particles) / sim_data._sim_params._num_particles)  - (1/2))
        init_ydata[particle_to_follow] = 0
        init_ydata[sim_data._sim_params._num_particles // 2] = -PLOT_YLIM
        offsets[:,1] = init_ydata

        # marker_size = int(np.round((2000 / 3) * PLOT_DOT_SIZE / sim_data._sim_params._num_particles))
        marker_size = max(10, np.pow((2 * PLOT_YLIM) / sim_data._sim_params._num_particles, 2))
    else:
        marker_size = PLOT_DOT_SIZE

    # make marker_size array

    marker_size = np.array([marker_size] * sim_data._sim_params._num_particles)
    marker_size[particle_to_follow] = marker_size[0] * 2

    if y_offset:
        marker_size[particle_to_follow] = marker_size[0] * 5



    up_rgba = mplc.colorConverter.to_rgba(SPIN_UP_COLOR, alpha=ALPHA)
    down_rgbpa = mplc.colorConverter.to_rgba(SPIN_DOWN_COLOR, alpha=ALPHA)
    following_rgbpa = mplc.colorConverter.to_rgba(FOLLOWING_COLOR, alpha=ALPHA)

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
        color_array = [up_rgba if x == 1 else down_rgbpa for x in sim_data.spins[frame]]
        color_array[particle_to_follow] = following_rgbpa
        sc.set_color(color_array)
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
            save_filepath = f"{in_file_prefix}_follow_{particle_to_follow}.mp4"
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



def make_mp4s_of_dir(dir_path:str, only_prepared_files:bool=False, fps:float=30, clear_new_files:bool=True, file_suffix=None, y_offset=False):
    """
    Saves all txt sim data files in a directory as .mp4 video files

    Set `only_prepared_files=False` to only generate an animation of data files that are already 1 segment and deserialized.
    If it is `True`, then prepare the data file and generate the animation.
    """

    if file_suffix is None:
        file_suffix = ""

    file_list = os.listdir(dir_path)
    for file_name in file_list:
        if file_name.endswith(".txt"):

            # if `only_prepared_files`, skip files that aren't collapsed or deserialized
            full_file_path = os.path.join(dir_path, file_name)
            (n_segments, is_compressed) = get_simfile_type(full_file_path)

            if only_prepared_files and (n_segments != 1 or is_compressed):
                continue
            else:
                prepared_file_path = prepare_simfile(full_file_path)

            print(f"Animating and Saving {prepared_file_path}...")
            
            try:
                sim_animate(prepared_file_path, show=False, save=True, fps=fps, y_offset=y_offset, delete_gif=True, file_suffix=file_suffix)

                # delete prepared file if new
                if clear_new_files and prepared_file_path != full_file_path:
                    os.remove(prepared_file_path)
            except:
                print(f"Could not animate sim at {file_name}.")

