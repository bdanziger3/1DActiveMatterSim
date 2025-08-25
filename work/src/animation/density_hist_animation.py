########################
# density_hist_animation.py
# Blake Danziger
# 1D Active Matter Sim
# MSc Theoretical Physics Dissertation (2025)

# Functions to plot a density histogram of simulation data
########################


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import moviepy

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.sim_structs import SimulationData, SimulationParameters
from file_management.simdata_files import loadsim, loadsim_n_lines, get_simfile_type, prepare_simfile, DataFileType
from utils.paths import fix_path
from animation.animation_helpers import print_save_progress

SPIN_UP_COLOR = "blue"
SPIN_DOWN_COLOR = "red"
COLORMAP = "coolwarm"

DEBUG_MODE_SHORT_NFRAMES = 5

CMAP = cm.coolwarm_r

SAVE = True
SHOW = False



def get_patch_polarity(positions:list[float], spins:list[int], bin_left_edge:int, bin_righ_tedge:int) -> float:
    # find average spin spins of all particles within given bounds

    # get indices of all particles in the bin
    inbin = (positions >= bin_left_edge) * (positions <= bin_righ_tedge)
    spins_inbin = spins[inbin]

    if len(spins_inbin) == 0:
        return 0
    else:
        return sum(spins_inbin)  / len(spins_inbin)


def plot_particle_density(simdata:SimulationData, time_index:int, particle_density_nbins:int=-1):

    positions = simdata.positions[time_index, :]
    spins = simdata.spins[time_index, :]

    if particle_density_nbins == -1:
        particle_density_nbins = int(round(simdata._sim_params._box_width))

        # except if it's 1
        if particle_density_nbins == 1:
            particle_density_nbins = 100
        

    interaction_str = ""
    if simdata._sim_params._interaction == "nointeraction":
        interaction_str = "no-interaction  "
    else:
        interaction_str = f"{simdata._sim_params._interaction} I={int(round(simdata._sim_params._interaction_fliprate))}  "

    ### Plotting

    # clear plot
    plt.clf()

    fig, ax = plt.subplots()
    
    ax.grid(True, zorder=0)
    n, bins, patches = ax.hist(positions, bins=particle_density_nbins, edgecolor="black", linewidth=.2, zorder=3, density=False)
    plt.xlabel("Particle Position")
    plt.ylabel("Particle Density")
    # plt.title(raw"Occurence of Particle Densities")
    # plt.ylabel(raw"Log of Probability Density (Log$P(n)$) of Bin with Particle Density $n$")
    plt.title(f"Particle Densities\n N={simdata._sim_params._num_particles}  Boxwidth={simdata._sim_params._box_width}  t={int(simdata._sim_params._total_time)}  {interaction_str}")
    

    # Calculate polarities of bins
    for i in range(particle_density_nbins):
        bin_left_edge = bins[i]
        bin_righ_tedge = bins[i + 1]

        # get bin polarity on scale from [0, 1]
        bin_polarity = (get_patch_polarity(positions, spins, bin_left_edge, bin_righ_tedge) / 2) + .5

        bin_color = CMAP(bin_polarity)
        patches[i].set_facecolor(bin_color)


    # Create a ScalarMappable for the colorbar
    norm = mplc.Normalize(vmin=-1, vmax=1)
    sm = cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])  # required for older matplotlib

    # Add the colorbar
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Average Polarity")
    cbar.set_ticks([-1, 1])  # labels on the ends
    cbar.ax.set_yticklabels(["-1", "+1"])

    plt.show()


def particle_density_animate(file_str:str, particle_density_nbins:int = -1, show:bool = True, save:bool = False, fps:float = 30, save_filepath=None, debug_mode=None, delete_gif=False):    
    
    # load data from file
    if debug_mode is None:
        sim_data:SimulationData = loadsim(file_str)
    elif debug_mode == "short":
        # load only a few frames for quick debug
        sim_data:SimulationData = loadsim_n_lines(file_str, 0, DEBUG_MODE_SHORT_NFRAMES, change_simparams=True)

    if particle_density_nbins == -1:
        particle_density_nbins = int(round(sim_data._sim_params._box_width))

        # except if it's 1
        if particle_density_nbins == 1:
            particle_density_nbins = 100
    
    # calculate bins
    bin_edges = np.linspace(-sim_data._sim_params._box_width / 2, sim_data._sim_params._box_width / 2, particle_density_nbins + 1)

    interaction_str = ""
    if sim_data._sim_params._interaction == "nointeraction":
        interaction_str = "no-interaction  "
    else:
        interaction_str = f"{sim_data._sim_params._interaction} I={int(round(sim_data._sim_params._interaction_fliprate))}  "


    # initialize plot
    plot_ylim = sim_data._sim_params._num_particles / 10

    # clear plot
    fig, ax = plt.subplots()

    times = sim_data._sim_params.get_save_times()
    time_label = ax.text((sim_data._sim_params._box_width/2) * 1.15, 1.1*plot_ylim, "t = 0", fontsize=12)

    def add_plot_labels():
        # add labels and axes limits
        ax.set_xlim(-sim_data._sim_params._box_width / 2, sim_data._sim_params._box_width / 2)
        ax.set_ylim(0, plot_ylim)
        ax.set_xlabel("Particle Position")
        ax.set_ylabel("Particle Density")
        title = ax.set_title(f"Particle Densities\n N={sim_data._sim_params._num_particles}  Boxwidth={sim_data._sim_params._box_width}  t={int(sim_data._sim_params._total_time)}  {interaction_str}")

        # Create a ScalarMappable for the colorbar
        norm = mplc.Normalize(vmin=-1, vmax=1)
        sm = cm.ScalarMappable(cmap=CMAP, norm=norm)
        sm.set_array([])  # required for older matplotlib

        # Add the colorbar
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Average Polarity")
        cbar.set_ticks([-1, 1])  # labels on the ends
        cbar.ax.set_yticklabels(["-1", "+1"])


    def hist_anim_init():
        ax.cla()  # clear current frame
        ax.grid(True, zorder=0)
        data = sim_data.positions[0, :]
        n, bins, bars = ax.hist(data, bins=bin_edges, edgecolor="black", linewidth=.2, zorder=3, density=False)
        add_plot_labels()
        time_label = ax.text((sim_data._sim_params._box_width/2) * 1.15, 1.1*plot_ylim, "t = 0", fontsize=12)

        # Calculate polarities of bins
        for i in range(particle_density_nbins):
            bin_left_edge = bins[i]
            bin_righ_tedge = bins[i + 1]

            # get bin polarity on scale from [0, 1]
            bin_polarity = (get_patch_polarity(data, sim_data.spins[i, :], bin_left_edge, bin_righ_tedge) / 2) + .5

            bin_color = CMAP(bin_polarity)
            bars[i].set_facecolor(bin_color)
        
        return bars

    def hist_anim_update(frame):
        # update histogram
        ax.cla()  # clear current frame
        ax.grid(True, zorder=0)
        data = sim_data.positions[frame, :]
        n, bins, bars = ax.hist(data, bins=bin_edges, edgecolor="black", linewidth=.2, zorder=3, density=False)
        ax.set_xlim(-sim_data._sim_params._box_width / 2, sim_data._sim_params._box_width / 2)
        ax.set_ylim(0, plot_ylim)
        ax.set_xlabel("Particle Position")
        ax.set_ylabel("Particle Density")
        title = ax.set_title(f"Particle Densities\n N={sim_data._sim_params._num_particles}  Boxwidth={sim_data._sim_params._box_width}  t={int(sim_data._sim_params._total_time)}  {interaction_str}")
        # add_plot_labels()
    
        # update bar colors
        # Calculate polarities of bins
        for i in range(particle_density_nbins):
            bin_left_edge = bins[i]
            bin_righ_tedge = bins[i + 1]

            # get bin polarity on scale from [0, 1]
            bin_polarity = (get_patch_polarity(data, sim_data.spins[i, :], bin_left_edge, bin_righ_tedge) / 2) + .5

            bin_color = CMAP(bin_polarity)
            bars[i].set_facecolor(bin_color)

        # update label
        ax.text((sim_data._sim_params._box_width/2) * 1.15, 1.1*plot_ylim, f"t = {np.round(times[frame], 5)}", fontsize=12)

        return bins, time_label

    # compress animation and save
    anim_length_s = sim_data._sim_params._total_time
    anim_fps = fps
    max_length = 10
    anim_save_step = max(1, int(np.round(sim_data._sim_params.get_nsaves() / anim_length_s) / anim_fps))
    ani = FuncAnimation(fig, hist_anim_update, frames=np.arange(0, sim_data._sim_params.get_nsaves(), anim_save_step, dtype=int),
                        init_func=hist_anim_init, blit=False, interval=1000/anim_fps)
    # hist_anim_init()
    if save:
        if save_filepath is None:
            in_file_prefix = file_str[:-4]
            save_filepath = f"{in_file_prefix}_density_hist_polarization.mp4"
        if save_filepath.endswith(".mp4"):
            save_filepath_pre = save_filepath[:-4]
            save_filepath_gif = f"{save_filepath_pre}_density_hist_polarization.gif"
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



def make_mp4s_of_dir(dir_path:str, only_prepared_files:bool=False):
    """
    Saves all txt sim data files in a directory as .mp4 video files

    Set `only_prepared_files=False` to only generate an animation of data files that are already 1 segment and deserialized.
    If it is `True`, then prepare the data file and generate the animation.
    """

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
                particle_density_animate(prepared_file_path, show=False, save=True, fps=100, y_offset=True, delete_gif=True)
            except:
                print(f"Could not animate sim at {file_name}.")

