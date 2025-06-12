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
from file_management.simdata_files import loadsim

SPIN_UP_COLOUR = "blue"
SPIN_DOWN_COLOUR = "red"
ALPHA = 0.3

SAVE = False
SHOW = True
INTERACTION = "None"
save_filepath = "N100_strong_interaction_norandflip.mp4"

N3_rand = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N3-alignsimple_rand_simdata.txt"
N100_rand = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N100-alignsimple_rand_simdata.txt"
N100_rand_strong = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N100-alignsimple-300_simdata.txt"
N100_rand_strong_norandflip = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N100-alignsimple-300-noflip_simdata.txt"

example_no_intearaction_sim100 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N100-align-simplelong_simdata.txt"
example_no_intearaction_sim2 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N2-align-simplelong_simdata.txt"
example_no_intearaction_sim3 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N3-nointeraction_simdata.txt"
example_no_intearaction_sim3_int = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N3-alignsimple_simdata.txt"
example_no_intearaction_sim3_strong_int = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N3-strong-alignsimple_simdata.txt"
example_no_intearaction_sim3_mid_int = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N3-alignsimple_1_simdata.txt"
example_no_intearaction_sim100_mid_int = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/6-6/5-6-N100-alignsimple_simdata.txt"
example_no_intearaction_sim = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/basic_N100_t100000.0_interaction_none_T0.3333333333333333_sim_simdata.txt"
example_align_intearaction_sim = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/basic_N100_t100000.0_interaction_alignsimple_T0.3333333333333333_sim_simdata.txt"
example_align_intearaction_sim2 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/16-5/16-5-N100-align-simplelong_simdata.txt"

new_rowwise = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/saveas_rwt.txt_simdata.txt"


def print_save_progress(current_frame: int, total_frames: int):
    pct_frames = total_frames // 100
    if current_frame % pct_frames == 0:
        print(f"\033[KRendering Progress: {np.round(100 * current_frame / total_frames, 2)} %", end="\r")
    if (total_frames - current_frame) / total_frames < .3: 
        print(f"\033[KSaving File...", end="\r")

def sim_animate(file_str:str, show:bool = True, save:bool = False, fps:float = 30):    
    
    # load data from file
    sim_data:SimulationData = loadsim(file_str)

    # calculate positions from data

    # Create a figure and axes
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ydata = np.zeros((1,len(xdata)))
    sc, = ax.plot([], [])

    times = sim_data._sim_params.get_times()
    print(times)

    # initialise xdata and ydata to al zeros
    offsets = np.zeros([sim_data._sim_params._num_particles, 2])

    init_xdata = offsets[:,0]
    init_ydata = offsets[:,1]

    up_rgba = mplc.colorConverter.to_rgba(SPIN_UP_COLOUR, alpha=ALPHA)
    down_rgbpa = mplc.colorConverter.to_rgba(SPIN_DOWN_COLOUR, alpha=ALPHA)

    sc = ax.scatter(np.transpose(init_xdata), np.transpose(init_ydata), s=300, c="k", marker='.', animated=True)
    title = ax.text(0, .15, "t = 0")
    label = ax.text(-.5, .15, f"Particles: {sim_data._sim_params._num_particles}\nStochastic Flip Rate: {np.round(sim_data._sim_params._flip_rate, 4)}\nInteraction: {sim_data._sim_params.interaction}\nInteraction Flip Rate: {np.round(sim_data._sim_params.interaction_fliprate, 4)}")

    def sim_init():
        ax.set_xlim(-sim_data._sim_params._box_width/2, sim_data._sim_params._box_width/2)
        ax.set_ylim(-.2, .2)
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
    anim_save_step = max(1, int(np.round(sim_data._sim_params.get_ntimes() / anim_length_s) / anim_fps))

    ani = FuncAnimation(fig, sim_update, frames=np.arange(0, sim_data._sim_params.get_ntimes(), anim_save_step, dtype=int),
                        init_func=sim_init, blit=True, interval=1000/anim_fps)
    
    if save:
        if save_filepath.endswith(".mp4"):
            save_filepath_pre = save_filepath[:-4]
            save_filepath_gif = f"{save_filepath_pre}.gif"
        elif save_filepath.endswith(".gif"):
            save_filepath_gif = save_filepath
        
        ani.save(save_filepath_gif, fps=anim_fps, progress_callback=print_save_progress)

        if save_filepath.endswith(".mp4"):
            clip = moviepy.VideoFileClip(save_filepath_gif)
            clip.write_videofile(save_filepath)



    if show:
        plt.show()



sim_animate(new_rowwise, SHOW, SAVE)
# sim_animate(example_align_intearaction_sim, SHOW, SAVE)