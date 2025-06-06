import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.animation import FuncAnimation
from work.simulation.sim_structs import SimulationData, SimulationParameters
from simdata_files import loadsim

SPIN_UP_COLOUR = "blue"
SPIN_DOWN_COLOUR = "red"
ALPHA = 0.3

SAVE = False
SHOW = True
COMPRESS = True
INTERACTION = "None"
save_filepath = "ee_test_out.gif"

# file_str = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/1basic_N100_t10000.0_align_T0.3333333333333333_sim_simdata.txt"
file_str = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/1basic_N100_t100000.0_align_T0.3333333333333333_sim_simdata.txt"
# file_str = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/1basic_N100_t100000.0_noint_T0.3333333333333333_sim_simdata.txt"
# file_str = loadsim("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/basic_N10_t1000.0_align_T1_sim_simdata.txt")
# sim_data:SimulationData = loadsim("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/basic_N10_t1000.0_align_T1_sim_simdata.txt")
sim_data:SimulationData = loadsim(file_str)


def print_save_progress(current_frame: int, total_frames: int):
    pct_frames = total_frames // 100
    if current_frame % pct_frames == 0:
        print(f"\033[KRendering Progress: {np.round(100 * current_frame / total_frames, 2)} %", end="\r")
    if (total_frames - current_frame) / total_frames < .3: 
        print(f"\033[KSaving File...", end="\r")

def sim_animate(sim_data:SimulationData):
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
    label = ax.text(-.5, .15, f"Particles: {sim_data._sim_params._num_particles}\nFlip Rate: {sim_data._sim_params._flip_rate}\nInteraction: {INTERACTION}")

    def sim_init():
        ax.set_xlim(-sim_data._sim_params._box_width/2, sim_data._sim_params._box_width/2)
        ax.set_ylim(-.2, .2)
        ax.get_yaxis().set_visible(False)
        ax.set_xlabel("Particle Position")
        return sc, title

    def sim_update(frame):
        offsets[:,0] = sim_data.wrapped_positions[frame]
        sc.set_offsets(offsets)
        sc.set_color([up_rgba if x == 1 else down_rgbpa for x in sim_data.spins[frame]])
        title.set_text(f"t = {np.round(times[frame], 5)}")

        return sc, title

    # compress animation and save
    if COMPRESS:
        anim_length_s = 100
        anim_fps = 30
        anim_save_step = max(1, int(np.round(sim_data._sim_params.get_ntimes() / anim_length_s) / anim_fps))
    else:
        anim_fps = 1 / sim_data._sim_params._dt
        anim_save_step = 1

    ani = FuncAnimation(fig, sim_update, frames=np.arange(0, sim_data._sim_params.get_ntimes(), anim_save_step, dtype=int),
                        init_func=sim_init, blit=True, interval=1000/anim_fps)
    
    if SAVE:
        ani.save(save_filepath, fps=anim_fps, progress_callback=print_save_progress)

    if SHOW:
        plt.show()



sim_animate(sim_data)