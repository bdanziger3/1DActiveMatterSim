import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.animation import FuncAnimation
from sim_structs import SimulationData, SimulationParameters
from simdata_files import loadsim, DataFileType
from utils.print_tools import ProgressBar
import time


SAVE = False
SHOW = True
INTERACTION = "None"
save_filepath = "no_interaction.gif"

example_no_intearaction_sim2 = "work/16-5/16-5-N100-no-interaction-2_simdata.txt"
example_no_intearaction_sim = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/basic_N100_t100000.0_interaction_none_T0.3333333333333333_sim_simdata.txt"
example_align_intearaction_sim = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/basic_N100_t100000.0_interaction_alignsimple_T0.3333333333333333_sim_simdata.txt"
example_align_intearaction_sim2 = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/16-5/16-5-N100-align-simplelong_simdata.txt"


def simdata_unwrap(sim_data:SimulationData):
    # default behavior uses `jump_thresh` of half the box width
    if jump_thresh <= 0:
        jump_thresh = sim_data._sim_params._box_width / 2

    # return error if jump_thresh box width is too small
    particle_dx = sim_data._sim_params._v0 * sim_data._sim_params._dt
    if particle_dx >= jump_thresh:
        raise ValueError("The threshold to detect particle wrapping is smaller than the expected distance particles move between frames, making it hard to detect when particles wrap.")

    nframes = sim_data._sim_params.get_ntimes()

    pb = ProgressBar(nframes, "Calculating wraps...", "Done")

    nwraps = np.zeros_like(sim_data.wrapped_positions)

    # count which 'screen' the particles are on by counting the number of times they wrap around the box edges
    for i in range(1, nframes):
        pb.sparse_progress(i)
        # particle advances one screen for each large "jump" in the positions (in the opposite direction of the jump)
        dxs = sim_data.wrapped_positions[i] - sim_data.wrapped_positions[i-1]
        nwraps[i, :] = nwraps[i-1, :] + (-np.sign(dxs) * (np.abs(dxs) >= jump_thresh))

    unwrapped_positions = sim_data.wrapped_positions + (nwraps * sim_data._sim_params._box_width)

    return unwrapped_positions

def sim_unwrap_preload(file_str:str, jump_thresh:float = -1, save:bool = False):
    # load data file and unwrap the positions
    sim_data:SimulationData = loadsim(file_str, file_type=DataFileType.SEQUENTIAL_TEXT_ABS)
    unwrapped_positions = sim_unwrap_preload(sim_data)

    return unwrapped_positions
    


def sim_wrap_buffer(file_str:str, jump_thresh:float = -1, save:bool = False):
    
    sim_data:SimulationData = loadsim(file_str, file_type=DataFileType.SEQUENTIAL_TEXT_ABS)

    # default behavior uses `jump_thresh` of half the box width
    if jump_thresh <= 0:
        jump_thresh = sim_data._sim_params._box_width / 2

    # return error if jump_thresh box width is too small
    particle_dx = sim_data._sim_params._v0 * sim_data._sim_params._dt
    if particle_dx >= jump_thresh:
        raise ValueError("The threshold to detect particle wrapping is smaller than the expected distance particles move between frames, making it hard to detect when particles wrap.")

    nframes = sim_data._sim_params.get_ntimes()

    pb = ProgressBar(nframes, "Calculating wraps...", "Done")

    nwraps = np.zeros_like(sim_data.wrapped_positions)

    # count which 'screen' the particles are on by counting the number of times they wrap around the box edges
    for i in range(1, nframes):
        pb.sparse_progress(i)
        # particle advances one screen for each large "jump" in the positions (in the opposite direction of the jump)
        dxs = sim_data.wrapped_positions[i] - sim_data.wrapped_positions[i-1]
        nwraps[i, :] = nwraps[i-1, :] + (-np.sign(dxs) * (np.abs(dxs) >= jump_thresh))

    all_calc_pos = sim_data.wrapped_positions + (nwraps * sim_data._sim_params._box_width)

    return all_calc_pos



t0 = time.time()
sim_wrap_preload(example_align_intearaction_sim2, -1, SHOW, SAVE)
print(f"total time: {time.time() - t0}")
# sim_animate(example_align_intearaction_sim, SHOW, SAVE)