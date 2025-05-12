import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sim_structs import SimulationData, SimulationParameters
from simdata_files import loadsim

# file_str = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/1basic_N100_t10000.0_align_T0.3333333333333333_sim_simdata.txt"
# file_str = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/1basic_N100_t100000.0_align_T0.3333333333333333_sim_simdata.txt"
file_str = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/1basic_N100_t100000.0_noint_T0.3333333333333333_sim_simdata.txt"
# sim_data:SimulationData = loadsim("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/basic_N10_t1000.0_align_T1_sim_simdata.txt")
# sim_data:SimulationData = loadsim("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/basic_N10_t1000.0_align_T1_sim_simdata.txt")
sim_data:SimulationData = loadsim(file_str)


# # Create a figure and axes
# fig, ax = plt.subplots()
# xdata, ydata = [], []
# # ydata = np.zeros([1, sim_data._sim_params._num_particles])
# ln, = ax.plot([], [])

# # Create a figure and axes
# fig, ax = plt.subplots()
# xdata, ydata = [], []
# ln, = plt.plot([], [], 'b-', animated=True)

# def init():
#     ax.set_xlim(-sim_data._sim_params._box_width, sim_data._sim_params._box_width)
#     ax.set_ylim(-1.1, 1.1)
#     return ln

# def update(frame):
#     xdata.append(frame)
#     ydata.append(np.sin(frame))
#     ln.set_data(xdata, ydata)
#     return ln,

# # Generate x values over time (step size controls speed)
# ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128),
#                     init_func=init, blit=True, interval=50)

# plt.show()


# Create a figure and axes
fig, ax = plt.subplots()
xdata, ydata = [], []
# # ydata = np.zeros([1, sim_data._sim_params._num_particles])
ln, = plt.plot([], [], markersize=30, c="r", marker='.', linewidth=0, alpha=0.3, animated=True)

def init():
    ax.set_xlim(-sim_data._sim_params._box_width/2, sim_data._sim_params._box_width/2)
    ax.set_ylim(-.2, .2)
    return ln,

def update(frame):
    # print(frame)
    xdata = sim_data.wrapped_positions[frame]
    # print(sim_data.positions[frame])
    ydata = np.zeros([1, sim_data._sim_params._num_particles])
    # append(frame)
    # ydata.append(np.sin(frame))
    # print(np.concatenate((np.reshape(xdata, (-1,1)), np.reshape(ydata, (-1,1))), 1))
    # sc.set_offsets(np.concatenate((np.reshape(xdata, (-1,1)), np.reshape(ydata, (-1,1))), 1))
    ln.set_data(xdata, ydata)
    # ln.set_color(['red' if x == 1 else 'blue' for x in sim_data.spins[frame]])
    # sc.set_cmap('bwr')
    # sc.set_clim(-1, 1)
    # sc.set_array(sim_data.spins[frame])
    return ln,

# Generate x values over time (step size controls speed)
ani = FuncAnimation(fig, update, frames=np.arange(0, sim_data._sim_params.get_ntimes(), dtype=int),
                    init_func=init, blit=True, interval=sim_data._sim_params._dt*1000)

plt.show()