import random as rd
import numpy as np
import matplotlib.pyplot as plt


# test how chacarcteristic time relates to fliprate

def rand_lin_flip(dt, flip_rate):
    random_val = rd.random()
    return random_val < dt * flip_rate

flip_rate = .5
dt = .001
total_time = 10000

t = dt
i = 0
flips = []
while t <= total_time:
    if rand_lin_flip(dt, flip_rate):
        flips.append(i)

    t += dt
    i += 1

flip_times = np.diff(flips) * dt
fliprates = 1 / flip_times

print(len(flips))
print(flip_times)
print(fliprates)
print(np.mean(flip_times))
print(np.mean(fliprates))


plt.hist(flip_times, bins = 20, alpha=.5)
# plt.hist(fliprates, bins = 20, alpha=.5)
plt.show()




