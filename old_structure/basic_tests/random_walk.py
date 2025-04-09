import numpy as np
import random as rd
import matplotlib.pyplot as plt

def random_walk_1d(steps, dx=1, print_steps=False):

    x = 0
    for _ in range(steps):
        x += int((2 * rd.randint(0,1)) - 1)
        
        if print_steps:
            print(x)
    
    return x

def random_walk_2d(steps, dx=1, dy=1, theta_window=2*np.pi, print_steps=False):

    x = 0
    y = 0
    for _ in range(steps):
        theta = theta_window * (rd.random() - 0.5)

        x += dx * np.cos(theta)
        y += dy * np.sin(theta)
        
        if print_steps:
            print(x, y)
    
    return (x, y)


def random_walk_sim_1d(max_steps, runs_per_n):
    mean_x = np.zeros(max_steps - 1)
    mean_x_sq = np.zeros(max_steps - 1)
    ns  = []

    for n in range(1, max_steps):
        ns.append(n)
        res = []
        
        for _ in range(runs_per_n):
            res.append(random_walk_1d(n))

        res = np.array(res)
        res_sq = np.array(res) ** 2
        res_abs = np.abs(np.array(res))

        mean_x[n-1] = np.mean(res)
        mean_x_sq[n-1] = np.mean(res_abs)


    return mean_x, mean_x_sq


def random_walk_sim_2d(max_steps, runs_per_n, theta_window=2*np.pi):
    mean_x = np.zeros(max_steps - 1)
    mean_x_sq = np.zeros(max_steps - 1)

    mean_y = np.zeros(max_steps - 1)
    mean_y_sq = np.zeros(max_steps - 1)

    mean_r = np.zeros(max_steps - 1)
    mean_r_sq = np.zeros(max_steps - 1)

    ns  = []

    for n in range(1, max_steps):
        ns.append(n)
        res_x = []
        res_y = []
        res_r = []
        
        for _ in range(runs_per_n):
            x, y = random_walk_2d(n, theta_window=theta_window)

            res_x.append(x)
            res_y.append(y)
            res_r.append(np.sqrt(x**2 + y**2))

        res_x = np.array(res_x)
        res_x_sq = np.array(res_x) ** 2
        res_y = np.array(res_y)
        res_y_sq = np.array(res_y) ** 2
        res_r = np.array(res_r)
        res_r_sq = np.array(res_r) ** 2

        mean_x[n-1] = np.mean(res_x)
        mean_x_sq[n-1] = np.mean(res_x_sq)

        mean_y[n-1] = np.mean(res_y)
        mean_y_sq[n-1] = np.mean(res_y_sq)

        mean_r[n-1] = np.mean(res_r)
        mean_r_sq[n-1] = np.mean(res_r_sq)


    return mean_x, mean_x_sq, mean_y, mean_y_sq, mean_r, mean_r_sq


## 2D SIM
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_0 = random_walk_sim_2d(100, 100)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_1 = random_walk_sim_2d(100, 100, theta_window=np.pi)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_3 = random_walk_sim_2d(100, 100, theta_window=np.pi / 3)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_6 = random_walk_sim_2d(100, 100, theta_window=np.pi / 4)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_4 = random_walk_sim_2d(100, 100, theta_window=np.pi / 6)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_8 = random_walk_sim_2d(100, 100, theta_window=np.pi / 8)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_12 = random_walk_sim_2d(100, 100, theta_window=np.pi / 12)
# mean_x_2d, mean_x_sq_2d, mean_y_2d, mean_y_sq_2d, mean_r_2d, mean_r_sq_2d_12 = random_walk_sim_2d(100, 100, theta_window=0)
# ns = np.arange(1, len(mean_x_2d) + 1)
# # plt.plot(ns, mean_x_2d)
# # plt.plot(ns, mean_x_sq_2d)
# # plt.plot(ns, mean_y_2d)
# # plt.plot(ns, mean_y_sq_2d)
# # plt.plot(ns, mean_r_2d)
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_0))
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_1))
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_3))
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_4))
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_6))
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_8))
# plt.plot(np.log(ns), np.log(mean_r_sq_2d_12))
# # plt.plot(ns, ns)
# # plt.plot(ns, ns / 2)
# # plt.plot(ns, np.sqrt(ns))
# plt.legend(["2pi", "pi", "pi/3", "pi/4", "pi/6", "pi/8", "pi/16", "0"])
# plt.show()


### 1D SIM
mean_x_1d, mean_x_sq_1d = random_walk_sim_1d(100, 500)
ns = np.arange(1, len(mean_x_1d) + 1)
plt.plot(ns, mean_x_1d)
plt.plot(ns, mean_x_sq_1d / np.sqrt(ns))
# plt.plot(ns, (1 / np.sqrt(2)) * np.sqrt(ns))
plt.plot(ns, (4/5) * np.sqrt(ns))
# # plt.plot(ns, (ns) * np.sqrt(ns))
# plt.plot(ns,  (5/6) * np.sqrt(ns))
plt.plot(ns, mean_x_sq_1d)
plt.show()