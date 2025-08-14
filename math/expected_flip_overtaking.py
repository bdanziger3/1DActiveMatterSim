import random as rnd
import numpy as np
import matplotlib.pyplot as plt

# experiment to see how many overtakings it takes for two colliding clusters of each size N to become one

n_max = 30
TRIALS = 1000

results = np.zeros(n_max)
for n in range (1, n_max + 1):

    trials_res = np.zeros(TRIALS)
    for i in range(TRIALS):
        total = 0
        rounds = 0

        while abs(total) < n:
            rounds += 1
            total += (2 * rnd.randint(0, 1) - 1)

        trials_res[i] = rounds
    
    print(max(trials_res))
    results[n-1] = np.mean(trials_res)


print(results)

plt.plot(np.arange(1, n_max + 1), results)
plt.xlabel("Number of particles in each cluster")
plt.ylabel("Expected number of trials")

plt.show()