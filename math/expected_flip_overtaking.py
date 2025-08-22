import random as rnd
import numpy as np
import matplotlib.pyplot as plt

# experiment to see how many overtakings it takes for two colliding clusters of each size N to become one
def num_trials():
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


# experiment to see the probability of winning by starting position

def win_pct(n_max):
    TRIALS = 3000

    results = np.zeros(n_max-1)
    wins_each = np.zeros(n_max-1)
    # for n in range (1, n_max):
    for n in range (1, int(np.ceil(n_max/2))+1):

        trials_rounds = np.zeros(TRIALS)
        wins = 0
        for i in range(TRIALS):
            # start at n
            position = n
            rounds = 0

            while 0 < position < n_max:
                rounds += 1
                position += (2 * rnd.randint(0, 1) - 1)

            trials_rounds[i] = rounds
            wins += (position == n_max)
        
        results[n-1] = np.mean(trials_rounds)
        wins_each[n-1] = wins

    wins_each = wins_each / TRIALS
    print(results)
    print(wins_each)

    # plt.plot(np.arange(1, n_max), wins_each)
    plt.plot(np.arange(1, n_max) / (n_max - np.arange(1, n_max)), wins_each, ".")
    plt.plot(np.arange(1, n_max) / (n_max - np.arange(1, n_max)), 1-wins_each, ".")
    plt.xlabel("Number of particles in my cluster to start")
    plt.ylabel("Win Percentage")
    plt.tick_params("both", direction="in")
    plt.xlim([0, 1])
    plt.show()


    plt.plot(np.arange(1, n_max), results)
    plt.xlabel("Number of particles in my cluster to start")
    plt.ylabel("Average Number of rounds before someone wins")
    plt.tick_params("both", direction="in")
    plt.show()

win_pct(50)