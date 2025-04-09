import numpy as np
import random as rd
import matplotlib.pyplot as plt

def expected_n(n, runs=100):
    s_i = 0
    for i in range(n):
        s_i += int((2 * rd.randint(0,1)) - 1)




# def random_walk_1d(steps, dx=1, print_steps=False):

#     x = 0
#     for _ in range(steps):
#         x += int((2 * rd.randint(0,1)) - 1)
        
#         if print_steps:
#             print(x)
    
#     return x