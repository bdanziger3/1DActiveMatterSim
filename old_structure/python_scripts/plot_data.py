import numpy as np
import matplotlib.pyplot as plt

f = open("data1.txt")

data_s = f.readline()

f.close()

data_array = [float(x) for x in data_s.replace(" ", "").split(",")]

print(data_array)

plt.plot(data_array, range(len(data_array)))
plt.show()