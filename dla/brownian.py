import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import random


def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    if 0 < x < (Lp - 1) and 0 < y < (Lp - 1):  # if the particle is not at the edge
        direction = random.randint(0, 3)
    else:
        options = np.arange(4).tolist()  # array of the possible directions
        # remove the illegal move depending on the particle's edge case
        if x == 0:
            options.pop(3)
        elif y == 0:
            options.pop(1)
        elif x == (Lp - 1):
            options.pop(2)
        elif y == (Lp - 1):
            options.pop(0)
        direction = options[random.randint(0, 2)]

    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y


font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)

Lp = 101  # size of domain
Nt = 5000  # number of time steps

centre_point = (Lp - 1) // 2  # middle point of domain

# arrays to record the trajectory of the particle
x_timeline = [centre_point]
y_timeline = [centre_point]

for i in range(Nt):
    xp, yp = nextmove(x_timeline[-1], y_timeline[-1])
    x_timeline.append(xp)
    y_timeline.append(yp)

plt.plot(x_timeline, y_timeline)
plt.xlim([-1, Lp])
plt.ylim([-1, Lp])
plt.title(f"Trajectory of the random walk particle")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig(f"Trajectory of the random walk particle")
plt.show()
