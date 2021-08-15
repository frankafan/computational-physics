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
N = 100  # number of particles

anchored = np.zeros((Lp, Lp), dtype=int)  # array to represent whether each grid-point has an anchored particle
anchored_points = [[], []]  # list to represent x and y positions of anchored points
centre_point = (Lp - 1) // 2  # middle point of domain


def notstuck(x, y, anchored_, Lp_=101):
    """Return True if the particle of the given coordinates is not at the edge
    nor beside any anchored particles"""
    if 0 < x < (Lp_ - 1) and 0 < y < (Lp_ - 1):  # if the particle is not at the edge
        if anchored_[yp + 1][xp] != 1 and anchored_[yp - 1][xp] != 1:  # if there are no anchored particles above or below the given particle
            if anchored_[yp][xp + 1] != 1 and anchored_[yp][xp - 1] != 1:  # if there are no anchored particles left or right of the given particle
                return True
    return False


x_last = [centre_point]
y_last = [centre_point]
for i in range(N):
    xp = centre_point
    yp = centre_point
    while notstuck(xp, yp, anchored, Lp):
        xp, yp = nextmove(xp, yp)
        if i == N - 1:  # if this is the last particle to iterate
            # store trajectory of the last particle
            x_last.append(xp)
            y_last.append(yp)
    anchored_points[0].append(xp)
    anchored_points[1].append(yp)
    anchored[yp][xp] = 1

plt.plot(anchored_points[0], anchored_points[1], '.r', label="Stuck particle")
plt.plot(x_last, y_last, label="Last particle trajectory")
plt.legend()
plt.xlim([-1, Lp])
plt.ylim([-1, Lp])
plt.title(f"Trajectory of the {N}th particle")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig(f"Trajectory of the {N}th particle")
plt.show()
