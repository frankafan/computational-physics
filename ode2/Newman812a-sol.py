# Author: Nicolas Grisouard, University of Toronto
# Solution to Newman 8.12, ""Orbit of the Earth", part (a), computing orbit

from scipy.constants import G, au  # I should be using the same number of
# SigFigs quoted in the text but YOLO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from time import time


def f(r_):
    """ Right-hand-side for velocity equation: gravitational force
    INPUT: r_ = [x, y, vx, vy]
    OUTPUT: Fx, Fy: along-x and -y graviational forces
    """
    x, y, vx, vy = r_[0], r_[1], r_[2], r_[3]
    from scipy.constants import G  # computationally inefficient but short
    M = 1.9891e30  # [kg] Mass of the Sun
    r = (x ** 2 + y ** 2) ** .5  # [m] distance to Sun
    prefac = -M * G / r ** 3
    return np.array([vx, vy, x * prefac, y * prefac], float)


ftsz = 16
font = {'family': 'normal', 'size': ftsz}  # font size
rc('font', **font)

M = 1.9891e30  # [kg] Mass of the Sun
PH = 1.4710e11  # [m] perihelion of the Earth
VP = 3.0287e4  # [m/s] velocity of Earth at perihelion
h = 3600.  # [s] time step

Nrevs = 5  # []  # of revolutions around the Sun
year = 365.25 * 24 * 3600.  # [s] duration of a year
T = Nrevs * year  # [s] duration of integration
Nsteps = int(T / h)  # [] number of time steps

# initialization of positions and  velocities
pos = np.empty((2, Nsteps), float)  # 1st index: x, y
vel = np.empty((2, Nsteps), float)  # idem for v

# Jump-start with RK2 --------------------------------------------------------|
# I define the x, y axes as: Earth start at perielion, along x>0, with along-y
# positive velocity
v0 = np.array([0., VP])  # initial velocity components
pos[:, 0] = np.array([PH, 0.])
r = (pos[0, 0] ** 2 + pos[1, 0] ** 2) ** .5
vel[:, 0] = v0 - h * G * M * pos[:, 0] / r ** 3

start = time()
# Verlet iterations ----------------------------------------------------------|
for tt in range(1, Nsteps):
    pos[:, tt] = pos[:, tt - 1] + h * vel[:, tt - 1]  # updating positions
    r = (pos[0, tt] ** 2 + pos[1, tt] ** 2) ** .5  # distance from Sun squared
    vel[:, tt] = vel[:, tt - 1] - h * G * M * pos[:,
                                              tt] / r ** 3  # updating velocities
end = time()
print('Verlet runtime:', end - start)

# Plot -----------------------------------------------------------------------|
pos_AU = pos / au  # for clearer axes
plt.figure(dpi=100)
plt.plot(pos_AU[0, :], pos_AU[1, :])
plt.axvline(0.)
plt.axhline(0.)
plt.grid()
plt.xlabel('$x$ (AU)')
plt.ylabel('$y$ (AU)')
plt.axis('equal')
plt.tight_layout()
plt.show()
