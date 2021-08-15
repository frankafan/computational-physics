import numpy as np
from numpy import e, dot
from random import random, randrange, randint
import matplotlib.pyplot as plt


def energyfunction(J_, dipoles):
    """ function to calculate energy of a 2D system """
    rows = dipoles.shape[0]  # Define number of rows
    cols = dipoles.shape[1]  # Number of columns
    energy = 0  # initialize energy
    for i in range(rows - 1):  # loop over all rows except the last one
        energy += dot(dipoles[i, :],
                      dipoles[i + 1, :])  # dot product of two adjacent rows
    for j in range(cols - 1):  # loop over all columns except the last one
        energy += dot(dipoles[:, j],
                      dipoles[:, j + 1])  # dot product of two adjacent columns

    energy = -J_ * energy
    return energy


def acceptance(Enew_, Eold_, p_):
    """ Function for acceptance probability """
    result = False
    if (Enew_ - Eold_) <= 0:
        result = True
    if ((Enew_ - Eold_) > 0) and p_ > random():
        result = True

    return result  # result is True of False


def getrandM(r, c):
    """Function generating a r by c matrix with randomly assigned elements
    of either 1 or -1 """
    M = np.zeros([r, c], int)
    values = np.array([1, -1])
    for i in range(r):
        for j in range(c):
            M[i, j] = values[randint(0, 1)]
    return M


# define constants
kB = 1.0
T = 3.0
J = 1.0
rows_dipoles = 20
cols_dipoles = 20
N = 100000

# generate array of dipoles and initialize diagnostic quantities
dipoles = getrandM(rows_dipoles, cols_dipoles)
energy = []  # empty list; to add to it, use energy.append(value)
magnet = []  # empty list; to add to it, use magnet.append(value)

for i in range(N + 1):
    rpicked = randrange(rows_dipoles)  # choose a row victim
    cpicked = randrange(cols_dipoles)  # choose a col victim
    Eold = energyfunction(J, dipoles)  # compute energy of the current state
    dipoles[rpicked, cpicked] *= -1  # propose to flip the victim
    Enew = energyfunction(J, dipoles)  # compute Energy of proposed new state
    p = e ** (-(Enew - Eold) / kB * T)  # compute boltzmann distribution
    # calculate acceptance probability
    accepted = acceptance(Enew, Eold, p)
    if not accepted:
        dipoles[rpicked, cpicked] *= -1  # if rejected, flip it back

    # store energy and magnetization
    E = energyfunction(J, dipoles)
    energy.append(E)
    magnet.append(np.sum(dipoles))
    # print(dipoles)

    if i % 10000 == 0:
        plt.figure()
        plt.imshow(dipoles)
        cbar = plt.colorbar()
        plt.clim(-1,1)
        cbar.set_label("Spin")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"System of dipoles with T = {T} at time step {i}")
        plt.savefig(f"System of dipoles with T = {T} at time step {i}.png")

# plot energy, magnetization
# plt.figure()
# plt.plot(energy, label="energy")
# plt.plot(magnet, label="magnet")
# plt.legend()

# plt.figure()
# plt.plot(magnet)
# plt.xlabel("Time")
# plt.ylabel("M")
# plt.title(f"Total magnetization of the system over {N} time steps")
# plt.savefig(f"Total magnetization of the system over {N} time steps")

plt.show()
