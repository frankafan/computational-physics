from math import sqrt, exp
from numpy import empty
from random import random, randrange, seed
import matplotlib.pyplot as plt

N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e4


# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0] ** 2 + x[1] ** 2)


# Function to calculate the total length of the tour
def distance():
    s = 0.0
    for i in range(N):
        s += mag(r[i + 1] - r[i])
    return s


# Choose N city locations and calculate the initial distance
r = empty([N + 1, 2], float)
seed(1)
for i in range(N):
    r[i, 0] = random()
    r[i, 1] = random()
r[N] = r[0]
D = distance()

seed_lst = [2, 5, 8, 11, 14]
for s in seed_lst:
    seed(s)

    # Main loop
    t = 0
    T = Tmax
    while T > Tmin:

        # Cooling
        t += 1
        T = Tmax * exp(-t / tau)

        # Update the visualization every 100 moves
        if t % 100 == 0:
            r_new = r

        # Choose two cities to swap and make sure they are distinct
        i, j = randrange(1, N), randrange(1, N)
        while i == j:
            i, j = randrange(1, N), randrange(1, N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i, 0], r[j, 0] = r[j, 0], r[i, 0]
        r[i, 1], r[j, 1] = r[j, 1], r[i, 1]
        D = distance()
        deltaD = D - oldD

        # If the move is rejected, swap them back again
        if random() > exp(-deltaD / T):
            r[i, 0], r[j, 0] = r[j, 0], r[i, 0]
            r[i, 1], r[j, 1] = r[j, 1], r[i, 1]
            D = oldD

    plt.figure()
    plt.plot(*zip(*r), '.', label='Cities')
    plt.plot(*zip(*r_new), 'r', label='Path')
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Path with scheduling time {tau} and random seed {s} has total distance {D}")
    plt.savefig(f"Path with scheduling time {tau} and random seed {s} has total distance {D}.png")

    # plt.figure()
    # plt.plot(*zip(*r), '.', label='Cities')
    # plt.legend()
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.title("Locations of randomly generated cities")

plt.show()
