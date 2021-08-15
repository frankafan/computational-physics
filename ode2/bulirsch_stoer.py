import numpy as np
from scipy.constants import G, au
import matplotlib.pyplot as plt
from time import time

M = 1.9891e30  # [kg] Mass of the Sun
P_E = 1.4710e11  # [m] perihelion of the Earth
V_E = 3.0287e4  # [m/s] velocity of Earth at perihelion
H = 604800.  # [s] time step

N_REVS = 1.1  # []  # of revolutions around the Sun
YEAR = 365.25 * 24 * 3600.  # [s] duration of a year
T = N_REVS * YEAR  # [s] duration of integration
N_STEPS = int(T / H) + 1  # [] number of time steps

ACCURACY = 1000 / YEAR  # Required position accuracy per unit time

# initialization of positions and  velocities
y_pos = []
x_pos = []
x_vel = []
y_vel = []
t_lst = []

# P_E = P_E / au
# V_E = V_E / au

def f(r_):
    x, y, vx, vy = r_[0], r_[1], r_[2], r_[3]
    r_ = (x ** 2 + y ** 2) ** .5  # [m] distance to Sun
    a_x = -G * M * x / r_ ** 3
    a_y = -G * M * y / r_ ** 3
    return np.array([vx, vy, a_x, a_y], float)


r = [P_E, 0, 0, V_E]

start = time()
for t in range(N_STEPS):
    t_lst.append(t)
    x_pos.append(r[0])
    y_pos.append(r[1])
    x_vel.append(r[2])
    y_vel.append(r[3])

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5 * H * f(r)
    r2 = r + H * f(r1)

    # The array R1 stores the first row of the
    # extrapolation table, which contains only the single
    # modified midpoint estimate of the solution at the
    # end of the interval
    R1 = np.empty([1, 4], float)
    R1[0] = 0.5 * (r1 + r2 + 0.5 * H * f(r2))

    # Now increase n until the required accuracy is reached
    error = 2 * H * ACCURACY
    while error > H * ACCURACY:

        n += 1
        h = H / n

        # Modified midpoint method
        r1 = r + 0.5 * h * f(r)
        r2 = r + h * f(r1)
        for i in range(n - 1):
            r1 += h * f(r2)
            r2 += h * f(r1)

        # Calculate extrapolation estimates.  Arrays R1 and R2
        # hold the two most recent lines of the table
        R2 = R1
        R1 = np.empty([n, 4], float)
        R1[0] = 0.5 * (r1 + r2 + 0.5 * h * f(r2))
        for m in range(1, n):
            epsilon = (R1[m - 1] - R2[m - 1]) / ((n / (n - 1)) ** (2 * m) - 1)
            R1[m] = R1[m - 1] + epsilon
        error = abs(epsilon[0])

    # Set r equal to the most accurate estimate we have,
    # before moving on to the next big step
    r = R1[n - 1]
end = time()
print("Bulirsch-Stoer runtime:", end - start)

for i in range(len(x_pos)):
    x_pos[i] = x_pos[i] / au
    y_pos[i] = y_pos[i] / au

# Plot the results
plt.plot(x_pos, y_pos)
# plt.plot(x_vel, y_vel)
plt.xlabel('x coordinate (AU)')
plt.ylabel('y coordinate (AU)')
plt.title('Orbit of the Earth using Bulirsch-Stoer method')
plt.savefig('Orbit of the Earth using Bulirsch-Stoer method')
plt.show()
