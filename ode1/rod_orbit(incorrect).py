import math
import numpy as np
import matplotlib.pyplot as plt

G = 1.0
M = 10.0
L = 2.0
T_START = 0.0
T_END = 10.0
T_STEPS = 1000
X_START = 1.0
Y_START = 0.0
V_X0 = 0.0
V_Y0 = 1.0

x = [X_START]
y = [Y_START]
r_x = np.array([X_START, V_X0])
r_y = np.array([Y_START, V_Y0])
h = abs(T_START - T_END) / T_STEPS


def f(r):
    radius = math.sqrt(x[-1] ** 2 + y[-1] ** 2)
    x_or_y = r[1]
    v = - (G * M * r[0]) / (radius ** 2 * math.sqrt(radius ** 2 + L ** 2 / 4))
    return np.array([x_or_y, v])


for t in range(T_STEPS):
    k1_x = h * f(r_x)
    k2_x = h * f(r_x + 0.5 * k1_x)
    k3_x = h * f(r_x + 0.5 * k2_x)
    k4_x = h * f(r_x + k3_x)
    r_x += (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6

    k1_y = h * f(r_y)
    k2_y = h * f(r_y + 0.5 * k1_y)
    k3_y = h * f(r_y + 0.5 * k2_y)
    k4_y = h * f(r_y + k3_y)
    r_y += (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6

    x.append(r_x[0])
    y.append(r_y[0])

# plt.plot(x, y)
# plt.xlabel('x coordinate (m)')
# plt.ylabel('y coordinate (m)')
# plt.xlim(-1.1, 1.1)
# plt.ylim(-1.1, 1.1)
# plt.title('Orbit of object around rod')
# plt.savefig('Orbit of object around rod')
# plt.show()
