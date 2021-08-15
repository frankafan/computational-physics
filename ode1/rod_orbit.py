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
ACCURACY = 10e-6

x_lst = [X_START]
y_lst = [Y_START]
r = np.array([X_START, V_X0, Y_START, V_Y0])
h = abs(T_START - T_END) / T_STEPS


def f(r):
    x = r[0]
    y = r[2]
    radius = np.sqrt(x ** 2 + y ** 2)
    v_x = - (G * M * x) / (radius ** 2 * np.sqrt(radius ** 2 + L ** 2 / 4))
    v_y = - (G * M * y) / (radius ** 2 * np.sqrt(radius ** 2 + L ** 2 / 4))
    return np.array([r[1], v_x, r[3], v_y], float)


def next_r(h, r):
    k1 = h * f(r)
    k2 = h * f(r + 0.5 * k1)
    k3 = h * f(r + 0.5 * k2)
    k4 = h * f(r + k3)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return r

for t in range(T_STEPS):
    r = next_r(h, r)

    x_lst.append(r[0])
    y_lst.append(r[2])

plt.plot(x_lst, y_lst)
plt.xlabel('x coordinate (m)')
plt.ylabel('y coordinate (m)')
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)
plt.title('Orbit of object around rod')
# plt.savefig('Orbit of object around rod')
plt.show()
