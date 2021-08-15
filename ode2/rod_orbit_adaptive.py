import numpy as np
import matplotlib.pyplot as plt
from time import time

G = 1.0
M = 10.0
L = 2.0
T_START = 0.0
T_END = 10.0
H = 0.01
H_NA = 0.001
X_START = 1.0
Y_START = 0.0
V_X0 = 0.0
V_Y0 = 1.0
ACCURACY = 10e-6

x_lst = [X_START]
y_lst = [Y_START]
v_lst = []
r = np.array([X_START, V_X0, Y_START, V_Y0])

r_na = r[:]
x_na = x_lst[:]
y_na = y_lst[:]

t_lst = []
h_lst = []


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
    r = r + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return r


def accurate_h(r_1, r_2, h, acc):
    if r_1[0] - r_2[0] == 0 or r_1[2] - r_2[2] == 0:
        return 2 * h
    else:
        # print('r_1:', r_1)
        # print('r_2:', r_2)
        e_x = 1 / 30 * (r_1[0] - r_2[0])
        e_y = 1 / 30 * (r_1[2] - r_2[2])
        # print('e_x:', e_x)
        # print('e_y:', e_y)
        p = (30 * h * acc) / (e_x ** 2 + e_y ** 2) ** 0.5
        # return min([h * p ** 0.25, 2 * h])
        return h * p ** 0.25


t = T_START
start_t = time()
while t < T_END:
    r_1 = next_r(H, next_r(H, r))
    r_2 = next_r(2 * H, r)
    h_accurate = accurate_h(r_1, r_2, H, ACCURACY)
    r = next_r(h_accurate, r)

    t_lst.append(t)
    h_lst.append(h_accurate)

    t += h_accurate
    x_lst.append(r[0])
    y_lst.append(r[2])
    v_lst.append((r[1] ** 2 + r[3] ** 2) ** 0.5)
end_t = time()
print('Adaptive step size runtime:', abs(end_t - start_t))

t = T_START
start_t = time()
while t < T_END:
    r_na = next_r(H_NA, r_na)

    t += H_NA
    x_na.append(r_na[0])
    y_na.append(r_na[2])
end_t = time()
print('Non-adaptive step size runtime:', abs(end_t - start_t))

plt.figure()
plt.plot(x_lst, y_lst, 'k.', label="Adaptive steps")
plt.plot(x_na, y_na, label="Non-adaptive orbit")
plt.legend()
plt.xlabel('x coordinate (m)')
plt.ylabel('y coordinate (m)')
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)
plt.title('Orbit of object around rod with adaptive time steps')
plt.savefig('Orbit of object around rod with adaptive time steps')
plt.show()

plt.figure()
plt.plot(t_lst, h_lst)
plt.xlabel('Time (s)')
plt.ylabel('Step size (s)')
plt.title('Size of adaptive time steps over time')
plt.savefig('Size of adaptive time steps over time')
plt.show()

fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Step size (s)', color=color)
ax1.plot(t_lst, h_lst, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:red'
ax2.set_ylabel('Velocity (m/s)', color=color)  # we already handled the x-label with ax1
ax2.plot(t_lst, v_lst, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Size of adaptive time steps and velocity over time')
plt.savefig('Size of adaptive time steps and velocity over time')
plt.show()
