import numpy as np
import random
import matplotlib.pyplot as plt

SIGMA = 1
Z = 0.5
TMAX = 1
TMIN = 1e-2
TAU = 1e5

random.seed(1)


def f(x, y):
    return x ** 2 - np.cos(4 * np.pi * x) + (y - 1) ** 2


def r(z):
    return (-2 * SIGMA ** 2 * np.log(1 - z)) ** 0.5


def cartesian(r_, theta):
    return [r_ * np.cos(theta), r_ * np.sin(theta)]


def random_pair(z):
    return cartesian(r(z), random.random() * 2 * np.pi)


t = 0
T = TMAX
current_coord = [2, 2]
current_minimum = f(current_coord[0], current_coord[1])
x_lst = [current_coord[0]]
y_lst = [current_coord[1]]
f_lst = [current_minimum]
t_lst = [t]
while T > TMIN:
    t += 1
    T = TMAX * np.exp(-t / TAU)

    delta = random_pair(Z)
    new_coord = [current_coord[0] + delta[0], current_coord[1] + delta[1]]
    new_minimum = f(new_coord[0], new_coord[1])

    if (new_minimum < current_minimum) or (
            random.random() < np.exp(-(new_minimum - current_minimum) / T)):
        current_coord = new_coord
        current_minimum = f(new_coord[0], new_coord[1])

    x_lst.append(current_coord[0])
    y_lst.append(current_coord[1])
    f_lst.append(current_minimum)
    t_lst.append(t)

    # print(current_coord)
    # print(current_minimum)

print(current_coord)
print(current_minimum)

plt.figure()
plt.plot(t_lst, x_lst)
plt.xlabel("t")
plt.ylabel("x")
plt.title("Candidate x coordinate over time")
plt.savefig("Candidate x coordinate over time")


plt.figure()
plt.plot(t_lst, y_lst)
plt.xlabel("t")
plt.ylabel("y")
plt.title("Candidate y coordinate over time")
plt.savefig("Candidate y coordinate over time")


plt.figure()
plt.plot(t_lst, f_lst)
plt.xlabel("t")
plt.ylabel("f(x,y)")
plt.title("Candidate global minimum over time")
plt.savefig("Candidate global minimum over time")

plt.show()
