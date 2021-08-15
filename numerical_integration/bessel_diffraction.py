import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import jv

n_simpson = 1000


def f(theta, m, x):
    """Returns the function to be integrated."""
    return math.cos(m * theta - x * math.sin(theta))


def j_simpson(n, m, x):
    """Returns an approximation of the dawson function's output using the
    Simpson's method, given the number of slices (n), as well as the m and x
    values."""
    h = math.pi / n
    odd_sum = 0
    even_sum = 0
    for i in range(1, n):
        if i % 2 == 1:
            odd_sum += 4 * f(i * h, m, x)
        elif i % 2 == 0:
            even_sum += 2 * f(i * h, m, x)
    return h / 3 * (
                f(0, m, x) + f(math.pi, m, x) + odd_sum + even_sum) / math.pi


def j_simpson_plot(n):
    """Return a list of two dictionaries used for plotting, where the x_dict
    dictionary has keys corresponding to the m values of the Bessel function
    and values representing the x values, and the second dictionary has the same
    keys as the x_dict but with values giving the output of the Bessel function
    with their respective x values approximated using Simpson's method."""
    x_dict = {}
    y_dict = {}
    for i in range(3):
        x_lst = []
        y_lst = []
        for j in np.arange(0, 20.1, 0.1):
            x_lst.append(j)
            y_lst.append(j_simpson(n, i, j))
        x_dict[i] = x_lst
        y_dict[i] = y_lst
    return [x_dict, y_dict]


def j_scipy_plot():
    """Return a list of two dictionaries used for plotting, where the x_dict
    dictionary has keys corresponding to the m values of the Bessel function
    and values representing the x values, and the second dictionary has the same
    keys as the x_dict but with values giving the output of the Bessel function
    with their respective x values calculated using SciPy's jv function."""
    x_dict = {}
    y_dict = {}
    for i in range(3):
        x_lst = []
        y_lst = []
        for j in np.arange(0, 20.1, 0.1):
            x_lst.append(j)
            y_lst.append(jv(i, j))
        x_dict[i] = x_lst
        y_dict[i] = y_lst
    return [x_dict, y_dict]


def i(r):
    """Returns the function representing the intensity of the light in the
    diffraction pattern, given the distance in the focal plane from the center
    of the diffraction pattern (r)."""
    wavelength = 500
    k = 2 * math.pi / wavelength
    return (jv(1, k * r) / (k * r)) ** 2


# def i_plot():
#     """Returns a list of two lists for plotting, where r_list is a list of r
#      representing the distance in the focal plane from the center of the
#      diffraction pattern, and i_list is the intensity given the radius."""
#     max_r = 1000
#     r_lst = []
#     i_lst = []
#     for r in range(1, max_r):
#         r_lst.append(r)
#         i_lst.append(i(r))
#     return [r_lst, i_lst]


def density_plot():
    """Return a matrix array for plotting, where each array in y_lst is a row,
    and each term the x_lst is the intensity at that column."""
    max_r = 1000
    y_lst = []
    # for y in range(1, max_r):
    #     x_lst = []
    #     for x in range(1, max_r):
    #         r = (x ** 2 + y ** 2) ** 0.5
    #         x_lst.append(i(r))
    #     y_lst.append(x_lst)
    #     print(y / (max_r / 100), "% done")
    for y in range(-max_r, max_r):
        x_lst = []
        for x in range(-max_r, max_r):
            r = (abs(x) ** 2 + abs(y) ** 2) ** 0.5
            if r != 0:
                x_lst.append(i(r))
            elif r == 0:
                x_lst.append(0.25)
        y_lst.append(x_lst)
    return y_lst


plt.figure()
plt.plot(j_simpson_plot(n_simpson)[0][0], j_simpson_plot(n_simpson)[1][0], label='$J_0$')
plt.plot(j_simpson_plot(n_simpson)[0][1], j_simpson_plot(n_simpson)[1][1], label='$J_1$')
plt.plot(j_simpson_plot(n_simpson)[0][2], j_simpson_plot(n_simpson)[1][2], label='$J_2$')
plt.title("Bessel function estimated with Simpson's rule")
plt.xlabel("x")
plt.ylabel("$J_m$(x)")
plt.legend()
plt.savefig("Bessel function estimated with Simpson's rule")

plt.figure()
plt.plot(j_scipy_plot()[0][0], j_scipy_plot()[1][0], label='$J_0$')
plt.plot(j_scipy_plot()[0][1], j_scipy_plot()[1][1], label='$J_1$')
plt.plot(j_scipy_plot()[0][2], j_scipy_plot()[1][2], label='$J_2$')
plt.title("Bessel function plotted by SciPy jv function")
plt.xlabel("x")
plt.ylabel("$J_m$(x)")
plt.legend()
plt.savefig("Bessel function plotted by SciPy jv function")

plt.figure()
plt.imshow(density_plot(), 'hot', vmax=0.005)
plt.title("Density plot of diffraction pattern intensity")
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
plt.colorbar().set_label('Intensity')
plt.savefig("Density plot of diffraction pattern intensity")

# plt.figure()
# i_plot_lst = i_plot()
# plt.plot(i_plot_lst[0], i_plot_lst[1])
plt.show()
