import numpy as np
from scipy.special import dawsn
import matplotlib.pyplot as plt
from gaussxw import gaussxwab
from math import exp

x_at = 4


def f(t):
    """Returns the function to be integrated."""
    return exp(t ** 2)


def d_trapezoidal(n, x_at):
    """Returns an approximation of the dawson function's output using the
    trapezoidal method, given the number of slices (n) and x value."""
    h = x_at / n
    s = 0.5 * f(0) + 0.5 * f(x_at)
    for i in range(1, n):
        s += f(i * h)
    return h * s * exp(-x_at ** 2)


def d_simpson(n, x_at):
    """Returns an approximation of the dawson function's output using the
    Simpson's method, given the number of slices (n) and x value."""
    if n % 2 != 0:
        return None
    h = x_at / n
    odd_sum = 0
    even_sum = 0
    for i in range(1, n):
        if i % 2 == 1:
            odd_sum += 4 * f(i * h)
        elif i % 2 == 0:
            even_sum += 2 * f(i * h)
    return h / 3 * (f(0) + f(x_at) + odd_sum + even_sum) * exp(-x_at ** 2)


def d_gaussian(n, x_at):
    x, w = gaussxwab(n, 0, x_at)
    s = 0.0
    for k in range(n):
        s += w[k] * f(x[k])
    return s.item() * exp(-x_at ** 2)


def gaussian_error(n):
    """Returns the practical estimation of error by Gaussian approximation
    given the number of sample points."""
    I1 = d_gaussian(n, x_at)
    I2 = d_gaussian(n * 2, x_at)
    return abs(I1 - I2)


num_sample = []

trapezoidal_lst = []
simpson_lst = []
gaussian_lst = []

trapezoidal_errors = []
simpon_errors = []
gaussian_errors = []

trap_errors_prac = []
simp_errors_prac = []
gauss_errors_prac = []

for i in range(3, 11 + 1):
    num_sample.append(2 ** i)

    trapezoidal = d_trapezoidal(2 ** i, x_at)
    simpson = d_simpson(2 ** i, x_at)
    gaussian = d_gaussian(2 ** i, x_at)

    trapezoidal_lst.append(trapezoidal)
    simpson_lst.append(simpson)
    gaussian_lst.append(gaussian)

    trapezoidal_errors.append(abs(dawsn(x_at) - trapezoidal) / dawsn(x_at))
    simpon_errors.append(abs(dawsn(x_at) - simpson) / dawsn(x_at))
    gaussian_errors.append(abs(dawsn(x_at) - gaussian) / dawsn(x_at))
    gauss_errors_prac.append(gaussian_error(2 ** i))

num_sample2 = []

gaussian_errors2 = []
gauss_errors_prac2 = []

for i in range(8, 2048 + 1, 32):
    num_sample2.append(i)
    gaussian_errors2.append(
        abs(dawsn(x_at) - d_gaussian(i, x_at)) / dawsn(x_at))
    gauss_errors_prac2.append(gaussian_error(i))

plt.figure()
plt.plot(num_sample, trapezoidal_lst, label="Trapezoidal method")
plt.plot(num_sample, simpson_lst, label="Simpson's method")
plt.plot(num_sample, gaussian_lst, label="Gaussian approximation")
plt.title("Dawson function estimated with three different methods")
plt.xlabel("Number of slices / sample points")
plt.ylabel("D(x)")
plt.legend()
plt.savefig("Dawson function estimated with three different methods")

plt.figure()
plt.plot(num_sample, trapezoidal_lst, label="Trapezoidal method")
plt.plot(num_sample, simpson_lst, label="Simpson's method")
plt.plot(num_sample, gaussian_lst, label="Gaussian approximation")
plt.title("Dawson function estimated with three different methods (zoomed in)")
plt.xlabel("Number of slices / sample points")
plt.ylabel("D(x)")
plt.legend()
plt.xlim(0, 150)
plt.savefig(
    "Dawson function estimated with three different methods (zoomed in)")

plt.figure()
plt.plot(num_sample, trapezoidal_errors, label="Trapezoidal method error")
plt.plot(num_sample, simpon_errors, label="Simpson's method error")
plt.plot(num_sample, gaussian_errors, label="Gaussian approximation error")
plt.plot(num_sample, gauss_errors_prac,
         label="Gaussian approximation error (practical estimate)")
plt.title(
    "Relative error from true value of Dawson function (exponential steps)")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of slices / sample points")
plt.ylabel("Relative error")
plt.legend()
plt.savefig(
    "Relative error from true value of Dawson function (exponential steps)")

plt.figure()
plt.plot(num_sample2, gaussian_errors2, label="Gaussian approximation error")
plt.plot(num_sample2, gauss_errors_prac2,
         label="Gaussian approximation error (practical estimate)")
plt.title("Relative error from true value of Dawson function (linear steps)")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of slices / sample points")
plt.ylabel("Relative error")
plt.legend()
plt.savefig("Relative error from true value of Dawson function (linear steps)")

plt.show()
