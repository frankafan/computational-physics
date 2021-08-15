import numpy as np
from scipy.special import dawsn
from time import time

n_trapezoidal = 8
n_simpson = 8
n_trapezoidal_o9 = 32660
n_simpson_o9 = 532
x_at = 4


def f(t):
    """Returns the function to be integrated."""
    return np.exp(t ** 2)


def d_trapezoidal(n, x):
    """Returns an approximation of the dawson function's output using the
    trapezoidal method, given the number of slices (n) and x value."""
    h = x / n
    s = 0.5 * f(0) + 0.5 * f(x)
    for i in range(1, n):
        s += f(i * h)
    return h * s * np.exp(-x ** 2)


def d_simpson(n, x):
    """Returns an approximation of the dawson function's output using the
    Simpson's method, given the number of slices (n) and x value."""
    if n % 2 != 0:
        return None
    h = x / n
    odd_sum = 0
    even_sum = 0
    for i in range(1, n):
        if i % 2 == 1:
            odd_sum += 4 * f(i * h)
        elif i % 2 == 0:
            even_sum += 2 * f(i * h)
    return h / 3 * (f(0) + f(x) + odd_sum + even_sum) * np.exp(-x ** 2)


def trapezoidal_error():
    """Returns the practical estimation of error by the trapezoidal method."""
    I1 = d_trapezoidal(32, 4)
    I2 = d_trapezoidal(64, 4)
    return abs(I1 - I2) / 3


def simpson_error():
    """Returns the practical estimation of error by the Simpson's method."""
    I1 = d_simpson(32, 4)
    I2 = d_simpson(64, 4)
    return abs(I1 - I2) / 15


print("SciPy function:", dawsn(x_at))
print("Trapezoidal method:", d_trapezoidal(n_trapezoidal, x_at))
print("Simpson's method:", d_simpson(n_simpson, x_at))

print()
print("Trapezoidal method error with N = 8:", abs(dawsn(x_at) - d_trapezoidal(n_trapezoidal, x_at)))
print("Simpson's method error with N = 8:", abs(dawsn(x_at) - d_simpson(n_simpson, x_at)))

print()
print("Trapezoidal method error with practical estimation:", trapezoidal_error())
print("Simpson's method error with practical estimation:", simpson_error())

print()
print("Number of slices (N) for the trapezoidal estimation to reach O(-9) error:", n_trapezoidal_o9)
start = time()
print("Trapezoidal method estimation with O(-9) error:", d_trapezoidal(n_trapezoidal_o9, x_at))
end = time()
print("runtime(s):", end - start)
print("Trapezoidal method error with N =", n_trapezoidal_o9, ":", abs(dawsn(x_at) - d_trapezoidal(n_trapezoidal_o9, x_at)))

print()
print("Number of slices (N) for the Simpson's estimation to reach O(-9) error:", n_simpson_o9)
start = time()
print("Simpson's method estimation with O(-9) error:", d_simpson(n_simpson_o9, x_at))
end = time()
print("runtime(s):", end - start)
print("Simpson's method error with N =", n_simpson_o9, ":", abs(dawsn(x_at) - d_simpson(n_simpson_o9, x_at)))
