from math import exp, pi
import matplotlib.pyplot as plt
import scipy.constants as spc

ROOT = 0


def root_f(x):
    """Function equal to 0"""
    return 5 * exp(-x) + x - 5


def root_df_dx(x):
    """Derivative of the function equal to 0"""
    return 1 - 5 * exp(-x)


def f(x):
    """Function equal to x"""
    return 5 - 5 * exp(-x)


def df_dx(x):
    """Derivative of the function equal to x"""
    return 5 * exp(-x)


def plot_equation(root, step_size, points):
    x_lst = []
    y_lst = []
    root_lst = []
    for i in range(points):
        x_lst.append(i * step_size)
        y_lst.append(root_f(i * step_size))
        root_lst.append(root)

    plt.plot(x_lst, y_lst)
    plt.plot(x_lst, root_lst)
    plt.title("Nonlinear equation")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.savefig("Q3(c) nonlinear equation")
    plt.show()


def binary_search_recursive(x_low, x_high, root):
    if x_low <= x_high:
        x_mid = (x_low + x_high) // 2
        if root_f(x_mid) == root:
            return x_mid
        elif root_f(x_mid) > root:
            return binary_search_recursive(x_low, x_mid, root)
        elif root_f(x_mid) < root:
            return binary_search_recursive(x_mid, x_high, root)
    else:
        return -1


def binary_search_iterative(x_low, x_high, root):
    i = 0
    while x_low <= x_high:
        i += 1
        x_mid = (x_low + x_high) / 2
        if root_f(x_mid) < root:
            x_low = x_mid
        elif root_f(x_mid) > root:
            x_high = x_mid
        else:
            print("x =", x_mid)
            print("Iterations:", i)
            return x_mid
    return -1


def relaxation(x_start, accuracy):
    x1 = x_start
    error = 1.0
    i = 0
    while error > accuracy:
        i += 1
        x1, x2 = f(x1), x1
        error = abs((x1 - x2) / (df_dx(x2)))
    print("x =", x1)
    print("Iterations:", i)
    return x1


def newton_raphson(x_start):
    x = x_start
    i = 0
    while root_f(x) != 0:
        i += 1
        x = x - root_f(x) / root_df_dx(x)
    print("x =", x)
    print("Iterations:", i)
    return x


def b_constant(x):
    return (spc.h * spc.c) / (spc.k * x)


def intensity(wavelength, x):
    wavelength = wavelength * 1e-9
    temperature = b_constant(x) / wavelength
    numerator = 2 * pi * spc.h * (spc.c ** 2) * (wavelength ** -5)
    denominator = exp((spc.h * spc.c) / (wavelength * spc.k * temperature)) - 1
    print(numerator / denominator)
    return numerator / denominator


plot_equation(ROOT, step_size=0.01, points=1000)
print("Binary search method:")
binary_search_iterative(0.01, 100, ROOT)
print()
print("Relaxation method:")
relaxation(x_start=100.0, accuracy=1e-6)
print()
print("Newton's method:")
x = newton_raphson(100)
print()
print("Surface temperature of the Sun (K):", b_constant(x) / 5.02e-7)
