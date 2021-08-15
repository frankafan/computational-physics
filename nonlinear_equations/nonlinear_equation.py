from math import exp
from numpy import linspace
import matplotlib.pyplot as plt


def f(x, c):
    """Function equal to x"""
    return 1 - exp(-c * x)


def df_dx(x, c):
    """Derivative of the function equal to x"""
    return c * exp(-c * x)


# relaxation method
def relaxation(c, x_start, accuracy):
    x1 = x_start
    error = 1.0
    i = 0
    while error > accuracy:
        i += 1
        x1, x2 = f(x1, c), x1
        error = abs((x1 - x2) / (df_dx(x2, c)))
    print("x =", x1)
    print("Error:", error)
    print("Iterations:", i)


# plotting varying c values
def c_plot(c_max, accuracy, step_size, points):
    y_lst = []
    c_lst = linspace(step_size, c_max, points)

    for c in c_lst:
        x1 = 1.0
        error = 1.0

        while error > accuracy:
            x1, x2 = f(x1, c), x1
            error = abs((x1 - x2) / (df_dx(x2, c)))
        y_lst.append(x1)

    plt.plot(c_lst, y_lst)
    plt.title("Illustration of percolation transition with varying c values")
    plt.xlabel("c")
    plt.ylabel("x")
    plt.savefig("Illustration of percolation transition with varying c values")
    plt.show()


# over-relaxation method
def over_relaxation(omega, c, x_start, accuracy):
    x1 = x_start
    error = 1.0
    i = 0
    while error > accuracy:
        i += 1
        x1, x2 = (1 + omega) * f(x1, c) - omega * x1, x1
        error = abs((x1 - x2) / (df_dx(x2, c)))
    print("x =", x1)
    print("Error:", error)
    print("Iterations:", i)


relaxation(c=2.0, x_start=1.0, accuracy=1e-6)
over_relaxation(omega=0.5, c=2.0, x_start=1.0, accuracy=1e-6)
c_plot(c_max=3.0, step_size=0.01, points=1000, accuracy=1e-6,)
