import numpy as np
from scipy import integrate
from time import time


def f(x):
    return np.exp(x ** 2)


def simpson(n, a, b, func):
    h = abs(b - a) / n
    x_lst = []
    y_lst = []
    start = time()
    for i in range(n + 1):
        x_lst.append(i * h)
        y_lst.append(func(i * h))
    end = time()
    print("runtime:", end - start)
    return integrate.simps(y_lst, x_lst)


def simpson_(n, a, b, func):
    h = abs(b - a) / n
    odd_sum = 0
    even_sum = 0
    start = time()
    for i in range(1, n):
        if i % 2 == 1:
            odd_sum += 4 * func(i * h)
        elif i % 2 == 0:
            even_sum += 2 * func(i * h)
    end = time()
    print("runtime:", end - start)
    return h / 3 * (func(a) + func(b) + odd_sum + even_sum)


print(simpson(1000, 0, 5, f))
print(simpson_(1000, 0, 5, f))
