#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:23:32 2020

@author: Chen Zhang
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import empty, zeros, e, mean
from pylab import plot, xlabel, ylabel, show

# Constants
L = 1  # The horizontal distance of water surface in meters
J = 50  # Number of divisions in grid
dx = L / J  # Grid spacing
dt = 0.01  # Time-step in second
g = 9.81  # gravitational constant
eta_b = 0  # bottom Topography
H = 0.01  # Water column height at rest
epsilon = dt / 1000


# Functions
def getF1(u_, eta_):
    """
    The first component of force function for surface wave.
    """
    return (1 / 2) * (u_ ** 2) + g * eta_


def getF2(u_, eta_):
    """
    The second component of force function for surface wave.
    """
    return (eta_ - eta_b) * u_


# Constants for initial n condition
x = np.linspace(0, L, J + 1)  # coordinate of x in meters
A = 0.002  # in meters
mu = 0.5  # in meters
sigma = 0.05  # in meters

# Time(seconds) to plot graphs
t1 = 0.01
t2 = 1.0
t3 = 5.0
# t4 = 10.0
tend = t3 + epsilon

# main program

# Create arrays for u and eta with initial conditions
eta = H + A * e ** (-(x - mu) ** 2 / sigma ** 2) - mean(
    A * e ** (-(x - mu) ** 2 / sigma ** 2))
u = zeros(J + 1, float)  # u vector consisting of u and eta.
# Main loop
t = 0.0
while t < tend:
    # Calculate the new values of T
    # initial arrays for n+1
    eta_new = zeros(J + 1, float)
    u_new = zeros(J + 1, float)
    # updataing eta at end points separately
    eta_new[0] = eta[0] - (dt / dx) * (
                getF2(u[1], eta[1]) - getF2(u[0], eta[0]))
    eta_new[J] = eta[J] - (dt / dx) * (
                getF2(u[J], eta[J]) - getF2(u[J - 1], eta[J - 1]))
    # updating all values inbetween, not that u at 0, and J remain 0

    for j in np.arange(1, J):
        u_new[j] = u[j] - (dt / (2 * dx)) * (
                    getF1(u[j + 1], eta[j + 1]) - getF1(u[j - 1], eta[j - 1]))
        eta_new[j] = eta[j] - (dt / (2 * dx)) * (
                    getF2(u[j + 1], eta[j + 1]) - getF2(u[j - 1], eta[j - 1]))

    eta = np.copy(eta_new)
    u = np.copy(u_new)
    t += dt

    # Make plots at the given times
    if abs(t - t1) < epsilon:
        plot(x, eta, label="t=0.1s")
    if abs(t - t2) < epsilon:
        plot(x, eta, label="t=1s")
    if abs(t - t3) < epsilon:
        plot(x, eta, label="t=4s")

xlabel("x(meters)")
ylabel("eta")
plt.title("Water surface using FTCS method")
plt.legend()
show()
