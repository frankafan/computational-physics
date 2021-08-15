import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h

L = 1e-8
M = 9.109e-31
SIGMA = L / 25
K = 500 / L
P = 1024  # number of cells in the x direction
a = L / P  # spatial step size
TAU = 1e-18  # time step size
ACCURACY = 10e-6

X0 = L / 5
N = 3000  # number of time steps
T = N * TAU  # duration

# constants for potential calculation
OMEGA = 3e15
V0 = 6e-17
X1 = L / 4

SAVEFIG = True
LOADOUPUT = False


def V(x):
    """Returns potential at point x"""
    # return V0 * (x ** 2 / X1 ** 2 - 1) ** 2
    # return 0.5 * M * OMEGA ** 2 * x ** 2
    return 0


def A():
    """Returns matrix A for the Hamiltonian"""
    return - h ** 2 / (2 * M * a ** 2)


def B(p):
    """Returns matrix B for the Hamiltonian"""
    return V(p * a - L / 2) - 2 * A()


def HD():
    """Returns the Hamiltonian"""
    H_ = np.zeros((P - 1, P - 1))
    for x in range(P - 1):
        for y in range(P - 1):
            if x == y:
                H_[x][y] = B(x)
            elif abs(x - y) == 1:
                H_[x][y] = A()
            else:
                H_[x][y] = 0
    return H_


H = HD()


def L_matrix():
    """Returns matrix L for the Crank-Nicolson scheme"""
    return np.identity(P - 1) + 1j * TAU / (2 * h) * H


def R_matrix():
    """Returns matrix R for the Crank-Nicolson scheme"""
    return np.identity(P - 1) - 1j * TAU / (2 * h) * H


def psi0_at_x(x, psi00_=1 / 10):
    """Returns the wave function at point x at time 0"""
    return np.exp(- (x - X0) ** 2 / (4 * SIGMA ** 2) + 1j * K * x) * psi00_


def psi0(psi00_=1 / 10):
    """Returns the wave function during at time 0"""
    psi_x = []
    for x_index in range(int(- P / 2) + 1, int(P / 2)):
        psi_at_x = psi0_at_x(x_index * a, psi00_)
        psi_x.append(psi_at_x)
    return psi_x


def psi_next(psi_now):
    """Returns the next wave function given the current wave function"""
    v = np.matmul(R_matrix(), psi_now)
    return np.linalg.solve(L_matrix(), v)


def energy(psi):
    """Returns the energy of the particle given the wave function"""
    return np.real(np.sum(np.matmul(np.matmul(np.conj(psi), H), psi)) * a)


def probability(psi):
    """Returns the normalization value, expected to be 1"""
    return np.real(np.matmul(np.conj(psi), psi))


def probability_density(psi):
    """Returns the probability density of the given wave function"""
    return (np.conj(psi)) * psi


def normalize_psi00(accuracy=10e-6):
    """Returns the normalization constant"""
    psi00_ = 1 / 10
    probability0 = 0
    while abs(1 - probability0) > accuracy:
        psi00_ -= accuracy
        probability0 = probability(psi0(psi00_))
    return psi00_


def psi_real(psi):
    """Returns a wave function with only real parts"""
    psi_real_ = []
    position = [- P / 2 * a]
    for x_index in range(len(psi)):
        psi_real_.append(psi[x_index].real)
        if x_index < (len(psi) - 1):
            position.append(position[-1] + a)
    return [psi_real_, position]


def plot_psi_real(psi_real_, position):
    """Plot real part of wave function at time 0"""
    plt.figure()
    plt.plot(position, psi_real_)
    plt.xlabel("X (m)")
    plt.ylabel("Real part of $\Psi$")
    plt.title("Real part of $\Psi$ at time=0")
    if SAVEFIG:
        plt.savefig("Real part of psi at time=0")


def psi_timeline(psi00_):
    """Return the array of wave functions over time"""
    psi_timeline_ = [psi0(psi00_)]
    for t_index in range(1, N):
        psi_timeline_.append(psi_next(psi_timeline_[-1]))
    return psi_timeline_


def plot_time_plots(psi_timeline_):
    """Plot energy and normalization"""
    time = []
    energy_timeline = []
    normalization_timeline = []
    for t_index in range(N):
        time.append(t_index * TAU)
        energy_timeline.append(energy(psi_timeline_[t_index]))
        normalization_timeline.append(probability(psi_timeline_[t_index]))

    plt.figure()
    plt.plot(time, normalization_timeline)
    plt.ylim((0, 2))
    plt.xlabel("Time (s)")
    plt.ylabel("Probability")
    plt.title("Electron probability across length L over time")
    if SAVEFIG:
        plt.savefig("Electron probability across length L over time")

    plt.figure()
    plt.plot(time, energy_timeline)
    plt.ylim((0, 1e-26))
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.title("Electron energy over time")
    if SAVEFIG:
        plt.savefig("Electron energy over time")


def plot_probability_density(psi_timeline_):
    """Plot probability densities over time"""
    probability_density0 = probability_density(psi_timeline_[0])
    probability_density1 = probability_density(
        psi_timeline_[int(len(psi_timeline_) / 4)])
    probability_density2 = probability_density(
        psi_timeline_[int(len(psi_timeline_) / 2)])
    probability_density3 = probability_density(psi_timeline_[-1])

    plt.figure()
    plt.plot(probability_density0, label="t = 0")
    plt.plot(probability_density1, label="t = T/4")
    plt.plot(probability_density2, label="t = T/2")
    plt.plot(probability_density3, label="t = T")
    plt.legend()
    plt.xlabel("X (m)")
    plt.ylabel("Probability")
    plt.title("Electron probability density across length L with square well")
    if SAVEFIG:
        plt.savefig("Electron probability density across length L with square well")


def plot_psi_real_timeline(psi_timeline_):
    """Plot the real part of wave function over time"""
    psi_real0 = psi_real(psi_timeline_[0])
    psi_real1 = psi_real(psi_timeline_[int(len(psi_timeline_) / 4)])
    psi_real2 = psi_real(psi_timeline_[int(len(psi_timeline_) / 2)])
    psi_real3 = psi_real(psi_timeline_[-1])
    plt.figure()
    plt.plot(psi_real0[1], psi_real0[0], label="t = 0")
    plt.plot(psi_real1[1], psi_real1[0], label="t = T/4")
    plt.plot(psi_real2[1], psi_real2[0], label="t = T/2")
    plt.plot(psi_real3[1], psi_real3[0], label="t = T")
    plt.legend()
    plt.xlabel("X (m)")
    plt.ylabel("Real part of $\Psi$")
    plt.title("Real part of $\Psi$ with square well at different times")
    if SAVEFIG:
        plt.savefig("Real part of psi at with square well different times")


# psi_real0_info = psi_real(psi0(psi00))
# plot_psi_real(psi_real0_info[0], psi_real0_info[1])


if LOADOUPUT:
    npzfile = np.load('psi_timeline.npz')
    psi_timeline = npzfile['psi']
else:
    psi00 = normalize_psi00()
    psi_timeline = psi_timeline(psi00)
    np.savez('psi_timeline', psi=psi_timeline)

# plot_time_plots(psi_timeline)
plot_probability_density(psi_timeline)
plot_psi_real_timeline(psi_timeline)

plt.show()
