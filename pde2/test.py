import matplotlib.pyplot as plt
import numpy as np

L = 1e-8
m = 9.109e-31
sigma = L / 25
kappa = 500 / L
x0 = -L / 2
xP = L / 2
P = 1024
tau = 1e-18
x_0 = L / 5
N = 300
a = L / P
hbar = 6.626e-34
A = -((hbar) ** 2) / (2 * m * a ** 2)
u = np.ones(P - 1)
I = np.diag(u, k=0)
omega = 3e15


def gaussian(x):
    return (np.exp(-((x - x_0) / (2 * sigma)) ** 2)) * (
        np.exp(complex(0, kappa * x)))


def gaussiansquared(x):
    return np.exp(-2 * ((x - x_0) / (2 * sigma)) ** 2)


def norm(x):
    d = 0
    for i in range(0, P - 1):
        d += (np.conj(x[i])) * x[i]
    psi0 = (1 / (d * a)) ** 0.5
    return psi0 * x


def normconstant(x):
    d = 0
    for i in range(0, P - 1):
        d += (np.conj(x[i])) * x[i]
    psi0 = (1 / (d * a)) ** 0.5
    return np.real(psi0)


def trapezoidal(x, y):
    d = 0
    for i in range(0, P - 1):
        d += x[i] * y[i]
    return d * a


xpoints = np.zeros(P - 1)
for i in range(0, P - 1):
    xpoints[i] = x0 + (i + 1) * a

tpoints = np.zeros(N)
for i in range(N):
    tpoints[i] = i * tau

vec_diag1 = np.ones(P - 1, complex)
for i in range(len(vec_diag1)):
    vec_diag1[i] = -2 * A

D1 = np.diag(vec_diag1, k=0)
Sup = A * (np.eye(P - 1, k=1))
Sub = A * (np.eye(P - 1, k=-1))
H1 = D1 + Sub + Sup
R1 = I - (complex(0, tau / (2 * hbar)) * H1)
M1 = I + (complex(0, tau / (2 * hbar)) * H1)

PSI = []
Normalization = np.zeros(N)
for i in range(N):
    PSI.append([])

for m in range(N):
    if m == 0:
        psi = np.zeros(P - 1, complex)
        for i in range(P - 1):
            psi[i] = gaussian(xpoints[i])

        s = 0
        for i in range(0, P - 1):
            s += gaussiansquared(x0 + (i + 1) * a)
        psi0 = (1 / (s * a)) ** 0.5
        PSI[m] = psi0 * psi
        Normalization[m] = psi0
    else:
        V = np.dot(R1, PSI[m - 1])
        X = np.linalg.solve(M1, V)
        PSI[m] = norm(X)
        Normalization[m] = normconstant(X)


def psi_real(psi):
    psi_real_ = []
    position = [- P / 2 * a]
    for x_index in range(len(psi)):
        psi_real_.append(psi[x_index].real)
        if x_index < (len(psi) - 1):
            position.append(position[-1] + a)
    return [psi_real_, position]


def plot_psi_real(psi_real_, position):
    plt.figure()
    plt.plot(position, psi_real_)
    plt.xlabel("X (m)")
    plt.ylabel("Real part of $\Psi$")
    plt.title("Real part of $\Psi$ at time=0")
    plt.show()


def plot_psi_real_timeline(psi_timeline_):
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
    plt.title("Real part of $\Psi$ at different times")
    plt.show()


plot_psi_real_timeline(PSI)
