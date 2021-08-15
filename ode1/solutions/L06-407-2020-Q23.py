# Author: Nico Grisouard, Univ. of Toronto

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

ftsz = 16
font = {'family': 'normal', 'size': ftsz}  # font size
rc('font', **font)


def init_pos(nmolecules): return nmolecules, np.empty((2, nmolecules))


def LJ_pot(pstns, periodic=False):
    """ Computes the Lennard-Jones potential at all locations for non-periodic
    boundary conditions.
    NOTE: too similar to rhs, should be modularized.
    """

    s = 1.  # sigma
    m = 1.  # mass
    eps = 1.  # epsilon

    NP = pstns.shape[1]  # number of molecules
    range_prtcls = list(range(NP))
    V = np.zeros(NP)

    if periodic:
        # replicate particle positions
        repeat_pstns = np.empty((2, 9*NP))  # 9 tiles => 8 replicas
        til = 0  # the tile index
        for xx in [0., -Lx, Lx]:  # starting from zero ensure the first tile
            for yy in [0., -Ly, Ly]:  # contains the original positions
                repeat_pstns[0, til*NP:(til+1)*NP] = pstns[0, :] + xx
                repeat_pstns[1, til*NP:(til+1)*NP] = pstns[1, :] + yy
                til += 1  # next tile
        range_other_prtcls = list(range(9*NP))  # accelerating molecules
    else:
        repeat_pstns = 1*pstns

    for n in range_prtcls:  # loop over molecules
        if periodic:
            other_molecules = range_prtcls[:n] + range_other_prtcls[n+1:]
        else:
            other_molecules = range_prtcls[:n] + range_prtcls[n+1:]

        for op in other_molecules:
            dist = pstns[:, n] - repeat_pstns[:, op]
            r = (dist[0]**2 + dist[1]**2)**.5
            V[n] += 4*eps/m*(s**12/r**12 - s**6/r**6)

    return V


def rhs(pstns, periodic=False):
    """ The right-hand-side of the velocity equations
    INPUT:
    pstns = [[x1, y1], [x2, y2], ...]  2*N array of floats
    note: no explicit dependence on time
    OUTPUT:
    [[Fx1, Fy1], [Fx2, Fy2], ...]  component of the sum of forces on each
    particle """
    # Keeping the quantities just in case
    s = 1.  # sigma
    m = 1.  # mass
    eps = 1.  # epsilon

    NP = pstns.shape[1]  # number of molecules
    range_prtcls = list(range(NP))
    F = 0*pstns

    if periodic:
        # replicate particle positions
        repeat_pstns = np.empty((2, 9*NP))  # 9 tiles => 8 replicas
        til = 0  # the tile index
        for xx in [0., -Lx, Lx]:  # starting from zero ensure the first tile
            for yy in [0., -Ly, Ly]:  # contains the original positions
                repeat_pstns[0, til*NP:(til+1)*NP] = pstns[0, :] + xx
                repeat_pstns[1, til*NP:(til+1)*NP] = pstns[1, :] + yy
                til += 1  # next tile
        range_other_prtcls = list(range(9*NP))  # accelerating molecules
    else:
        repeat_pstns = 1*pstns

    for n in range_prtcls:  # loop over molecules
        if periodic:
            other_molecules = range_prtcls[:n] + range_other_prtcls[n+1:]
        else:
            other_molecules = range_prtcls[:n] + range_prtcls[n+1:]

        for op in other_molecules:
            dist = pstns[:, n] - repeat_pstns[:, op]
            r = (dist[0]**2 + dist[1]**2)**.5
            # print(r**14)
            PreFac = 4*eps/m*(12*s**12/r**14 - 6*s**6/r**8)
            F[:, n] += PreFac * dist

    return F


dt = 0.01  # time step
case = 'v'  # one case for each set of initial conditions

if case == 'i':  # Q2b
    Nt = 100  # number of time steps
    N, rs = init_pos(2)  # number of molecules, inital positions
    periodic = False  # True or False
    rs[:, 0] = np.array([4., 4.])
    rs[:, 1] = np.array([5.2, 4.])
elif case == 'ii':  # Q2b
    Nt = 100  # number of time steps
    N, rs = init_pos(2)  # number of molecules, initial positions
    periodic = False  # True or False
    rs[:, 0] = np.array([4.5, 4.])
    rs[:, 1] = np.array([5.2, 4.])
elif case == 'iii':  # Q2b
    Nt = 100  # number of time steps
    N, rs = init_pos(2)  # number of molecules, initial positions
    periodic = False  # True or False
    rs[:, 0] = np.array([2., 3.])
    rs[:, 1] = np.array([3.4, 4.4])
elif case == 'iv':  # Q3a
    Nt = 1000  # number of time steps
    N, rs = init_pos(16)
    periodic = False  # True or False
    Lx = 4.
    Ly = 4.
    dx = Lx/np.sqrt(N)
    dy = Ly/np.sqrt(N)
    x_grid = np.arange(.5*dx, Lx, dx)
    y_grid = np.arange(.5*dy, Ly, dy)
    xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    rs[0, :] = xx_grid.flatten()
    rs[1, :] = yy_grid.flatten()
elif case == 'v':  # Q3b
    Nt = 350  # number of time steps
    periodic = True  # True or False
    N, rs = init_pos(16)
    Lx = 4.
    Ly = 4.
    dx = Lx/np.sqrt(N)
    dy = Ly/np.sqrt(N)
    x_grid = np.arange(.5*dx, Lx, dx)
    y_grid = np.arange(.5*dy, Ly, dy)
    xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    rs[0, :] = xx_grid.flatten()
    rs[1, :] = yy_grid.flatten()
else:
    raise NameError('Which case of initial positions are we investigating?')

if rs.shape[1] != N:
    raise NameError('rs does not have the correct number of molecules')

pos = np.empty((2, Nt, N))  # 1st index: x, y
vel = np.empty((2, Nt, N))  # idem

# Jump-start with RK2
v0 = np.zeros((2, N))  # initial velocity components for each particle

# Below: need to make sure that rs and pos[:, i, :] have same shape
pos[:, 0, :] = np.array([rs[0], rs[1]], float)
vel[:, 0, :] = v0[:, :] + .5*dt*rhs(pos[:, 0, :], periodic)

KE = np.empty(Nt)  # this array will contain the Kin. energy w.r.t. time.
PE = np.empty(Nt)

KE[0] = 0.  # zero initial velocities
PE[0] = sum(LJ_pot(pos[:, 0, :], periodic)/2)

for tt in range(1, Nt):
    pos[:, tt, :] = pos[:, tt-1, :] + dt*vel[:, tt-1, :]
    if periodic:
        for pp in range(N):
            pos[0, tt, pp] = np.mod(pos[0, tt, pp], Lx)
            pos[1, tt, pp] = np.mod(pos[1, tt, pp], Ly)
    k = dt*rhs(pos[:, tt, :], periodic)
    vel[:, tt, :] = vel[:, tt-1, :] + k

    PE[tt] = sum(LJ_pot(pos[:, tt, :], periodic)/2)
    # above: note the division by two bc V is the PE of a couple of molecs.
    v_int = vel[:, tt-1, :] + .5*k
    KE[tt] = sum(0.5*(v_int[0, :]**2 + v_int[1, :]**2))

E = KE + PE  # total energy

plt.figure()
for n in range(N):
    plt.plot(pos[0, :, n], pos[1, :, n], '.', color='C{0}'.format(n % 10))
    # Above: 'C0', 'C1', ...,. C9 are the default colour lines.
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title('Trajectory of {0} molecules, case {1}.'.format(N, case),
          fontsize=ftsz)
plt.axis('equal')
plt.grid()
plt.tight_layout()
# plt.savefig('pos_{0}parts_case_{1}.png'.format(N, case), dpi=150)

plt.figure()
for n in range(N):
    plt.plot(pos[0, :, n], '.', color='C{0}'.format(n % 10))
    # Above: 'C0', 'C1', ...,. C9 are the default colour lines.
plt.xlabel("$t$ (time step)")
plt.ylabel("$x$")
plt.title('$x$-coordinates of {0} molecules, case {1}.'.format(N, case),
          fontsize=ftsz)
plt.grid()
plt.tight_layout()
# plt.savefig('x_{0}parts_case_{1}.png'.format(N, case), dpi=150)

if case == 'iv':

    print(max(E))
    print(min(E))
    print(np.mean(E))
    print('Ratio (max(E)-min(E))/mean(E) = ', (max(E)-min(E))/np.mean(E))

    plt.figure(figsize=(8, 6))
    plt.plot(E, label='total energy')
    plt.plot(KE, label='kinetic energy')
    plt.plot(PE, label='potential energy')
    plt.xlabel("$t$ (time steps)")
    plt.ylabel("$E(t)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    # plt.savefig('E_{0}parts_case_{1}.png'.format(N, case), dpi=150)

plt.show()
