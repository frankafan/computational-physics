import matplotlib.pyplot as plt

SIGMA = 1.0
EPSILON = 1.0
M = 1.0
V1_0X = 0.0
V1_0Y = 0.0
V2_0X = 0.0
V2_0Y = 0.0
STEPS = 100
H = 0.01

trajectory_names = [
    'Particle pair 1',
    'Particle pair 2',
    'Particle pair 3',
]


def a_1x(r1, r2):
    first_term = - 2 * (r2[0] - r1[0]) / (
            (r2[0] - r1[0]) ** 2 + (r2[1] - r2[1]) ** 2) ** 7
    second_term = (r2[0] - r1[0]) / (
            (r2[0] - r1[0]) ** 2 + (r2[1] - r2[1]) ** 2) ** 4
    return 12 * (first_term + second_term)


def a_2x(r1, r2):
    return - a_1x(r1, r2)


def a_1y(r1, r2):
    first_term = - 2 * (r2[1] - r1[1]) / (
            (r2[0] - r1[0]) ** 2 + (r2[1] - r2[1]) ** 2) ** 7
    second_term = (r2[1] - r1[1]) / (
            (r2[0] - r1[0]) ** 2 + (r2[1] - r2[1]) ** 2) ** 4
    return 12 * (first_term + second_term)


def a_2y(r1, r2):
    return - a_1y(r1, r2)


for i in range(len(trajectory_names)):
    if i == 0:
        r1 = [4.0, 4.0]
        r2 = [5.2, 4.0]
    elif i == 1:
        r1 = [4.4, 4.0]
        r2 = [5.2, 4.0]
    elif i == 2:
        r1 = [2.0, 3.0]
        r2 = [3.5, 4.4]
    trajectory_name = trajectory_names[i]

    t = []
    v1_x = [V1_0X]
    v1_y = [V1_0Y]
    v2_x = [V2_0X]
    v2_y = [V2_0Y]
    e_k = []
    e_p = []
    e_t = []

    trajectory_1x = [r1[0]]
    trajectory_1y = [r1[1]]
    trajectory_2x = [r2[0]]
    trajectory_2y = [r2[1]]

    for i in range(STEPS):
        if i == 0:
            v1_xf = v1_x[0] + 0.5 * H * a_1x(r1, r2)
            v1_yf = v1_y[0] + 0.5 * H * a_1y(r1, r2)
            v2_xf = v2_x[0] + 0.5 * H * a_2x(r1, r2)
            v2_yf = v2_y[0] + 0.5 * H * a_2y(r1, r2)

            r = ((trajectory_1x[0] - trajectory_2x[0]) ** 2 + (
                        trajectory_1y[0] - trajectory_2y[0]) ** 2) ** 0.5
            ep = 4 * EPSILON * ((SIGMA / r) ** 12 - (SIGMA / r) ** 6)
            ek = 1 * M * (v1_xf ** 2 + v1_yf ** 2)
            et = ep + ek

            v1_x.append(v1_xf)
            v1_y.append(v1_yf)
            v2_x.append(v2_xf)
            v2_y.append(v2_yf)
            e_k.append(ek)
            e_p.append(ep)
            e_t.append(et)
            t.append(i)
        else:
            r1[0] += H * v1_x[-1]
            r1[1] += H * v1_y[-1]
            r2[0] += H * v2_x[-1]
            r2[1] += H * v2_y[-1]
            trajectory_1x.append(r1[0])
            trajectory_2x.append(r2[0])
            trajectory_1y.append(r1[1])
            trajectory_2y.append(r2[1])

            k_1x = H * a_1x(r1, r2)
            k_1y = H * a_1y(r1, r2)
            k_2x = H * a_2x(r1, r2)
            k_2y = H * a_2y(r1, r2)

            v1_xf = v1_x[-1] + k_1x
            v1_yf = v1_y[-1] + k_1y
            v2_xf = v2_x[-1] + k_2x
            v2_yf = v2_y[-1] + k_2y

            r = ((trajectory_1x[-1] - trajectory_2x[-1]) ** 2 + (
                    trajectory_1y[-1] - trajectory_2y[-1]) ** 2) ** 0.5
            ep = 4 * EPSILON * ((SIGMA / r) ** 12 - (SIGMA / r) ** 6)
            ek = 1 * M * (v1_xf ** 2 + v1_yf ** 2)
            et = ep + ek

            v1_x.append(v1_xf)
            v1_y.append(v1_yf)
            v2_x.append(v2_xf)
            v2_y.append(v2_yf)
            e_k.append(ek)
            e_p.append(ep)
            e_t.append(et)
            t.append(i)

    plt.figure()
    plt.plot(trajectory_1x, trajectory_1y, '.', label='Particle 1')
    plt.plot(trajectory_2x, trajectory_2y, '.', label='Particle 2')
    plt.legend()
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.title(trajectory_name + ' trajectory')
    plt.savefig(trajectory_name + ' trajectory')
    # plt.show()

    plt.figure()
    plt.plot(t, trajectory_1x, '.', label='Particle 1')
    plt.plot(t, trajectory_2x, '.', label='Particle 2')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('X coordinate')
    plt.title(trajectory_name + ' x-dimensional trajectory')
    plt.savefig(trajectory_name + ' x-dimensional trajectory')
    # plt.show()

    plt.figure()
    plt.plot(t, e_k, label='Kinetic energy')
    plt.plot(t, e_p, label='Potential energy')
    # plt.plot(t, e_t, label='Total energy')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (J)')
    plt.title(trajectory_name + ' energy')
    plt.savefig(trajectory_name + ' energy')
    # plt.show()
