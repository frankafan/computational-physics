import numpy as np
import matplotlib.pyplot as plt
from SolveLinear import PartialPivot

w = 1000
r1, r3, r5 = 1e3, 1e3, 1e3
r2, r4, r6 = 2e3, 2e3, 2e3
c1, c2 = 1e-6, 0.5e-6
x_all = 3

M = [[complex(1 / r1 + 1 / r4, w * c1),
      complex(0, -w * c1),
      complex(0, 0)
      ],
     [complex(0, -w * c1),
      complex(1 / r2 + 1 / r5, w * c1 + w * c2),
      complex(0, -w * c2)
      ],
     [complex(0, 0),
      complex(0, -w * c2),
      complex(1 / r3 + 1 / r6, w * c2)
      ]]

b = [complex(x_all / r1, 0),
     complex(x_all / r2, 0),
     complex(x_all / r3, 0)
     ]

# x_vector = np.linalg.solve(M, b)
x_vector = PartialPivot(M, b)
x1 = x_vector[0]
x2 = x_vector[1]
x3 = x_vector[2]


def v(x, w, t):
    return x * np.exp(complex(0, w * t))


def t0():
    v1 = v(x1, w, 0)
    v2 = v(x2, w, 0)
    v3 = v(x3, w, 0)

    print("Amplitude of |V1|:", abs(v1))
    print("Amplitude of |V2|:", abs(v2))
    print("Amplitude of |V3|:", abs(v3))

    print("Phase of x1:", np.angle(v1) / np.pi * 180)
    print("Phase of x2:", np.angle(v2) / np.pi * 180)
    print("Phase of x3:", np.angle(v3) / np.pi * 180)


def v_plot():
    t_lst = []
    v1_lst = []
    v2_lst = []
    v3_lst = []
    v_all_lst = []

    for i in range(130):
        time = i * 0.0001
        t_lst.append(time)
        v1_lst.append(v(x1, w, time).real)
        v2_lst.append(v(x2, w, time).real)
        v3_lst.append(v(x3, w, time).real)

    plt.figure()
    plt.plot(t_lst, v1_lst, label='$V_1$')
    plt.plot(t_lst, v2_lst, label='$V_2$')
    plt.plot(t_lst, v3_lst, label='$V_3$')
    plt.ylim(-3, 3)
    plt.title("Voltage over time")
    plt.xlabel("time (s)")
    plt.ylabel("real voltage (V)")
    plt.legend()
    plt.savefig("Voltage over time")
    plt.show()


M_inductor = [[complex(1 / r1 + 1 / r4, w * c1),
               complex(0, -w * c1),
               complex(0, 0)
               ],
              [complex(0, -w * c1),
               complex(1 / r2 + 1 / r5, w * c1 + w * c2),
               complex(0, -w * c2)
               ],
              [complex(0, 0),
               complex(0, -w * c2),
               complex(1 / r3, w * c2) + 1 / complex(0, r6)
               ]]

# x_vector_inductor = np.linalg.solve(M_inductor, b)
x_vector_inductor = PartialPivot(M_inductor, b)
x1_inductor = x_vector_inductor[0]
x2_inductor = x_vector_inductor[1]
x3_inductor = x_vector_inductor[2]


def t0_inductor():
    v1_inductor = v(x1_inductor, w, 0)
    v2_inductor = v(x2_inductor, w, 0)
    v3_inductor = v(x3_inductor, w, 0)

    print("Amplitude of |V1|:", abs(v1_inductor))
    print("Amplitude of |V2|:", abs(v2_inductor))
    print("Amplitude of |V3|:", abs(v3_inductor))

    print("Phase of x1:", np.angle(v1_inductor) / np.pi * 180)
    print("Phase of x2:", np.angle(v2_inductor) / np.pi * 180)
    print("Phase of x3:", np.angle(v3_inductor) / np.pi * 180)


def v_plot_inductor():
    t_inductor_lst = []
    v1_inductor_lst = []
    v2_inductor_lst = []
    v3_inductor_lst = []

    for i in range(130):
        time = i * 0.0001
        t_inductor_lst.append(time)
        v1_inductor_lst.append(v(x1_inductor, w, time).real)
        v2_inductor_lst.append(v(x2_inductor, w, time).real)
        v3_inductor_lst.append(v(x3_inductor, w, time).real)

    plt.figure()
    plt.plot(t_inductor_lst, v1_inductor_lst, label='$V_1$')
    plt.plot(t_inductor_lst, v2_inductor_lst, label='$V_2$')
    plt.plot(t_inductor_lst, v3_inductor_lst, label='$V_3$')
    plt.ylim(-3, 3)
    plt.title("Voltage over time with 2H inductor as $R_6$")
    plt.xlabel("time (s)")
    plt.ylabel("real voltage (V)")
    plt.legend()
    plt.savefig("Voltage over time with 2H inductor as R6")
    plt.show()


t0()
v_plot()
t0_inductor()
v_plot_inductor()
