import numpy as np
import matplotlib.pyplot as plt

def uexata(T0, k, T_s, t):
    return (T0 - T_s) * np.exp(-k * t) + T_s


def cooling(T0, k, T_s, t_end, dt, theta=0.5):
    npoints = int(t_end/dt)+1
    t = np.linspace(0, t_end, npoints)
    T = np.zeros(npoints)
    T[0] = T0

    for n in range(0, npoints-1):
        T[n+1] = (1.0-k*dt*(1.0-theta))/(1.0+k*dt*theta) * \
            T[n] + (k*dt*(1.0-theta)*T_s)/(1.0+k*dt*theta)

    return T, t


def main():
    T0 = 1.0
    k = 0.1
    T_s = 0.0
    t_end = 5.0
    dt = 0.1

    T, t = cooling(T0, k, T_s, t_end, dt, theta=0.5)
    tt = np.linspace(0, t_end, 200)
    plt.plot(t, T, '-*', tt, uexata(T0, k, T_s, tt))
    plt.legend(['aproximada', 'exata'])
    plt.show()
    T, t = cooling(T0, k, T_s, t_end, dt, theta=0.0)
    tt = np.linspace(0, t_end, 200)
    plt.plot(t, T, '-*', tt, uexata(T0, k, T_s, tt))
    plt.legend(['aproximada', 'exata'])
    plt.show()
    T, t = cooling(T0, k, T_s, t_end, dt, theta=1.0)
    tt = np.linspace(0, t_end, 200)
    plt.plot(t, T, '-*', tt, uexata(T0, k, T_s, tt))
    plt.legend(['aproximada', 'exata'])
    plt.show()


if __name__ == '__main__':
    main()
