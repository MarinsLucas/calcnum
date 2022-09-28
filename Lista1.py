import matplotlib.pyplot as plt
import numpy as np
import cmath as math


def cooling(T0, k, T_s, t_end, dt, theta, npoints):
    T = np.zeros(npoints)
    T[0] = T0
    for i in range(npoints-1):
        T[i+1] = (1.0 - k*dt*(1.0-theta))/(1.0 + k*dt*theta) * \
            T[i] + (k*dt*(1.0-theta)*T_s)/(1.0 + k*dt*theta)
    return T


def main():
    T0 = 99
    k = 0.0358
    T_s = 27
    t_ini = 0
    t_end = 100
    refin = 5

    for it in range(0, refin):
        npoints = 4**(it+1)+1
        dt = (t_end-t_ini)/(npoints-1)
        ''' T, t = cooling(T0, k, T_s, t_end, dt, 0, npoints)
        plt.plot(t, T, '-*')
        T, t = cooling(T0, k, T_s, t_end, dt, 0.5, npoints)
        plt.plot(t, T, '-o') '''
        t = np.linspace(t_ini,t_end,npoints)
        T = cooling(T0, k, T_s, t_end, dt, 1, npoints)
        plt.plot(t, T, '-^')

        tt = np.linspace(t_ini, t_end, 200)
        plt.plot(tt, ((T0 - T_s)*np.exp(-k*tt)) + T_s, '--')
        plt.legend(['Forwards Euller', 'Exact'])

        plt.xlabel('tempo (s)')
        plt.ylabel('temperatura (C)')
        plt.title('Resfriamento de Newton')
        plt.grid()
        plt.show()


if __name__ == "__main__":
    main()
