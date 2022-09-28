import matplotlib.pyplot as plt
import numpy as np
import cmath as math


def cooling(T0, k, T_m, dt, npoints, type): # Explicit Euler, Implicit Euler and Crank-Nicolson
    k *= -1
    T = np.zeros(npoints)
    T[0] = T0
    for i in range(npoints-1):
        if type == 1:
            T[i+1] = dt*(k*(T[i] - T_m)) + T[i] # Explicit Euler
        
        if type == 2:
            T[i+1] = (T[i] - k * dt * T_m)/(1-k*dt) # Implicit Euler
        
        if(type == 3): # Crank-Nicolson
            T[i+1] = T[i] + dt*k*(T[i] + (1/2)*dt*k*(T[i] - T_m) - T_m) # Crank-Nicolson
        
    return T


def main():
    T0 = 99
    k = 0.0358
    T_m = 27
    t_ini = 0
    t_end = 100
    npoints = 10

    t = np.linspace(t_ini,t_end,npoints)
    dt = (t_end-t_ini)/(npoints-1)
    
    T= cooling(T0, k, T_m, dt, npoints, 1)
    plt.plot(t, T, '-*')
    T = cooling(T0, k, T_m, dt, npoints, 2)
    plt.plot(t, T, '-o')
    T = cooling(T0, k, T_m, dt, npoints, 3)
    plt.plot(t, T, '-^')

    tt = np.linspace(t_ini, t_end, 200)
    plt.plot(tt, ((T0 - T_m)*np.exp(-k*tt)) + T_m, '--')
    plt.legend(['Explicit Euler', 'Implicit Euler', 'Crank-Nicolson', 'Exact'])

    plt.xlabel('tempo (s)')
    plt.ylabel('temperatura (C)')
    plt.title('Resfriamento de Newton')
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
