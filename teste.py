'''
Formulate a θ-rule for the three schemes (Forward Euler, Backward Euler and Crank-Nicolson) for Newton's law of cooling such that you can get the three schemes from a single formula by varying the θ parameter. Implement the θ scheme in a function cooling(T0, k, T_s, t_end, dt, theta=0.5), where T0 is the initial temperature, k is the heat transfer coefficient, T_s is the temperature of the surroundings, t_end is the end time of the simulation, dt is the time step, and theta corresponds to θ. The cooling function should return the temperature as an array T of values at the mesh points and the time mesh t. Construct verification examples to check that the implementation works.

For verification, use the exact solution T(t) = (T0 - T_s) * e^(-kt) + T_s.
Plot one graph for the approximations and the exact solution using pyplot.

Filename: cooling.py.
'''

import numpy as np
from matplotlib import pyplot as plt

def cooling(T0, k, T_s, t_end, dt, theta=0.5):
    t = np.arange(0, t_end, dt)
    N_t = len(t)
    T = np.zeros(N_t)
    T[0] = T0
    for n in range(N_t-1):
        if theta == 0:
            T[n+1] = T[n] - dt*k*(T[n] - T_s)
        elif theta == 1:
            T[n+1] = T[n] + dt*k*(T_s - T[n])
        else:
            T[n+1] = T[n] + theta*dt*k*(T_s - T[n]) + (1-theta)*dt*k*(T[n] - T_s)
    return T, t

T0 = 50
k = 0.1
T_s = 20
t_end = 10
dt = 0.1

T, t = cooling(T0, k, T_s, t_end, dt)
T_exact = (T0 - T_s)*np.exp(-k*t) + T_s

plt.plot(t, T, 'r-o', label='dT/dt = k*(T_s - T)')
plt.plot(t, T_exact, 'b-o', label='exact solution')
plt.xlabel('t')
plt.ylabel('T')
plt.legend()
plt.title('T(t) for k = 0.1, T0 = 50, T_s = 20')
plt.show()