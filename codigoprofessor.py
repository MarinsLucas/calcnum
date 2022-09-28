'''
Formulate a θ-rule for the three schemes (Forward Euler, Backward Euler and Crank-Nicolson) for Newton's law of cooling such that you can get the three schemes from a single formula by varying the θ parameter. Implement the θ scheme in a function cooling(T0, k, T_s, t_end, dt, theta=0.5), where T0 is the initial temperature, k is the heat transfer coefficient, T_s is the temperature of the surroundings, t_end is the end time of the simulation, dt is the time step, and theta corresponds to θ. The cooling function should return the temperature as an array T of values at the mesh points and the time mesh t. Construct verification examples to check that the implementation works.

For verification, use the exact solution T(t) = (T0 - T_s) * e^(-kt) + T_s.
Plot one graph for the approximations and the exact solution using pyplot.

Hint: Use the following code as basis.
'''
import numpy as np
import matplotlib.pyplot as plt

def uexata(a,b,t):
  return b*np.exp(-a*t)

ti = 0.0
tf = 1.0
refin = 5
erro = np.zeros(refin)
dtt = np.zeros(refin)
for k in range(0,refin):
  npoints = 4**(k+1)+1
  dt = (tf-ti)/(npoints-1)
  a = 1.0
  # condicao inicial
  b = 1.0
  # coordenadas
  t = np.linspace(ti,tf,npoints)
  u = np.zeros(npoints)
  u[0] = b

  for n in range(0,npoints-1):
    u[n+1] = (1.0/(1.0+a*dt))*u[n]
    #u[n+1] = ((2.-a*dt)/(2.0+a*dt))*u[n]

  tt = np.linspace(ti,tf,200)
  plt.plot(t,u,'-*',tt,uexata(a,b,tt))
  plt.legend(['aproximada','exata'])
  plt.show()
  erro[k] = np.max(np.abs(u-uexata(a,b,t)))
  dtt[k] = dt
  print(dt)
plt.plot(-np.log(dtt),np.log(erro))
plt.show()
'''

print((np.log(erro[-1])-np.log(erro[0]))/(np.log(dtt[-1])-np.log(dtt[0])))


Filename: cooling.py.
'''

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
    T[n+1] = (1.0-k*dt*(1.0-theta))/(1.0+k*dt*theta)*T[n] + (k*dt*(1.0-theta)*T_s)/(1.0+k*dt*theta)

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
  plt.legend(['aproximada','exata'])
  plt.show()
  T, t = cooling(T0, k, T_s, t_end, dt, theta=0.0)
  tt = np.linspace(0, t_end, 200)
  plt.plot(t, T, '-*', tt, uexata(T0, k, T_s, tt))
  plt.legend(['aproximada','exata'])
  plt.show()
  T, t = cooling(T0, k, T_s, t_end, dt, theta=1.0)
  tt = np.linspace(0, t_end, 200)
  plt.plot(t, T, '-*', tt, uexata(T0, k, T_s, tt))
  plt.legend(['aproximada','exata'])
  plt.show()

if __name__ == '__main__':
  main()