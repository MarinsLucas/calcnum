from telnetlib import TM
import matplotlib.pyplot as plt
import numpy as np
import cmath as math

def create_plot(ptype):
    x = np.arange(0, 500, dt)
    if ptype == 'resfriamentoNewtonExata':
        y = ((thetaInicial - thetaM)*np.exp(-k*x)) + thetaM        
    return(x, y)

thetaInicial = 99
thetaM = 27
dt = 0.01
k = 0.035871952
x,y = create_plot('resfriamentoNewtonExata')

u = np.zeros(x.size, np.int16)
u[0] = y[0]
for n in range(0, x.size-1):
    u[n+1] = u[n]*(1-dt*u[0])

plt.plot(x, u)
plt.plot(x, y)
plt.show()