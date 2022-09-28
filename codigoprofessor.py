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

print((np.log(erro[-1])-np.log(erro[0]))/(np.log(dtt[-1])-np.log(dtt[0])))
