import numpy as np
import matplotlib.pyplot as plt

n = 1089 #dimensao da matriz
k = int(np.sqrt(n)) # dimensao da matriz quadrada
A = np.zeros((n,n)) # matriz de zeros 
b = 10*np.ones(n) # vetor de n elementos com valor 10
u = np.zeros(n) # vetor de n elementos com valor 0

for i in range(0,n):
  A[i,i] = -4
for i in range(0,n-1):
  A[i+1,i] = 1
  A[i,i+1] = 1
  if((i+1)%k ==0):
    A[i+1,i] = 0
    A[i,i+1] = 0
for i in range(k,n):
  A[i-k,i] = 1
  A[i,i-k] = 1

A = - A*(k-1)**2
# print(A)
u = np.linalg.solve(A,b)

x = np.linspace(0,1,k)
y = np.linspace(0,1,k)

X,Y = np.meshgrid(x,y)

U = np.zeros((k,k))

for i in range(0,k):
  for j in range(0,k):
    U[i,j] = u[i+j*k]

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_surface(X, Y, U)
plt.show()