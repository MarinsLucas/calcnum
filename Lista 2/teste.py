import pprint
import numpy as np
import scipy
import scipy.linalg   # SciPy Linear Algebra Library

n = int(input('Insira o numero de equacoes: '))

# Initialize A as a n x n scipy.array of zeros
A = np.zeros((n, n + 1), dtype=np.float16)
B = np.zeros((n), dtype=np.float16)


for i in range(n):
    for j in range(n):
        A[i][j] = np.float16(1/(i + j + 1))

for i in range(n):
    A[i][n] = np.float16(1/(i + n + 1))

# B is the last column of A
for i in range(n):
    B[i] = A[i][n]

P, L, U = scipy.linalg.lu(A)

print("A:")
pprint.pprint(A)

print("P:")
pprint.pprint(P)

print("L:")
pprint.pprint(L)

print("U:")
pprint.pprint(U)

V = np.zeros((n, n), dtype=np.float16)

for i in range(n):
    for j in range(n):
        if i>j:
            V[i][j] = float(L[i][j])
        else:
            V[i][j] = float(U[i][j])
            
print("V:")
pprint.pprint(V)

#Ly=b
y = np.zeros(n, dtype=float)
for i in range(n):
    s = 0
    for j in range(0, i):
        s += V[i][j]*y[j]
    y[i] = (B[i] - s)/ float(A[i][i])

#Ux=y
x = np.zeros(n, dtype=float)
for i in range(n-1, -1 , -1):
    s = 0
    for j in range(i+1, n):
        s+= V[i][j]*x[j]
    x[i] = y[i] - s

print("x:")
pprint.pprint(x)
