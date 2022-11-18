import pprint
import numpy as np
import scipy
import scipy.linalg   # SciPy Linear Algebra Library

""" n = int(input('Insira o numero de equacoes: '))

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
    B[i] = A[i][n] """

n = 4
A = np.array([[1, 1, 0, 3], [2, 1, -1, 1], [3, -1, -1, 2], [-1, 2, 3, -1]], dtype=np.float16)
B = np.array([4, 1, -3, 4], dtype=np.float16)
solucao = np.zeros(n)

P, L, U = scipy.linalg.lu(A)

print("A:")
pprint.pprint(A)

""" print("P:")
pprint.pprint(P) """

print("L:")
pprint.pprint(L)

print("U:")
pprint.pprint(U)

LU = np.zeros((n, n), dtype=np.float16)

for i in range(n):
    for j in range(n):
        if i>j:
            LU[i][j] = float(L[i][j])
        else:
            LU[i][j] = float(U[i][j])
            
print("LU:")
pprint.pprint(LU)


