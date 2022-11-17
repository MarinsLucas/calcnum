import pprint
import numpy as np
import scipy
import scipy.linalg   # SciPy Linear Algebra Library

n = int(input('Insira o numero de equacoes: '))

# Initialize A as a n x n scipy.array of zeros
A = np.zeros((n, n + 1), dtype=np.float16)

for i in range(n):
    for j in range(n):
        A[i][j] = np.float16(1/(i + j + 1))

for i in range(n):
    A[i][n] = np.float16(1/(i + n + 1))


P, L, U = scipy.linalg.lu(A)

print("A:")
pprint.pprint(A)

print("P:")
pprint.pprint(P)

print("L:")
pprint.pprint(L)

print("U:")
pprint.pprint(U)
