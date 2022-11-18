import pprint
import numpy as np
import scipy
import scipy.linalg   # SciPy Linear Algebra Library

def decompoeLU(A, n):
    for j in range(n):
        for i in range(n):
            # Parte L da matriz
            if(i >= j):
                s = 0;
                for k in range(0, j):
                    s += A[i][k] * A[k][j]
                A[i][j] -= s 
                
            # Parte U da matriz
            else: 
                s = 0
                for k in range(0, i):
                    s += A[i][k] * A[k][j]
                    
                A[i][j] = (A[i][j] - s) / float(A[i][i])
                
    return A

def resolveLU(LU, B, n):
    
    # Resolve Ly = b
    y = np.zeros(n, dtype=float)
    
    for i in range(n):
        s = 0
        for j in range(0, i):
            s += LU[i][j] * y[j]
        y[i] = (B[i] - s) / float(LU[i][i])
        
    # Resolve Ux = y
    x = np.zeros(n, dtype=float)
    
    for i in range(n-1, -1, -1):
        s = 0
        for j in range(i+1, n):
            s += LU[i][j] * x[j]
        x[i] = (y[i] - s)
    
    return x

def lu(A, B, n):
    A = decompoeLU(A, n)
    return resolveLU(A, B, n)


def main():
    n = 3
    A = np.array([[2, 4, -2], [4, 9, -3], [-2, -3, 7]], dtype=float)
    B = np.array([2, 8, 10], dtype=float)
    solucao = np.zeros(n)

    A_aumentada = np.zeros((n, n + 1), dtype=float)
    for i in range(n):
        for j in range(n):
            A_aumentada[i][j] = A[i][j]
        A_aumentada[i][n] = B[i]
    
    solucao = lu(A, B, n)
    pprint.pprint(solucao)
    
    
    
if __name__ == "__main__":
    main()





