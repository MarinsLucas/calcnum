import numpy as np
import matplotlib.pyplot as plt


def create_matrix(n):
    A = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            A[i][j] = 1.0 / (float)(i + j + 1)

    return A


def main():
    arr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=int)
    n = 10
    det_n = np.zeros(n, dtype=float)
    det_n_growth = np.zeros(n, dtype=float)
    det_n_growth[0] = 1

    det_n_growth_growth = np.zeros(n, dtype=float)
    det_n_growth_growth[0] = 1

    for i in range(n):
        A = create_matrix(i)
        det_n[i] = np.linalg.det(A)

    for i in range(1, det_n.size):
        if (det_n[i-1] != 0):
            det_n_growth[i] = det_n[i] / det_n[i - 1]

    for i in range(1, det_n_growth.size):
        if (det_n_growth[i-1] != 0):
            det_n_growth_growth[i] = det_n_growth[i] / det_n_growth[i - 1]

    print('det_n: ', det_n)
    print('det_n_growth: ', det_n_growth)
    print('det_n_growth_growth: ', det_n_growth_growth)

    # plot a log scale graph for det_n with the "Valores de N" on the x axis and "Determinante da matriz NxN" on the y axis and the title "Relação entre o valor de N e o determinante da matriz NxN"
    plt.plot(arr, det_n)
    plt.yscale('log')
    plt.xlabel('Valores de N')
    plt.ylabel('Determinante da matriz NxN')
    plt.title('Relação entre o valor de N e o determinante da matriz NxN')
    plt.show()

    # plot a log scale graph for det_n_growth with the "Valores de N" on the x axis and "Determinante da matriz NxN" on the y axis and the title "Relação entre o valor de N e o determinante da matriz NxN"
    plt.plot(arr, det_n_growth)
    plt.yscale('log')
    plt.xlabel('Valores de N')
    plt.ylabel('Decrescimento da matriz NxN')
    plt.title('Relação de derescimento do determinante da matriz NxN')
    plt.show()


if __name__ == '__main__':
    main()
