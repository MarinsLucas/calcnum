
    # Inicializa a matriz L e U
    L = np.zeros((n, n), dtype=np.float16)
    U = np.zeros((n, n), dtype=np.float16)

    # Cria a matriz pivô P e a matriz multiplicada PA
    P = matriz_pivo(A)
    PA = matriz_multi(P, A)

    # Decomposição LU
    for j in range(n):
        # Os termos da diagonal de L são 1
        L[j][j] = 1.0

        # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}
        for i in range(j + 1):
            s1 = np.float16(sum(U[k][j] * L[i][k] for k in range(i)))
            U[i][j] = np.float16(PA[i][j] - s1)

        # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik})
        for i in range(j, n):
            s2 = np.float16(sum(U[k][j] * L[i][k] for k in range(j)))
            L[i][j] = np.float16((PA[i][j] - s2) / U[j][j])

    
    pprint.pprint(L)
    pprint.pprint(U)


    # Substituição regressiva
    x[n - 1] = np.float16(L[n - 1][n]/L[n - 1][n - 1])

    for i in range(n - 2, -1, -1):
        x[i] = np.float16(L[i][n])

        for j in range(i + 1, n):
            x[i] = np.float16(x[i] - L[i][j] * x[j])

        x[i] = np.float16(x[i]/L[i][i])

    print('Matriz L:')
    pprint.pprint(L)
    
    print('Matriz U:')
    pprint.pprint(U)

    return x




def cholesky(matriz, n, solucao):
    # Decomposição de Cholesky
    for i in range(n):
        for j in range(i + 1):
            soma = sum(matriz[i][k] * matriz[j][k] for k in range(j))

            if (i == j):  # se for na diagonal
                if (matriz[i][i] - soma <= 0):  # verifica se é positiva definida
                    sys.exit('Matriz nao e positiva definida!')

                matriz[i][i] = np.sqrt(matriz[i][i] - soma)
            else:
                matriz[i][j] = np.float16((matriz[i][j] - soma)/matriz[j][j])

    # Substituição progressiva
    for i in range(n):
        if (matriz[i][i] == 0.0):
            sys.exit('Divisao por zero detectada!')

        for j in range(i + 1, n):
            razao = np.float16(matriz[j][i]/matriz[i][i])

            for k in range(n + 1):
                matriz[j][k] = np.float16(matriz[j][k] - razao * matriz[i][k])

    # Substituição regressiva
    solucao[n - 1] = np.float16(matriz[n - 1][n]/matriz[n - 1][n - 1])

    for i in range(n - 2, -1, -1):
        solucao[i] = np.float16(matriz[i][n])

        for j in range(i + 1, n):
            solucao[i] = np.float16(solucao[i] - matriz[i][j] * solucao[j])

        solucao[i] = np.float16(solucao[i]/matriz[i][i])

    return solucao