import numpy as np
import pprint
import sys

# -----------------------------------------------#
# Funções auxiliares


def troca_linha(matriz, i, j, n):
    for k in range(n + 1):
        temp = matriz[i][k]
        matriz[i][k] = matriz[j][k]
        matriz[j][k] = temp


def matriz_pivo(M):
    m = len(M)

    # Cria a matriz identidade com valores float16
    matriz_id = np.identity(m, dtype=np.float16)

    for k in range(m):
        # Inicializa o maior elemento para índice e pivô
        i_max = k
        pivo_max = M[i_max][k]

        # Encontra o maior elemento da coluna e associa ao pivô
        for i in range(k + 1, m):
            if (abs(M[i][k]) > pivo_max):
                pivo_max = M[i][k]
                i_max = i

        if not M[k][i_max]:  # Se o pivô for zero
            sys.exit('Divisao por zero detectada!')

        # Troca a linha do maior elemento com a linha atual
        if (i_max != k):
            troca_linha(matriz_id, k, i_max, m)

    return matriz_id


def matriz_multi(M, N):
    m = len(M)
    n = len(N[0])
    p = len(M[0])
    matriz = np.zeros((m, n), dtype=np.float16)

    for i in range(m):
        for j in range(n):
            for k in range(p):
                matriz[i][j] += M[i][k] * N[k][j]

    return matriz


def eliminacao_progressiva(M, n):
    for k in range(n):
        # Inicializa o maior elemento para índice e pivô
        i_max = k
        pivo_max = M[i_max][k]

        # Encontra o maior elemento da coluna e associa ao pivô
        for i in range(k + 1, n):
            if (abs(M[i][k]) > pivo_max):
                pivo_max = M[i][k]
                i_max = i

        if not M[k][i_max]:  # Se o pivô for zero
            sys.exit('Divisao por zero detectada!')

        # Troca a linha do maior elemento com a linha atual
        if (i_max != k):
            troca_linha(M, k, i_max, n)

        for i in range(k + 1, n):
            # fatora f para zerar os elementos abaixo do pivô
            f = np.float16(M[i][k]/M[k][k])

            # subtrai a linha do pivô multiplicada pelo fator f
            for j in range(k + 1, n + 1):
                M[i][j] = np.float16(M[i][j] - f * M[k][j])

            # preenche a matriz triangular inferior com zeros
            M[i][k] = 0

    return M
# -----------------------------------------------#

# Método de Gauss sem pivoteamento


def gauss(matriz, n, solucao):

    # Eliminação progressiva
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


# Método de Gauss com pivoteamento
def gauss_pivoteamento(matriz, n, solucao):

    # Eliminação progressiva
    for k in range(n):
        # Inicializa o maior elemento para índice e pivô
        i_max = k
        pivo_max = matriz[i_max][k]

        # Encontra o maior elemento da coluna e associa ao pivô
        for i in range(k + 1, n):
            if (abs(matriz[i][k]) > pivo_max):
                pivo_max = matriz[i][k]
                i_max = i

        if not matriz[k][i_max]:  # Se o pivô for zero
            sys.exit('Divisao por zero detectada!')

        # Troca a linha do maior elemento com a linha atual
        if (i_max != k):
            troca_linha(matriz, k, i_max, n)

        for i in range(k + 1, n):
            # fatora f para zerar os elementos abaixo do pivô
            f = np.float16(matriz[i][k]/matriz[k][k])

            # subtrai a linha do pivô multiplicada pelo fator f
            for j in range(k + 1, n + 1):
                matriz[i][j] = np.float16(matriz[i][j] - f * matriz[k][j])

            # preenche a matriz triangular inferior com zeros
            matriz[i][k] = 0

    # Substituição regressiva
    for i in range(n - 1, -1, -1):
        # Inicializa o vetor solução com a última coluna da matriz
        solucao[i] = matriz[i][n]

        # Subtrai os elementos da solução já encontrados
        # j é i + 1 pois os elementos da diagonal principal já foram zerados
        for j in range(i + 1, n):
            solucao[i] = np.float16(solucao[i] - matriz[i][j] * solucao[j])

        # Divide o elemento da solução pelo elemento da diagonal principal
        solucao[i] = np.float16(solucao[i]/matriz[i][i])

    return solucao

# Método da Decomposição LU


def decomposicao_LU(A, n, x):
    

# Método de Cholesky


def cholesky(matriz, n, solucao):
    return solucao


def main():

    # Lê o número de equações
    n = int(input('Insira o numero de equacoes: '))

    # Cria e inicializa a matriz ampliada
    matriz = np.zeros((n, n + 1), dtype=np.float16)

    # Cria e inicializa o vetor solução
    solucao = np.zeros(n, dtype=np.float16)

    # Insere os coeficientes na matriz ampliada
    for i in range(n):
        for j in range(n):
            matriz[i][j] = np.float16(1/(i + j + 1))

    for i in range(n):
        matriz[i][n] = np.float16(1/(i + n + 1))

    # Imprime a matriz ampliada
    print('\nMatriz ampliada: ')
    for i in range(n):
        for j in range(n + 1):
            print('%0.2f' % matriz[i][j], end='\t')
        print('')

    # Pergunta se o usuário quer usar o método de gauss com ou sem pivoteamento, e dependendo da resposta chama a função necessária
    metodo = int(input(
        '\n1 - Gauss sem pivoteamento\n2 - Gauss com pivoteamento\n3 - Decomposicao LU\n4 - Cholesky\n\nInsira o metodo que deseja usar: '))
    if metodo == 1:
        solucao = gauss(matriz, n, solucao)
    elif metodo == 2:
        solucao = gauss_pivoteamento(matriz, n, solucao)
    elif metodo == 3:
        solucao = decomposicao_LU(matriz, n, solucao)
    elif metodo == 4:
        solucao = cholesky(matriz, n, solucao)
    else:
        sys.exit('Metodo invalido!')

    # Imprime a solução
    """ print('\nSolucao: ')
    for i in range(n):
        print('X%d = %0.2f' % (i, solucao[i]), end='\t') """


if __name__ == "__main__":
    main()
