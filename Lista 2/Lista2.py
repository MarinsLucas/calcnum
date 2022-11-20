import numpy as np
import pprint
import sys
import math
import matplotlib.pyplot as plt
import time


#TODO: - Calculo de erro
#TODO: (1) comparar o tempo de execução de todos para cada n em uma tabela
#TODO: (1) tabela de erros
#TODO: (1) conclusões sobre o determinante da matriz
#TODO: (2) tabela comparando os tempos de execução dos métodos da questão 1 e da 2 adotando direfentes valores de epslon e usando
#          n = 81, 289, 1089, 4225, 16641
# -----------------------------------------------#
# Funções auxiliares
def transposta(A, n):
    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            B[i][j] = A[j][i]
    return B


def troca_linha(A, i, j, n):
    for k in range(n + 1):
        temp = A[i][k]
        A[i][k] = A[j][k]
        A[j][k] = temp


def normaMaximo(x):
    size = len(x)
    maximo = abs(x[0])

    for i in range(size):
        temp = abs(x[i])
        if (temp > maximo):
            maximo = temp

    return maximo


def distanciaMaximo(x1, x2):
    if (len(x1) != len(x2)):
        print("O tamanho dos vetores x1 e x2 precisa ser o mesmo")
        return 0

    size = len(x1)
    dist = abs(x1[0] - x2[0])

    for i in range(size):
        temp = abs(x1[i] - x2[i])
        if (temp > dist):
            dist = temp

    return dist


def calculaErro(x_prox, x_atual):
    return distanciaMaximo(x_prox, x_atual) / normaMaximo(x_prox)


def substituicao_regressiva(M, B, n):
    if (len(M) == 0):
        return 0

    ordem = n
    temp = 0
    passos = 0

    #cria array de solucao
    sol = np.array([0] * ordem, np.float128)

    #checa se e superior
    isSuperior = (M[ordem - 1][0] == 0)

    #percorre as linhas da A
    for n in range(0, ordem):

        #define o i se for superior ou inferior
        if (isSuperior):
            i = ordem - n - 1
        else:
            i = n

        #valor do vetor solucao da linha correspondente
        temp = B[i]

        #soma os valores exceto o valor do pivo
        for j in range(0, ordem):
            if (j != i):
                temp += sol[j] * M[i][j] * -1
                passos += 1

        #calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i]

    return list(sol)


def substituicao_progressiva(M, B, n):
    if (len(M) == 0):
        return 0

    ordem = len(M[0])
    temp = 0
    passos = 0

    #cria array de solucao
    sol = np.array([0] * ordem, np.float128)

    #checa se e superior
    isSuperior = (M[ordem - 1][0] == 0)

    #percorre as linhas da A
    for n in range(0, ordem):

        #define o i se for superior ou inferior
        if (isSuperior):
            i = n
        else:
            i = ordem - n - 1

        #valor do vetor solucao da linha correspondente
        temp = B[i]

        #soma os valores exceto o valor do pivo
        for j in range(0, ordem):
            if (j != i):
                temp += sol[j] * M[i][j] * -1
                passos += 1

        #calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i]

    return sol


def matriz_aumentada(A, B, n):
    # Cria uma nova matriz com n linhas e n + 1 colunas e adiciona a coluna B à matriz A
    M = np.zeros((n, n + 1), dtype=np.float128)

    for i in range(n):
        for j in range(n):
            M[i][j] = A[i][j]

    for i in range(n):
        M[i][n] = B[i]

    return M


def decompoeLU(M, n):

    ordem = len(M[0])

    L = np.zeros((n, n), dtype=np.float128)
    U = np.zeros((n, n), dtype=np.float128)

    #print(type(M[0]))

    for j in range(ordem):
        U[0][j] = M[0][j]
    for i in range(ordem):
        L[i][0] = M[i][0] / U[0][0]

    for i in range(ordem):
        #Calcula L
        for j in range(i + 1):
            soma = 0.0
            for k in range(j):
                soma += L[i][k] * U[k][j]
            L[i][j] = M[i][j] - soma
            
        """ print("L: ")
        pprint.pprint(L) """

    #Calcula U
        for j in range(i, ordem):
            soma = 0.0
            for k in range(i):
                soma += L[i][k] * U[k][j]
            U[i][j] = (M[i][j] - soma) / L[i][i]
            
        """ print("U: ")
        pprint.pprint(U) """
    
    return (L, U)


def decomposicao_LU(M, B, n):
    if (len(M) == 0):
        return 0

    ordem = n

    y = [0.0 for i in range(ordem)]
    x = [0.0 for i in range(ordem)]
    (L, U) = decompoeLU(M, n)

    # Resolução L * y = b

    y[0] = B[0] / L[0][0]
    for i in range(1, ordem):
        soma = 0.0
        for j in range(i):
            soma += L[i][j] * y[j]
        y[i] = (B[i] - soma) / L[i][i]

    # Resolucao U * x = y
    x[ordem - 1] = y[ordem - 1] / U[ordem - 1][ordem - 1]
    for i in range(ordem - 1, -1, -1):
        soma = y[i]
        for j in range(i + 1, ordem):
            soma = soma - U[i][j] * x[j]
        x[i] = soma / U[i][i]
    return x


def erro(A, B, x, n):
    X = np.zeros(n)

    for i in range(n):
        for j in range(n):
            X[i] = abs(A[i][j] * x[i] - B[i])
    return X


def maxNorma(A, B, x, n):
    max = 0

    for i in range(n):
        for j in range(n):
            temp = abs(A[i][j] * x[i] - B[i])
            if (temp > max):
                max = temp

    return max


# Fim de funções auxiliares
# -----------------------------------------------#


# Método de Gauss sem pivoteamento
def gauss(A, B, n):
    solucao = np.zeros(n, dtype=np.float128)

    A = matriz_aumentada(A, B, n)

    # Substituição progressiva
    for i in range(n):
        if (A[i][i] == 0.0):
            sys.exit('Divisao por zero detectada!')

        for j in range(i + 1, n):
            razao = np.float128(A[j][i] / A[i][i])

            for k in range(n + 1):
                A[j][k] = np.float128(A[j][k] - razao * A[i][k])

    # Substituição regressiva
    solucao[n - 1] = np.float128(A[n - 1][n] / A[n - 1][n - 1])

    for i in range(n - 2, -1, -1):
        solucao[i] = np.float128(A[i][n])

        for j in range(i + 1, n):
            solucao[i] = np.float128(solucao[i] - A[i][j] * solucao[j])

        solucao[i] = np.float128(solucao[i] / A[i][i])

    return solucao


# Método de Gauss com pivoteamento
def gauss_pivoteamento(A, B, n):

    solucao = np.zeros(n, dtype=np.float128)

    A = matriz_aumentada(A, B, n)

    # Eliminação progressiva
    for k in range(n):
        # Inicializa o maior elemento para índice e pivô
        i_max = k
        pivo_max = A[i_max][k]

        # Encontra o maior elemento da coluna e associa ao pivô
        for i in range(k + 1, n):
            if (abs(A[i][k]) > pivo_max):
                pivo_max = A[i][k]
                i_max = i

        if not A[k][i_max]:  # Se o pivô for zero
            sys.exit('Divisao por zero detectada!')

        # Troca a linha do maior elemento com a linha atual
        if (i_max != k):
            troca_linha(A, k, i_max, n)

        for i in range(k + 1, n):
            # fatora f para zerar os elementos abaixo do pivô
            f = np.float128(A[i][k] / A[k][k])

            # subtrai a linha do pivô multiplicada pelo fator f
            for j in range(k + 1, n + 1):
                A[i][j] = np.float128(A[i][j] - f * A[k][j])

            # preenche a matriz triangular inferior com zeros
            A[i][k] = 0

    # Substituição regressiva
    for i in range(n - 1, -1, -1):
        # Inicializa o vetor solução com a última coluna da matriz
        solucao[i] = A[i][n]

        # Subtrai os elementos da solução já encontrados
        # j é i + 1 pois os elementos da diagonal principal já foram zerados
        for j in range(i + 1, n):
            solucao[i] = np.float128(solucao[i] - A[i][j] * solucao[j])

        # Divide o elemento da solução pelo elemento da diagonal principal
        solucao[i] = np.float128(solucao[i] / A[i][i])

    return solucao


# Método de Cholesky
def cholesky(A, B, n):

    L = [[0.0] * len(A) for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(A[i][i] - s) if (i == j) else (1.0 / L[j][j] *
                                                               (A[i][j] - s))

    #transposta da matriz L
    Lt = transposta(L, n)

    #vetor X de resposta
    X = substituicao_regressiva(Lt, substituicao_regressiva(L, B, n), n)

    return X


#Método de Gauss Seidel
def gauss_seidel(M, B, u, E, max_iteracoes):

    ordem = len(M[0])
    X = list(u)
    Xerro = list(u)  #vetor pra calcular o erro
    passos = 0

    print("ordem da matriz", ordem)

    for k in range(max_iteracoes):

        #percorre a matriz
        for i in range(ordem):
            #comeca a soma pelo termo do vetor fonte
            soma = B[i]
            div = 0
            for j in range(ordem):
                #separa o divisor
                if (i == j):
                    div = M[i][j]
                else:
                    soma += M[i][j] * X[j] * -1.0
            #cria vetor de solucoes para proxima iteracao com resultados da linha
            X[i] = soma / div

        if (k % 100 == 0):
            print("Gauss-Seidel fez " + str(k) + " iteracoes...")

        #se atingir o criterio de parada, interrompe e retorna os resultados
        erro = calculaErro(X, Xerro)

        if (erro < E):
            print("Terminou Gauss Seidel com erro de: ", erro)
            return X
        Xerro = list(X)

    print("Terminou Gauss Seidel com erro de: ", erro)
    print(
        "Gauss Seidel nao convergiu ou precisa de mais iteracoes para convergir"
    )
    return X


#Método Jacobi
def jacobi(M, B, u, E, max_iteracoes):

    ordem = len(M[0])
    X = list(u)
    Xerro = list(u)  #vetor pra calcular o erro
    passos = 0

    print("ordem da matriz", ordem)

    for k in range(max_iteracoes):

        #percorre a matriz
        for i in range(ordem):
            #comeca a soma pelo termo do vetor fonte
            soma = B[i]
            div = 0
            for j in range(ordem):
                passos += 1
                #separa o divisor
                if (i == j):
                    div = M[i][j]
                else:
                    soma += M[i][j] * Xerro[j] * -1.0
            #cria vetor de solucoes para proxima iteracao com resultados da linha
            X[i] = soma / div

        if (k % 100 == 0):
            print("Jacobi fez " + str(k) + " iteracoes...")

        #se atingir o criterio de parada, interrompe e retorna os resultados
        erro = calculaErro(X, Xerro)

        if (erro < E):
            print("Terminou Jacobi com erro de: ", erro)
            return X

        Xerro = list(X)

    print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
    return X


def main():
    questao = int(input("Qual a questão?"))
    n = int(input('Digite o numero de equacoes: '))
    A = np.zeros((n, n), dtype=np.float128)
    B = np.zeros(n, dtype=np.float128)
    solucao = np.zeros(n, dtype=np.float128)
    start = 0
    end = 0
    if questao == 1:
        for i in range(n):
            for j in range(n):
                A[i][j] = np.float128(1 / (i + j + 1))
            B[i] = np.float128(1 / (i + n + 1))
    elif questao == 2:
        k = int(np.sqrt(n))
        B = 10 * np.ones(n)
        for i in range(0, n):
            A[i, i] = -4
        for i in range(0, n - 1):
            A[i + 1, i] = 1
            A[i, i + 1] = 1
            if ((i + 1) % k == 0):
                A[i + 1, i] = 0
                A[i, i + 1] = 0
        for i in range(k, n):
            A[i - k, i] = 1
            A[i, i - k] = 1

        A = -A * (k - 1)**2

    x = int(input('Digite o metodo desejado:\n1 - Gauss com pivoteamento \n2 - Gauss sem pivoteamento \n3 - Decomposicao LU \n4 - Cholesky\n5 - Gauss Seidel\n6 - Jacobi\n'))
    if x == 1:
        print("Gauss pivoteado")
        start = time.time()
        solucao = gauss_pivoteamento(A, B, n)
        end = time.time()
    elif x == 2:
        print("Gauss sem pivoteado")
        start = time.time()
        solucao = gauss(A, B, n)
        end = time.time()

    elif x == 3:
        print("Decomposicao LU")
        start = time.time()
        solucao = decomposicao_LU(A, B, n)
        end = time.time()
    elif x == 4:
        print("Cholesky")
        start = time.time()
        solucao = cholesky(A, B, n)
        end = time.time()
    elif x == 5:
        start = time.time()
        solucao = gauss_seidel(A, B, [0.0] * n, 0.01, 1000)
        end = time.time()
    elif x == 6:
        start = time.time()
        solucao = jacobi(A, B, [0.0] * n, 0.01, 1000)
        end = time.time()
    else:
        print("valor de método incorreto")
        sys.exit(1)

    #xs = np.arange(0, n, 1)
    #plt.plot(xs, solucao, 'b')
    #plt.show()

    print('\nErro: ' + str(maxNorma(A, B, solucao, n)))
    print('\nElapsed Time ' + str(end - start))


if __name__ == "__main__":
    main()
