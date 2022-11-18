import numpy as np
import pprint
import sys
import math

# -----------------------------------------------#
# Funções auxiliares
def transposta(A, n):
    B = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            B[i][j] = A[j][i]
    return B

def troca_linha(A, i, j, n):
    for k in range(n + 1):
        temp = A[i][k]
        A[i][k] = A[j][k]
        A[j][k] = temp

def matriz_pivo(M):
    m = len(M)

    # Cria a A identidade com valores float16
    A_id = np.identity(m, dtype=np.float16)

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
            troca_linha(A_id, k, i_max, m)

    return A_id


def matriz_multi(M, N):
    m = len(M)
    n = len(N[0])
    p = len(M[0])
    A = np.zeros((m, n), dtype=np.float16)

    for i in range(m):
        for j in range(n):
            for k in range(p):
                A[i][j] += M[i][k] * N[k][j]

    return A

def substituicao_regressiva(M, B, n):
    if(len(M)==0):
        return 0
    
    ordem = n
    temp = 0
    passos = 0
    
    #cria array de solucao
    sol = np.array([0] * ordem, np.float64)
    
    #checa se e superior
    isSuperior = (M[ordem-1][0] == 0)
    
    #percorre as linhas da A 
    for n in range(0, ordem):  
        
        #define o i se for superior ou inferior
        if(isSuperior):
            i = ordem - n - 1
        else:
            i = n
        
        #valor do vetor solucao da linha correspondente
        temp = B[i]
        
        #soma os valores exceto o valor do pivo
        for j in range(0, ordem):
            if(j != i):
                temp += sol[j] * M[i][j] * -1
                passos += 1
                
            
        #calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i]
        

        
    return list(sol)

def substituicao_progressiva(M, B, n):
    if(len(M)==0):
        return 0
    
    ordem = len(M[0])
    temp = 0
    passos = 0
    
    #cria array de solucao
    sol = np.array([0] * ordem, np.float64)
    
    #checa se e superior
    isSuperior = (M[ordem-1][0] == 0)
    
    #percorre as linhas da A 
    for n in range(0, ordem):  
        
        #define o i se for superior ou inferior
        if(isSuperior):
            i = n
        else:
            i = ordem - n - 1
        
        #valor do vetor solucao da linha correspondente
        temp = B[i]
        
        #soma os valores exceto o valor do pivo
        for j in range(0, ordem):
            if(j != i):
                temp += sol[j] * M[i][j] * -1
                passos += 1
                
            
        #calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i]
        

        
    return sol

def matriz_aumentada(A, B, n):
    # Cria uma nova matriz com n linhas e n + 1 colunas e adiciona a coluna B à matriz A
    M = np.zeros((n, n + 1), dtype=np.float16)

    for i in range(n):
        for j in range(n):
            M[i][j] = A[i][j]

    for i in range(n):
        M[i][n] = B[i]

    return M

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

# Fim de funções auxiliares
# -----------------------------------------------#

# Método de Gauss sem pivoteamento
def gauss(A, B, n):
    solucao = np.zeros(n, dtype=np.float16)
    
    A = matriz_aumentada(A, B, n)
    
    # Substituição progressiva
    for i in range(n):
        if (A[i][i] == 0.0):
            sys.exit('Divisao por zero detectada!')

        for j in range(i + 1, n):
            razao = np.float16(A[j][i]/A[i][i])

            for k in range(n + 1):
                A[j][k] = np.float16(A[j][k] - razao * A[i][k])

    # Substituição regressiva
    solucao[n - 1] = np.float16(A[n - 1][n]/A[n - 1][n - 1])

    for i in range(n - 2, -1, -1):
        solucao[i] = np.float16(A[i][n])

        for j in range(i + 1, n):
            solucao[i] = np.float16(solucao[i] - A[i][j] * solucao[j])

        solucao[i] = np.float16(solucao[i]/A[i][i])

    return solucao


# Método de Gauss com pivoteamento
def gauss_pivoteamento(A, B, n):
    
    solucao = np.zeros(n, dtype=np.float16)
    
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
            f = np.float16(A[i][k]/A[k][k])

            # subtrai a linha do pivô multiplicada pelo fator f
            for j in range(k + 1, n + 1):
                A[i][j] = np.float16(A[i][j] - f * A[k][j])

            # preenche a matriz triangular inferior com zeros
            A[i][k] = 0

    # Substituição regressiva
    for i in range(n - 1, -1, -1):
        # Inicializa o vetor solução com a última coluna da matriz
        solucao[i] = A[i][n]

        # Subtrai os elementos da solução já encontrados
        # j é i + 1 pois os elementos da diagonal principal já foram zerados
        for j in range(i + 1, n):
            solucao[i] = np.float16(solucao[i] - A[i][j] * solucao[j])

        # Divide o elemento da solução pelo elemento da diagonal principal
        solucao[i] = np.float16(solucao[i]/A[i][i])

    return solucao

# Método da Decomposição LU


def decomposicao_LU(A, B, n):
    A = decompoeLU(A, n)
    return resolveLU(A, B, n)

# Método de Cholesky
def cholesky(A, B, n):
    G = np.zeros((n,n))

    passos = 0

    for i in range(n):
        soma = 0
        #gerando diagonal principal
        for j in range(i+1):
            passos+=1
            if(j != i):
                #soma dos quadrados
                soma += G[i][j]**2
            else:
                G[i][j] = math.sqrt(A[i][i] - soma)
        
    #gerando demais elementos
    for i in range(j+1):
        soma = 0
        for k in range(j-1):
            passos+=1
            soma += G[i][k]* G[j][k]

        G[i][j] = (A[i][j] - soma) / G[j][j]

    #gerando a A L lower
    L = [[0.0] * len(G) for _ in range(len(G))]
    for i in range(len(G)):
        for j in range(i+1):
            passos += 1
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if(i==j):
                L[i][j] = math.sqrt(G[i][i] - s)
            else:
                L[i][j] = (1.0 / L[j][j] * (G[i][j] - s))

    Y = substituicao_regressiva(L, B, n)
    L = transposta(L, n)
    sol = substituicao_regressiva(L, Y, n)

    return sol

def main():
    n = 3
    A = np.array([[2, 4, -2], [4, 9, -3], [-2, -3, 7]], dtype=np.float16)
    B = np.array([2, 8, 10], dtype=np.float16)
    solucao = np.zeros(3)
    
    """ print('A original:')
    pprint.pprint(A) """
    
    print('Solução do sistema pelo método de Gauss com pivoteamento:')
    solucao = gauss_pivoteamento(A, B, n)
    pprint.pprint(solucao)


if __name__ == "__main__":
    main()
