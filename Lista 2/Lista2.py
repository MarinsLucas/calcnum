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

def substituicao_regressiva(M, B):
    if(len(M)==0):
        return 0
    
    ordem = len(M[0])
    temp = 0
    passos = 0
    
    #cria array de solucao
    sol = np.array([0] * ordem, np.float64)
    
    #checa se e superior
    isSuperior = (M[ordem-1][0] == 0)
    
    #percorre as linhas da matriz (ta funcionando)
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

def substituicao_progressiva(M, n):
    solucao = np.zeros(n, dtype=np.float16)

    for i in range(n):
        # Inicializa o vetor solução com a última coluna da matriz
        solucao[i] = M[i][n]

        # Subtrai os elementos da solução já encontrados
        # j é i + 1 pois os elementos da diagonal principal já foram zerados
        for j in range(i):
            solucao[i] = np.float16(solucao[i] - M[i][j] * solucao[j])

        # Divide o elemento da solução pelo elemento da diagonal principal
        solucao[i] = np.float16(solucao[i]/M[i][i])

    return solucao

# Fim de funções auxiliares
# -----------------------------------------------#

# Método de Gauss sem pivoteamento
def gauss(matriz, n, solucao):
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


def decomposicao_LU(A, n, B):
    n = len(A)

    for j in range(n):
        for i in range(n): 

            #Parte L da matriz
            if i>= j:
                s = 0
                for k in range(0, j):
                    s+= A[i][j] * A[k][j]
                A[i][j] -= s
            else:
                s = 0
                for k in range(0, i):
                    s += A[i][k] * A[k][j]
                A[i][j] -= s

                A[i][j] = (A[i][j] - s)/float(A[i][i])
    #FIM DA DECOMPOSIÇÃO EM LU 


    #INÍCIO DA RESOLUÇÃO DA LU
    assert len(B) == n
    
    pprint.pprint(A)
    
    #Ly=b
    y = np.zeros(n, dtype=float)
    for i in range(n):
        s = 0
        for j in range(0, i):
            s += A[i][j]*y[j]
        y[i] = (B[i] - s)/ float(A[i][i])

    #Ux=y
    x = np.zeros(n, dtype=float)
    for i in range(n-1, -1 , -1):
        s = 0
        for j in range(i+1, n):
            s+= A[i][j]*x[j]
        x[i] = y[i] - s

    return x

# Método de Cholesky
def cholesky(A, n, B):
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

        G[i][j] = (A[i][j] - soma) / G[j][j] #TODO: será que precisa fazer a transporta da G?

    #gerando a matriz L lower
    L = [[0.0] * len(G) for _ in range(len(G))]
    for i in range(len(G)):
        for j in range(i+1):
            passos += 1
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if(i==j):
                L[i][j] = math.sqrt(G[i][i] - s)
            else:
                L[i][j] = (1.0 / L[j][j] * (G[i][j] - s))
    
    Y = substituicao_regressiva(L, B)
    L = transposta(L, n)
    retro2 = substituicao_regressiva(L,Y)

    return retro2

def main():

    # Lê o número de equações
    n = int(input('Insira o numero de equacoes: '))

    # Cria e inicializa a matriz ampliada
    matrizAmpliada = np.zeros((n, n + 1), dtype=np.float16)

    # Cria e inicializa o vetor solução
    solucao = np.zeros(n, dtype=np.float16)

    # Insere os coeficientes na matriz ampliada
    for i in range(n):
        for j in range(n):
            matrizAmpliada[i][j] = np.float16(1/(i + j + 1))

    for i in range(n):
        matrizAmpliada[i][n] = np.float16(1/(i + n + 1))
        
        
    termosIndependentes = np.zeros(n, dtype=np.float16)
    matriz = np.zeros((n, n), dtype=np.float16)
    
    for i in range(n):
        termosIndependentes[i] = matrizAmpliada[i][n]
        
    for i in range(n):
        for j in range(n):
            matriz[i][j] = matrizAmpliada[i][j]
            
    # Imprime a matriz ampliada
    print('\nMatriz ampliada: ')
    for i in range(n):
        for j in range(n + 1):
            print('%0.2f' % matrizAmpliada[i][j], end='\t')
        print('')

    # Pergunta se o usuário quer usar o método de gauss com ou sem pivoteamento, e dependendo da resposta chama a função necessária
    metodo = int(input(
        '\n1 - Gauss sem pivoteamento\n2 - Gauss com pivoteamento\n3 - Decomposicao LU\n4 - Cholesky\n\nInsira o metodo que deseja usar: '))
    if metodo == 1:
        solucao = gauss(matrizAmpliada, n, solucao)
    elif metodo == 2:
        solucao = gauss_pivoteamento(matrizAmpliada, n, solucao)
    elif metodo == 3:
        solucao = decomposicao_LU(matriz, n, termosIndependentes)
    elif metodo == 4:
        solucao = cholesky(matriz, n, termosIndependentes)
    else:
        sys.exit('Metodo invalido!')

    # Imprime a solução
    print('\nSolucao: ')
    for i in range(n):
        print('X%d = %0.2f' % (i, solucao[i]), end='\t')


if __name__ == "__main__":
    main()
