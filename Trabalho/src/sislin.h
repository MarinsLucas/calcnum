#include <iostream>
#include <vector>

using namespace std;

typedef long double flu;

typedef vector<vector<flu>> Matrix;
typedef vector<flu> Vector;

//FIXME: Acho que dá pra refazer essa função com ponteiros!!
void troca_linha(flu** A, int i, int j, int n)
{
    double temp_i, temp_j;
    for (int k = 0; k < n + 1; k++)
    {
        temp_i = A[i][k];
        temp_j = A[j][k];
        A[i][k] = temp_j;
        A[j][k] = temp_i;
    }
}

flu** matriz_aumentada(flu** A, flu* B, int n)
{
    // Cria uma nova matriz com n linhas e n + 1 colunas e adiciona a coluna B à matriz A
    flu** M = new flu*[n];

    for (int i = 0; i < n; i++)
    { 
        M[i] = new flu[n+1];
        for (int j = 0; j < n; j++)
        {
            M[i][j] = A[i][j];
        }
        M[i][n] = B[i];
    }

    return M;
}

Vector gauss_pivoteamento(flu** A, flu* B, int n)
{
    vector<flu> solucao(n, 0.0);

    A = matriz_aumentada(A, B, n);

    // Eliminação progressiva
    for (int k = 0; k < n; k++)
    {
        // Inicializa o maior elemento para índice e pivô
        int i_max = k;
        double pivo_max = A[i_max][k];

        // Encontra o maior elemento da coluna e associa ao pivô
        for (int i = k + 1; i < n; i++)
        {
            if (abs(A[i][k]) > pivo_max)
            {
                pivo_max = A[i][k];
                i_max = i;
            }
        }

        if (A[k][i_max] == 0)
        { // Se o pivô for zero
            std::cout << "Divisao por zero detectada!" << endl;
            exit(1);
        }

        // Troca a linha do maior elemento com a linha atual
        if (i_max != k)
        {
            troca_linha(A, k, i_max, n);
        }

        for (int i = k + 1; i < n; i++)
        {
            // fatora f para zerar os elementos abaixo do pivô
            double f = (double)A[i][k] / (double)A[k][k]; // TODO: verificar se não é necessário converter para float128

            // subtrai a linha do pivô multiplicada pelo fator f
            for (int j = k + 1; j < n + 1; j++)
            {
                A[i][j] = A[i][j] - f * A[k][j];
            }

            // preenche a matriz triangular inferior com zeros
            A[i][k] = 0;
        }
    }

    // Substituição regressiva
    for (int i = n - 1; i >= 0; i--)
    {
        // Inicializa o vetor solução com a última coluna da matriz
        solucao[i] = A[i][n];

        // Subtrai os elementos da solução já encontrados
        // j é i + 1 pois os elementos da diagonal principal já foram zerados
        for (int j = i + 1; j < n; j++)
        {
            solucao[i] = solucao[i] - A[i][j] * solucao[j];
        }

        // Divide o elemento da solução pelo elemento da diagonal principal
        solucao[i] = (double)solucao[i] / (double)A[i][i];
    }

    return solucao;
}

// Método de Gauss sem pivoteamento
Vector gauss(flu** A, flu* B, int n)
{
    Vector solucao(n, 0.0);

    A = matriz_aumentada(A, B, n);

    // Substituição progressiva
    for (int i = 0; i < n; i++)
    {
        if (A[i][i] == 0.0)
        {
            std::cout << "Divisao por zero detectada!" << endl;
            exit(1);
        }

        for (int j = i + 1; j < n; j++)
        {
            double razao = ((double)A[j][i] / (double)A[i][i]);

            for (int k = 0; k < n + 1; k++)
            {
                A[j][k] = (double)(A[j][k] - razao * A[i][k]);
            }
        }
    }

    // Substituição regressiva
    solucao[n - 1] = (double)(A[n - 1][n] / (double)A[n - 1][n - 1]);

    for (int i = n - 2; i >= 0; i--)
    {
        solucao[i] = A[i][n];

        for (int j = i + 1; j < n; j++)
        {
            solucao[i] = (double)(solucao[i] - A[i][j] * solucao[j]);
        }

        solucao[i] = ((double)solucao[i] / (double)A[i][i]);
    }

    return solucao;
}

void printVector(flu* a, int n)
{
    for(int i =0; i<n; i++)
    {
        cout<<a[i]<<endl;
    }
}

void printMatrix(flu**a, int n)
{
    for(int i=0; i<n; i++)
    {
        for(int j=0 ; j<n; j++)
        {
            cout<<a[j][i]<<";";
        }
        cout<<endl;
    }
}