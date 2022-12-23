#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <chrono>
#include "pbPlots.hpp"
#include "supportLib.hpp"

#define eps 1e-9

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;
typedef chrono::high_resolution_clock Clock;

// TODO: - Calculo de erro
// TODO: (1) comparar o tempo de execução de todos para cada n em uma tabela
// TODO: (1) tabela de erros
// TODO: (1) conclusões sobre o determinante da matriz
// TODO: (2) tabela comparando os tempos de execução dos métodos da questão 1 e da 2 adotando direfentes valores de epslon e usando
//           n = 81, 289, 1089, 4225, 16641
//  -----------------------------------------------#

// Funções auxiliares

void imprimeVetor(Vector A){
    for(int i = 0; i < A.size(); i++){
        cout << A[i] << " ";
    }
    cout << endl;
}

void DrawGraphic(Vector solucao, int n)
{
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
    Vector xs; 
    for(int i =0; i<n; i++)
    {
        xs.push_back(i);
    }
    bool success = DrawScatterPlot(imageReference, 1080, 720, &xs,&solucao, errorMessage);
    if(success){
		vector<double> *pngdata = ConvertToPNG(imageReference->image);
		WriteToFile(pngdata, "example1.png");
		DeleteImage(imageReference->image);
	}else{
		cerr << "Error: ";
		for(wchar_t c : *errorMessage->string){
			wcerr << c;
		}
		cerr << endl;
	}

	FreeAllocations();
}

void imprimeMatriz(Matrix A)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            std::cout << A[i][j] << " ";
        }
        std::cout << endl;
    }
}

Matrix transposta(Matrix A, int n)
{
    Matrix B(n, Vector(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            B[i][j] = A[j][i];
        }
    }

    return B;
}

void troca_linha(Matrix A, int i, int j, int n)
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

double normaMaximo(Vector x)
{
    double maximo = abs(x[0]);

    for (unsigned int i = 0; i < x.size(); i++)
    {
        if (abs(x[i]) > maximo)
        {
            maximo = abs(x[i]);
        }
    }

    return maximo;
}

double distanciaMaximo(Vector x1, Vector x2)
{
    if (x1.size() != x2.size())
    {
        std::cout << "O tamanho dos vetores x1 e x2 precisa ser o mesmo" << endl;
        return 0;
    }

    double dist = abs(x1[0] - x2[0]);

    for (unsigned int i = 0; i < x1.size(); i++)
    {
        if (abs(x1[i] - x2[i]) > dist)
        {
            dist = abs(x1[i] - x2[i]);
        }
    }

    return dist;
}

double calculaErro(Vector x_prox, Vector x_atual)
{
    return distanciaMaximo(x_prox, x_atual) / normaMaximo(x_prox);
}

Vector substituicao_regressiva(Matrix M, Vector B, int tam)
{
    if (tam == 0)
    {
        return {};
    }

    int ordem = tam;
    double temp = 0;

    // cria array de solucao
    Vector sol(ordem);

    // checa se e superior
    bool isSuperior = (M[ordem - 1][0] == 0);

    // percorre as linhas da A
    for (int n = 0; n < ordem; n++)
    {

        // define o i se for superior ou inferior
        int i = (isSuperior ? ordem - n - 1 : n);

        // valor do vetor solucao da linha correspondente
        temp = B[i];

        // soma os valores exceto o valor do pivo
        for (int j = 0; j < ordem; j++)
        {
            if (j != i)
            {
                temp += sol[j] * M[i][j] * -1;
            }
        }

        // calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = (double)temp / (double)M[i][i];
    }

    return sol;
}

Vector substituicao_progressiva(Matrix M, Vector B, int tam)
{
    if (tam == 0)
    {
        return {};
    }

    int ordem = tam;
    double temp = 0;

    // cria array de solucao
    Vector sol(ordem);

    // checa se e superior
    bool isSuperior = (M[ordem - 1][0] == 0);

    // percorre as linhas da A
    for (int n = 0; n < ordem; n++)
    {

        // define o i se for superior ou inferior
        int i = (isSuperior ? ordem - n - 1 : n);

        // valor do vetor solucao da linha correspondente
        temp = B[i];

        // soma os valores exceto o valor do pivo
        for (int j = 0; j < ordem; j++)
        {
            if (j != i)
            {
                temp += sol[j] * M[i][j] * -1;
            }
        }

        // calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = (double)temp / (double)M[i][i];
    }

    return sol;
}

Matrix matriz_aumentada(Matrix A, Vector B, int n)
{
    // Cria uma nova matriz com n linhas e n + 1 colunas e adiciona a coluna B à matriz A
    Matrix M(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            M[i][j] = A[i][j];
        }
        M[i][n] = B[i];
    }

    return M;
}

vector<Matrix> decompoeLU(Matrix M, int n)
{

    int ordem = n;

    Matrix L = Matrix(ordem, vector<double>(ordem, 0));
    Matrix U = Matrix(ordem, vector<double>(ordem, 0));

    // percorre as linhas da A
    for (int j = 0; j < ordem; j++)
    {
        U[0][j] = M[0][j];
    }

    for (int i = 0; i < ordem; i++)
    {
        L[i][0] = (double)M[i][0] / (double)U[0][0];
    }

    for (int i = 0; i < ordem; i++)
    {
        // Calcula L
        for (int j = i + 1; j < ordem; j++)
        {
            double soma = 0.0;
            for (int k = j; k < ordem; k++)
            {
                soma += L[i][k] * U[k][j];
            }
            L[i][j] = M[i][j] - soma;
        }

        // Calcula U
        for (int j = i; j < ordem; j++)
        {
            double soma = 0.0;
            for (int k = i; k < ordem; k++)
            {
                soma += L[i][k] * U[k][j];
            }
            U[i][j] = (double)(M[i][j] - soma) / (double)L[i][i];
        }
    }

    vector<Matrix> LU = {L, U};

    return LU;
}

vector<double> decomposicao_LU(Matrix M, Vector B, int n)
{
    if (n == 0)
    {
        return {};
    }

    int ordem = n;

    vector<double> y(ordem);
    vector<double> x(ordem);
    Matrix L = Matrix(ordem, vector<double>(ordem, 0));
    Matrix U = Matrix(ordem, vector<double>(ordem, 0)); 

    // percorre as linhas da A
    for (int j = 0; j < ordem; j++)
    {
        U[0][j] = M[0][j];
    }

    for (int i = 0; i < ordem; i++)
    {
        L[i][0] = (double)M[i][0] / (double)U[0][0];
    }

    for (int i = 0; i < ordem; i++)
    {
        // Calcula L
        for (int j = 0; j < i+1; j++)
        {
            double soma = 0.0;
            for (int k = 0; k < j; k++)
            {
                soma += L[i][k] * U[k][j];
            }
            L[i][j] = M[i][j] - soma;
        }

        // Calcula U
        for (int j = i; j < ordem; j++)
        {
            double soma = 0.0;
            for (int k = 0; k < i; k++)
            {
                soma += L[i][k] * U[k][j];
            }
            U[i][j] = (double)(M[i][j] - soma) / (double)L[i][i];
        }
        //imprimeMatriz(L);
    }

    // Resolução L * y = b

    y[0] = (double)B[0] / (double)L[0][0];
    for (int i = 1; i < ordem; i++)
    {
        double soma = 0.0;
        for (int j = 0; j < i; j++)
        {
            soma += L[i][j] * y[j];
        }
        y[i] = (double)(B[i] - soma) / (double)L[i][i];
    }

    // Resolucao U * x = y
    x[ordem - 1] = (double)y[ordem - 1] / (double)U[ordem - 1][ordem - 1];
    for (int i = ordem - 1; i > -1; i--)
    {
        double soma = y[i];
        for (int j = i + 1; j < ordem; j++)
        {
            soma = soma - U[i][j] * x[j];
        }
        x[i] = (double)soma / (double)U[i][i];
    }
    return x;
}

Vector erro(Matrix A, double *B, double *x, int n)
{
    Vector X(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            X[i] = abs(A[i][j] * x[i] - B[i]);
        }
    }
    return X;
}

double maxNorma(Matrix A, Vector B, Vector x, int n)
{
    double max = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double temp = abs(A[i][j] * x[i] - B[i]);
            if (temp > max)
            {
                max = temp;
            }
        }
    }

    return max;
}

// Fim de funções auxiliares

// Método de Gauss sem pivoteamento
Vector gauss(Matrix A, Vector B, int n)
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

// Método de Gauss com pivoteamento
Vector gauss_pivoteamento(Matrix A, Vector B, int n)
{
    Vector solucao(n, 0.0);

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

// Método de Cholesky
Vector cholesky(Matrix A, Vector B, int n)
{
    Matrix L(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double s = 0;
            for (int k = 0; k < j; k++)
            {
                s += L[i][k] * L[j][k];
            }
            L[i][j] = (i == j) ? sqrt(A[i][i] - s) : (1.0 / (double)L[j][j] * (double)(A[i][j] - s));
        }
    }

    // transposta da matriz L
    Matrix Lt = transposta(L, n);

    // vetor X de resposta
    Vector X = substituicao_regressiva(Lt, substituicao_regressiva(L, B, n), n);

    return X;
}

// Método de Gauss Seidel
Vector gauss_seidel(Matrix M, Vector B, Vector u, double E, int max_iteracoes)
{
    int ordem = M[0].size();
    double erro = 0.0;
    Vector X = u;
    Vector Xerro = u; // vetor pra calcular o erro

    std::cout << "ordem da matriz " << ordem << endl;

    for (int k = 0; k < max_iteracoes; k++)
    {
        // percorre a matriz
        for (int i = 0; i < ordem; i++)
        {
            // comeca a soma pelo termo do vetor fonte
            double soma = B[i];
            double div = 0;
            for (int j = 0; j < ordem; j++)
            {
                // separa o divisor
                if (i == j)
                {
                    div = M[i][j];
                }
                else
                {
                    soma += M[i][j] * X[j] * -1.0;
                }
            }
            // cria vetor de solucoes para proxima iteracao com resultados da linha
            X[i] = (double)soma / (double)div;
        }

        if (k % 100 == 0)
        {
            std::cout << "Gauss-Seidel fez " << k << " iteracoes..." << endl;
        }

        // se atingir o criterio de parada, interrompe e retorna os resultados
        erro = calculaErro(X, Xerro);

        if (erro < E)
        {
            std::cout << "Terminou Gauss Seidel com erro de: " << erro << endl;
            return X;
        }
        Xerro = X;
    }

    std::cout << "Terminou Gauss Seidel com erro de : " << erro << endl;
    std::cout << "Gauss Seidel nao convergiu ou precisa de mais iteracoes para convergir" << endl;
    return X;
}

// Método Jacobi
Vector jacobi(Matrix M, Vector B, Vector u, double E, int max_iteracoes)
{
    int ordem = M[0].size();
    Vector X = u;
    Vector Xerro = u; // vetor pra calcular o erro

    std::cout << "ordem da matriz " << ordem << endl;

    for (int k = 0; k < max_iteracoes; k++)
    {
        // percorre a matriz
        for (int i = 0; i < ordem; i++)
        {
            // comeca a soma pelo termo do vetor fonte
            double soma = B[i];
            double div = 0;
            for (int j = 0; j < ordem; j++)
            {
                // separa o divisor
                if (i == j)
                {
                    div = M[i][j];
                }
                else
                {
                    soma += M[i][j] * Xerro[j] * -1.0;
                }
            }
            // cria vetor de solucoes para proxima iteracao com resultados da linha
            X[i] = (double)soma / (double)div;
        }

        if (k % 100 == 0)
        {
            std::cout << "Jacobi fez " << k << " iteracoes..." << endl;
        }

        // se atingir o criterio de parada, interrompe e retorna os resultados
        double erro = calculaErro(X, Xerro);

        if (erro < E)
        {
            std::cout << "Terminou Jacobi com erro de: " << erro << endl;
            return X;
        }

        Xerro = X;
    }

    std::cout << "Jacobi nao convergiu ou precisa de mais iteracoes para convergir" << endl;
    return X;
}

int main()
{
    int questao;
    int n;

    std::cout << "Digite o numero da questao: " << endl;
    std::cin >> questao;

    std::cout << "Digite o numero de equacoes: " << endl;
    std::cin >> n;

    Matrix A(n, Vector(n, 0.0));
    Vector B(n, 0.0);
    Vector solucao(n, 0.0);
    std::chrono::_V2::system_clock::time_point start, end;

        if (questao == 1)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i][j] = ((double)1 / (double)(i + j + 1));
                }
                B[i] = ((double)1 / (double)(i + n + 1));
            }
        }
        else if (questao == 2)
        {
            int k = (int)sqrt(n);

            for (int i = 0; i < n; i++)
            {
                B[i] = 10;
            }

            for (int i = 0; i < n; i++)
            {
                A[i][i] = -4;
            }
            for (int i = 0; i < n - 1; i++)
            {
                A[i + 1][i] = 1;
                A[i][i + 1] = 1;
                if ((i + 1) % k == 0)
                {
                    A[i + 1][i] = 0;
                    A[i][i + 1] = 0;
                }
            }
            for (int i = k; i < n; i++)
            {
                A[i - k][i] = 1;
                A[i][i - k] = 1;
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i][j] = ((double)A[i][j] * (double)pow(k - 1, 2));
                }
            }
        }

        /* std::cout << "Matriz A: " << endl;
        imprimeMatriz(A); */

         std::cout << "Digite o metodo desejado:" << endl;
        std::cout << "1 - Gauss com pivoteamento" << endl;
        std::cout << "2 - Gauss sem pivoteamento" << endl;
        std::cout << "3 - Decomposicao LU" << endl;
        std::cout << "4 - Cholesky" << endl;
        std::cout << "5 - Gauss Seidel" << endl;
        std::cout << "6 - Jacobi" << endl; 

        int metodo;
        cin>>metodo;

        if (metodo == 1)
        {
            std::cout << "Gauss pivoteado" << endl;
            start = Clock::now();
            solucao = gauss_pivoteamento(A, B, n);
            end = Clock::now();
        }
        else if (metodo == 2)
        {
            std::cout << "Gauss sem pivoteado" << endl;
            start = Clock::now();
            solucao = gauss(A, B, n);
            end = Clock::now();
        }
        else if (metodo == 3)
        {
            std::cout << "Decomposicao LU" << endl;
            start = Clock::now();
            solucao = decomposicao_LU(A, B, n);
            end = Clock::now();
        }
        else if (metodo == 4)
        {
            std::cout << "Cholesky" << endl;
            start = Clock::now();
            solucao = cholesky(A, B, n);
            end = Clock::now();
        }
        else if (metodo == 5)
        {
            start = Clock::now();
            std::cout<<"Gauss Seidel"<<endl;
            solucao = gauss_seidel(A, B, solucao, 0.1, 1000);
            end = Clock::now();
        }
        else if (metodo == 6)
        {
            start = Clock::now();
            std::cout<<"Jacobi"<<endl;
            solucao = jacobi(A, B, solucao, 0.1, 1000);
            end = Clock::now();
        }
        else
        {
            std::cout << "valor de método incorreto" << endl;
            exit(1);
        }
    
    auto diff = chrono::duration_cast<chrono::microseconds>(end - start).count();
    double seconds = diff / 1e6;
    DrawGraphic(solucao, n);   
    
    std::cout << "\nErro: " << maxNorma(A, B, solucao, n) << endl;
    std::cout << "\nTempo percorrido: " << seconds << " segundos" << endl;
    

    

    
    return 0;
}


