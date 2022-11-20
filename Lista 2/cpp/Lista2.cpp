#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "pbPlots.hpp"
#include "supportLib.hpp"
#include <sstream>
#include <string>

//TODO: Consertar o Decomposição LU
//TODO: Confirmar se os métodos de gauss estão funcionando


using namespace std;

// -----------------------------------------------//
// Funções auxiliares

vector<vector<double>> transposta(vector<vector<double>> A, int n){
    vector<vector<double>> B(n, vector<double>(n));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            B[i][j] = A[j][i];
        }
    }
    return B;
}

void troca_linha(vector<vector<double>> &A, int i, int j, int n){
    for(int k = 0; k < n + 1; k++){
        double temp = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = temp;
    }
}

vector<double> substituicao_regressiva(vector<vector<double>> M, vector<double> B, int n){
    if(M.size() == 0){
        return {};
    }
    
    int ordem = n;
    double temp = 0;
    
    //cria array de solucao
    vector<double> sol(ordem);
    
    //checa se e superior
    bool isSuperior = (M[ordem-1][0] == 0);
    
    //percorre as linhas da A 
    for(int n = 0; n < ordem; n++){  
        
        //define o i se for superior ou inferior
        int i;
        if(isSuperior){
            i = ordem - n - 1;
        }else{
            i = n;
        }
        
        //valor do vetor solucao da linha correspondente
        temp = B[i];
        
        //soma os valores exceto o valor do pivo
        for(int j = 0; j < ordem; j++){
            if(j != i){
                temp += sol[j] * M[i][j] * -1;
            }
        }
            
        //calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i];
    }
        
    return sol;
}

void plotaGrafico(vector<double> &solucao, vector<double> &xs)
{
    bool success;
	StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
	RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();

    success = DrawScatterPlot(imageReference, 1920, 1080, &xs, &solucao, errorMessage);

    if(success){
		vector<double> *pngdata = ConvertToPNG(imageReference->image);
        string file_name;
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

vector<double> substituicao_progressiva(vector<vector<double>> M, vector<double> B, int n){
    if(M.size() == 0){
        return {};
    }
    
    int ordem = M[0].size();
    double temp = 0;
    
    //cria array de solucao
    vector<double> sol(ordem);
    
    //checa se e superior
    bool isSuperior = (M[ordem-1][0] == 0);
    
    //percorre as linhas da A 
    for(int n = 0; n < ordem; n++){  
        
        //define o i se for superior ou inferior
        int i;
        if(isSuperior){
            i = n;
        }else{
            i = ordem - n - 1;
        }
        
        //valor do vetor solucao da linha correspondente
        temp = B[i];
        
        //soma os valores exceto o valor do pivo
        for(int j = 0; j < ordem; j++){
            if(j != i){
                temp += sol[j] * M[i][j] * -1;
            }
        }
            
        //calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i];
    }
        
    return sol;
}

vector<vector<double>> matriz_aumentada(vector<vector<double>> A, vector<double> B, int n){
    // Cria uma nova matriz com n linhas e n + 1 colunas e adiciona a coluna B à matriz A
    vector<vector<double>> M(n, vector<double>(n + 1));

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            M[i][j] = A[i][j];
        }
    }

    for(int i = 0; i < n; i++){
        M[i][n] = B[i];
    }

    return M;
}

vector<vector<double>> decompoeLU(vector<vector<double>> A, int n){
    for(int j = 0; j < n; j++){
        for(int i = 0; i < n; i++){
            // Parte L da matriz
            if(i >= j){
                double s = 0;
                for(int k = 0; k < j; k++){
                    s += A[i][k] * A[k][j];
                }
                A[i][j] -= s;
                
            // Parte U da matriz
            }else{ 
                double s = 0;
                for(int k = 0; k < i; k++){
                    s += A[i][k] * A[k][j];
                }
                    
                A[i][j] = (A[i][j] - s) / A[i][i];
            }
        }
    }
    return A;
}

vector<double> resolveLU(vector<vector<double>> LU, vector<double> B, int n){
    
    // Resolve Ly = b
    vector<double> y(n);
    
    for(int i = 0; i < n; i++){
        double s = 0;
        for(int j = 0; j < i; j++){
            s += LU[i][j] * y[j];
        }
        y[i] = (B[i] - s) / LU[i][i];
    }
        
    // Resolve Ux = y
    vector<double> x(n);
    
    for(int i = n-1; i >= 0; i--){
        double s = 0;
        for(int j = i+1; j < n; j++){
            s += LU[i][j] * x[j];
        }
        x[i] = (y[i] - s);
    }
    
    return x;
}

double maxNorma(vector<vector<double>> A, vector<double> B, vector<double> x, int n){
    double max = abs(A[0][0] * x[0] - B[0]);
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            double temp = abs(A[i][j] * x[i] - B[i]);
            if(temp > max){
                max = temp;
            }
        }
    }
    return max;
}

// Fim de funções auxiliares
// -----------------------------------------------//

// Método de Gauss sem pivoteamento
vector<double> gauss(vector<vector<double>> A, vector<double> B, int n){
    vector<double> solucao(n);
    
    A = matriz_aumentada(A, B, n);
    
    // Substituição progressiva
    for(int i = 0; i < n; i++){
        if (A[i][i] == 0.0){
            cout << "Divisao por zero detectada!" << endl;
            return {};
        }

        for(int j = i + 1; j < n; j++){
            double razao = A[j][i]/A[i][i];

            for(int k = 0; k < n + 1; k++){
                A[j][k] = A[j][k] - razao * A[i][k];
            }
        }
    }

    // Substituição regressiva
    solucao[n - 1] = A[n - 1][n]/A[n - 1][n - 1];

    for(int i = n - 2; i >= 0; i--){
        solucao[i] = A[i][n];

        for(int j = i + 1; j < n; j++){
            solucao[i] = solucao[i] - A[i][j] * solucao[j];
        }

        solucao[i] = solucao[i]/A[i][i];
    }

    return solucao;
}

// Método de Gauss com pivoteamento
vector<double> gauss_pivoteamento(vector<vector<double>> A, vector<double> B, int n){
    
    vector<double> solucao(n);
    
    A = matriz_aumentada(A, B, n);
    
    // Eliminação progressiva
    for(int k = 0; k < n; k++){
        // Inicializa o maior elemento para índice e pivô
        int i_max = k;
        double pivo_max = A[i_max][k];

        // Encontra o maior elemento da coluna e associa ao pivô
        for(int i = k + 1; i < n; i++){
            if (abs(A[i][k]) > pivo_max){
                pivo_max = A[i][k];
                i_max = i;
            }
        }

        if(!A[k][i_max]){  // Se o pivô for zero
            cout << "Divisao por zero detectada!" << endl;
            return {};
        }

        // Troca a linha do maior elemento com a linha atual
        if (i_max != k){
            troca_linha(A, k, i_max, n);
        }

        for(int i = k + 1; i < n; i++){
            // fatora f para zerar os elementos abaixo do pivô
            double f = A[i][k]/A[k][k];

            // subtrai a linha do pivô multiplicada pelo fator f
            for(int j = k + 1; j < n + 1; j++){
                A[i][j] = A[i][j] - f * A[k][j];
            }

            // preenche a matriz triangular inferior com zeros
            A[i][k] = 0;
        }
    }

    // Substituição regressiva
    for(int i = n - 1; i >= 0; i--){
        // Inicializa o vetor solução com a última coluna da matriz
        solucao[i] = A[i][n];

        // Subtrai os elementos da solução já encontrados
        // j é i + 1 pois os elementos da diagonal principal já foram zerados
        for(int j = i + 1; j < n; j++){
            solucao[i] = solucao[i] - A[i][j] * solucao[j];
        }

        // Divide o elemento da solução pelo elemento da diagonal principal
        solucao[i] = solucao[i]/A[i][i];
    }

    return solucao;
}

// Método da Decomposição LU
vector<double> decomposicao_LU(vector<vector<double>> A, vector<double> B, int n){
    A = decompoeLU(A, n);
    return resolveLU(A, B, n);
}

// Método de Cholesky
vector<double> cholesky(vector<vector<double>> A, vector<double> B, int n){    
    int passos = 0;
    
    vector<vector<double>> L(A.size(), vector<double>(A.size()));
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < i+1; j++){
            passos += 1;
            double s = 0;
            for(int k = 0; k < j; k++){
                s += L[i][k] * L[j][k];
            }
            L[i][j] = (i == j) ? sqrt(A[i][i] - s) : (1.0 / L[j][j] * (A[i][j] - s));
        }
    }

    //transposta da matriz L
    vector<vector<double>> Lt = transposta(L, n);
    
    //vetor X de resposta
    vector<double> X = substituicao_regressiva(Lt, substituicao_regressiva(L, B, n), n);
    
    return X;
}

//Função para calcular a distância máxima entre dois vetores
double distanciaMaximo(vector<double> x1, vector<double> x2){
    if(x1.size() != x2.size()){
        cout << "O tamanho dos vetores x1 e x2 precisa ser o mesmo" << endl;
        return 0;
    }
    
    int size = x1.size();
    double dist = abs(x1[0] - x2[0]);  
    
    for(int i = 0; i < size; i++){
        double temp = abs(x1[i] - x2[i]);
        if(temp > dist){
            dist = temp;
        }
    }
    
    return dist;
}

//Função para calcular o erro
double calculaErro(vector<double> x_prox, vector<double> x_atual){
    return distanciaMaximo(x_prox, x_atual) / distanciaMaximo(x_prox, x_prox);
}

//Função para calcular a norma máxima de um vetor
double normaMaximo(vector<double> x){
    return distanciaMaximo(x, x);
}

//Função para imprimir um vetor
void printVetor(vector<double> x){
    for(int i = 0; i < x.size(); i++){
        cout << x[i] << " ";
    }
    cout << endl;
}

//Função para imprimir uma matriz
void printMatriz(vector<vector<double> > M){
    for(int i = 0; i < M.size(); i++){
        for(int j = 0; j < M[i].size
        (); j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

// Método de Gauss Seidel


vector<double> gauss_seidel(vector<vector<double> > M, vector<double> B, vector<double> u, double E, int max_iteracoes){
        
    int ordem = M[0].size();
    vector<double> X = u;
    vector<double> Xerro = u;//vetor pra calcular o erro
    int passos = 0;

    cout << "ordem da matriz " << ordem << endl;

    for(int k = 0; k < max_iteracoes; k++){
        
        //percorre a matriz
        for(int i = 0; i < ordem; i++){
            //comeca a soma pelo termo do vetor fonte
            double soma = B[i];
            double div = 0;
            for(int j = 0; j < ordem; j++){
                passos++;
                //separa o divisor
                if(i==j){
                    div = M[i][j];
                }
                else{
                    soma += M
                    [i][j] * X[j] * -1.0;
                }
            }
            //cria vetor de solucoes para proxima iteracao com resultados da linha
            X[i] = soma / div;
            
        }
        
        if(k % 100 == 0){
            cout << "Gauss-Seidel fez " << k << " iteracoes..." << endl;
        }
        
        //se atingir o criterio de parada, interrompe e retorna os resultados
        double erro = calculaErro(X,Xerro);
        
        if(erro < E){
            cout << "Terminou Gauss Seidel com erro de: " << erro << endl;
            return X;
        }

    }

    cout << "Gauss Seidel nao convergiu ou precisa de mais iteracoes para convergir" << endl;
    return X;
}

// Método Jacobi

vector<double> jacobi(vector<vector<double> > A, vector<double> B, vector<double> u, double E, int max_iteracoes){
    int ordem = A.size();
    vector<double> X = u;
    vector<double> Xp = u;
    int passos = 0;
    int iteracoes = 0;
    
    for(int k = 0; k < max_iteracoes; k++){
        iteracoes++;
        
        for(int i = 0; i < ordem; i++){
            double soma = 0.0;
            passos++;
            //passa linha pra um array
            vector<double> L = A[i];
            //zera o lugar onde seria o valor da diagonal principal
            L[i] = 0;
            //faz produto escalar entre os vetores
            for(int j = 0; j < ordem; j++){
                soma += L[j] * X[j];
            }
            
            Xp[i] = (1.0 / A[i][i])
            * (B[i] - soma);
            
        }
        
        //se atingir o criterio de parada, interrompe e retorna os resultados
        double erro = calculaErro(Xp, X); 

        if(erro < E){
            cout << "Terminou Jacobi com erro de: " << erro << endl;
            return Xp;
        }
        
        //prepara proxima iteracao com aproximacao da anterior
        X = Xp;
        
        //Imprime número de iterações de 100 em 100
        if(k % 100 == 0){
            cout << "Jacobi fez " << k << " iteracoes..." << endl;
        }
    
    }

    cout << "Jacobi nao convergiu ou precisa de mais iteracoes para convergir" << endl;
    return Xp;
}


int main(){
    int n, questao;
    clock_t begin, end;
    double time_spent;
    cout <<"Qual questão: ";
    cin>>questao;
    cout << "Digite o numero de equacoes: ";
    cin >> n;
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);
    vector<double> solucao(n);
    vector<double> xs; 
    int k = sqrt(n);
    
    if(questao == 1)
    {
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                A[i][j] =
                1.0/(i + j + 1);
            }
            B[i] = 1.0/(i + n + 1); 
            xs.push_back(1);
        }
    }
    else
    {
        //Todas as posições de B com 10 
        for(int i=0; i<n; i++)
        {
            A[i][i] = -4;
        }
        for(int i=0; i<n-1; i++)
        {
            A[i+1][i] = 1;
            A[i][i+1] = 1;
        
            if((i+1)%k == 0)
            {
                A[i+1][i] =0;
                A[i][i+1] =0; 
            }
        }
        for(int i=k; i<n; i++)
        {
            A[i-k][i] = 1;
            A[i][i-k] = 1; 
        }

        for(int i =0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                A[i][j] = -(A[i][j]*(k-1) * A[i][j]*(k-1));
            }
        }
    }

    int x;
    cout << "Digite o metodo desejado:\n1 - Gauss com pivoteamento \n2 - Gauss sem pivoteamento \n3 - Decomposicao LU \n4 - Cholesky\n5 - Gauss Seidel\n";
    cin >> x;

        switch (x){
        case 1:
            cout << "Gauss com pivoteamento" << endl;
            begin = clock();
            solucao = gauss_pivoteamento(A, B, n);
            end = clock();
            break;
        case 2:
            cout << "Gauss sem pivoteamento" << endl;
            begin = clock();
            solucao = gauss(A, B, n);
            end = clock();
            break;
        case 3:
            cout << "Decomposicao LU" << endl;
            begin = clock();
            solucao = decomposicao_LU(A, B, n);
            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            break;
        case 4:
            cout << "Cholesky" << endl;
            begin = clock();
            solucao = cholesky(A, B, n);
            end = clock();
            break;
        case 5:
            cout<<"Gauss Seidel"<<endl;
            solucao = gauss_seidel(A, B, xs, 0.001, 5000);
            break;
        case 6:
            cout<<"jacobi"<<endl;
            solucao = jacobi(A, B, xs, 0.001, 5000);
            break;

    xs.clear();
    /* cout << "Solucao: " << endl; */
    for(int j = 0; j < n; j++){
        /* cout << "x" << i
             << " = " << solucao[i]
             << "\n"; */
        xs.push_back(j);
    }

    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    /* cout << endl << "Erro: " << maxNorma(A, B, solucao, n) << endl; */
    cout << "Tempo de execucao: " << time_spent << endl;

    plotaGrafico(solucao, xs); 
    }
    

    return 0;
}