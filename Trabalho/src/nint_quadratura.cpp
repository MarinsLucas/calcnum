#include <iostream>
#include <cmath>
#include "sislin.h"

using namespace std;

#define MIN_SIMPSON_PART 1000
#define debug std::cout

// Definindo um nome para a variável de ponto flutuante para
// facilitar possíveis trocas de precisão de ponto flutuante
typedef long double flu;

// variáveis globais
int N;
flu epslon = 1e-8; 

// FIXME: essa função, está com precisão de 10e-20, porém era para estar exata. É erro de ponto flutuante!
flu simpson_3_8(flu a, flu b, int pot)
{
    flu I = 0;
    flu h = flu((b - a) / MIN_SIMPSON_PART);
    for (int i = 0; i < MIN_SIMPSON_PART; i += 3)
    {
        I += (flu(3 * h) / flu(8.0)) * (powl(a + i * h, pot) + 3 * powl(a + (i + 1) * h, pot) + 3 * powl(a + (i + 2) * h, pot) + powl(a + (i + 3) * h, pot));
    }

    //FIXME: TIRAR ESSE TREM QUANDO TIVER A EXATA
    if(I>0 && I<1E-19)
        I = 0; 
    return I;
}

//TODO: fazer integral exata
flu integralexata( flu min, flu max, int pot)
{
    return pow(max, pot+1)/(pot+1) - pow(min, pot+1)/(pot+1);
}

#pragma region Funções de cálculo de F
//! retorna um valor para f_j
flu* f(flu *w, flu *t, flu *g)
{
    flu *f = new flu[2*N];
    for(int j = 0; j<2*N; j++)
    {
        f[j] = 0; 
        for (int i = 0; i < N; i++)
        {
            f[j] += w[i] * pow(t[i], j);
        }
        f[j] -= g[j];
    }
    return f;
}

flu* f4(flu *w, flu *t, flu *g)
{
    flu *f = new flu[2*N];
    for(int j = 0; j<2*N; j++)
    {
        f[j] = 0; 
        for (int i = 0; i < N; i++)
        {
            f[j] += w[i] * pow(t[i], j);
        }
        f[j] -= g[j];
        f[j] *= -1; 
    }
    return f;
}

flu* f3(flu *w, flu *t, flu *g, int k)
{
    flu *f = new flu[2*N];
    for(int j = 0; j<2*N; j++)
    {
        f[j] = 0; 
        for (int i = 0; i < N; i++)
        {
            if(i == k)
            {
                f[j] += w[i] * pow(t[i] + epslon, j);
            }
            else
            {
                f[j] += w[i] * pow(t[i], j);
            }
        }
        f[j] -= g[j];
    }
    return f;
}

flu* f2(flu *w, flu *t, flu *g, int k)
{
    flu *f = new flu[2*N];
    for(int j = 0; j<2*N; j++)
    {
        f[j] = 0; 
        for (int i = 0; i < N; i++)
        {
            if(i == k)
            {
                f[j] += (w[i]+epslon)*pow(t[i], j);
            }
            else
            {
                f[j] += w[i] * pow(t[i], j);
            }
        }
        f[j] -= g[j];
    }
    return f;
}
#pragma endregion

// FIXME: função que cria matriz jacobiana
flu **Jacobiana(flu* w, flu* t, flu* g)
{
    flu** jaco = new flu*[2*N];
    for(int j =0 ; j< 2*N; j++)
    {
        jaco[j] = new flu[2*N];
        for(int i = 0; i<N; i++)
        {
            jaco[j][i] = f2(w, t, g, i)[j]/epslon - f(w, t, g)[j]/epslon;
        }
        for(int i = N; i<2*N; i++)
        {
            jaco[j][i] = f3(w,t,g,i-N)[j]/epslon - f(w, t, g)[j]/epslon;
        }
    }
    return jaco;
}

//! resolução de sistema linear pelo método de newton
void newtonSisNLin(flu *w, flu *t, flu *g) // J é a matriz jacobiana, e f é o sistema de funções não lineares
{
    int kmax = 1000;  //TODO: Definir tolerancia como 1e-8
    flu** j = new flu*[2*N]; 
    flu* fntion;
    vector<flu> solution(2*N);
    // Processo iterativo:
    for (int k = 0; k < kmax; k++)
    {
        fntion = f4(w, t, g);
        j = Jacobiana(w, t, g);
        solution = gauss(j, fntion, 2*N);

        for(int i = 0; i<N; i++)
        {
            w[i] = w[i] + solution[i];
            t[i] = t[i] + solution[N+i]; 
        } 
    }
    debug<<"resultado w: "<<endl;
    printVector(w, N);
    debug<<"resultado t: "<<endl;
    printVector(t,N);
}

int main()
{
    //! Definir o domínio de integração
    flu a, b;
    cout << "Definir domínio de integração [a,b]: ";
    cin >> a >> b;

    //! Definir o número de integração N (lembrando que N pontos = 2N incógnitas)
    cout << "Definir o número de pontos de integração N: ";
    cin >> N;

    //! Definir condição inicial w_0 e t_0, com N entradas cada uma, para o método iterativo de Newton (relações descritas no roteiro)
    flu *w0 = new flu[N];
    flu *t0 = new flu[N];

    for (int i = 0; i < N / 2; i++)
    {
        w0[i] = flu((b - a) * (i + 1)) / flu(2 * N);
        t0[i] = flu(a + (i + 1) * (w0[i] / 2));
        w0[N - 1 - i] = w0[i];
        t0[N - 1 - i] = (a + b) - t0[i];
    }

    if (N % 2 == 1)
    {
        t0[N / 2] = (a + b) / flu(2);
        w0[N / 2] = flu((b - a) * (N / flu(2))) / flu(2 * N); // FIXME: esse valor está correto?
    }

    //! criar vetor com integrais de x^i para criar um sistema linear bitelo
    flu *g = new flu[2 * N]; // g = vetor de integrais de x^i
    for (int i = 0; i < 2 * N; i++)
    {
        //g[i] = simpson_3_8(a, b, i);
        g[i] = integralexata(a, b, i);
    }

    //! Montar matriz Jacobiana formada pela derivada da função f(w,t) (descrito no roteiro)
    newtonSisNLin(w0, t0, g);
}
