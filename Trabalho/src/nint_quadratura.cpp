#include <iostream>
#include <cmath>
#include "sislin.h"

using namespace std;

#define MIN_SIMPSON_PART 1500
#define debug std::cout
#define epslon 1E-8;

// Definindo um nome para a variável de ponto flutuante para
// facilitar possíveis trocas de precisão de ponto flutuante
typedef long double flu;

// variáveis globais
int N;

// FIXME: essa função, está com precisão de 10e-20, porém era para estar exata. É erro de ponto flutuante!
flu simpson_3_8(flu a, flu b, int pot)
{
    flu I = 0;
    // TODO: calcular simpson para no mínimo 1000 partições
    flu h = flu((b - a) / MIN_SIMPSON_PART);
    for (int i = 0; i < MIN_SIMPSON_PART; i += 3)
    {
        I += (flu(3 * h) / flu(8.0)) * (powl(a + i * h, pot) + 3 * powl(a + (i + 1) * h, pot) + 3 * powl(a + (i + 2) * h, pot) + powl(a + (i + 3) * h, pot));
    }
    return I;
}

// retorna um valor para f_j
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
                flu aux = t[i] + epslon;
                f[j] = w[i] + pow(aux, j);
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
                f[j] += w[i];
                f[j] += epslon;
                f[j] *= pow(t[i], j);
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

// TODO: função que cria matriz jacobiana
flu **Jacobiana(flu* w, flu* t, flu* g)
{
    flu** jaco = new flu*[2*N];
    for (int j = 0; j < 2 * N; j++)
    {
        jaco[j] = new flu[2*N];
        //percorre as colunas
        for(int i = 0; i < N; i++)
        {
            //soma-se epslon na iesima posição de w
            jaco[j][i] = ((f2(w, t, g, i)[j] - f(w, t, g)[j]))/epslon;
            //debug<<jaco[j][i]<<endl;

        }
        for(int i = N; i<2*N ;i++)
        {
            //soma-se epslon na iesima posição de t
            jaco[j][i] = ((f3(w, t, g, i-N)[j] - f(w, t, g)[j]))/epslon; 
            //debug<<jaco[j][i]<<endl;

        }
    }
    return jaco;
}

// TODO: resolução de sistema linear pelo método de newton
void newtonSisNLin(flu *w, flu *t, flu *g) // J é a matriz jacobiana, e f é o sistema de funções não lineares
{
    int kmax = 100; 
    flu** j; 
    flu* fntion;
    vector<flu> solution;
    // Processo iterativo:
    for (int k = 0; k < kmax; k++)
    {
        fntion = f4(w, t, g);
        j = Jacobiana(w, t, g);
        solution = gauss_pivoteamento(j, fntion, 2*N);

        /* for(int i=0;i<2*N; i++)
        {
            for(int k=0 ; k<2*N; k++)
            {
                debug<<j[i][k]<<";";
            }
            debug<<endl;
        } */
        
        /* for(int i = 0; i<2*N; i++)
        {
            debug<<solution[i]<<endl;
        }    
        debug<<endl;  */

        for(int i = 0; i<N; i++)
        {
            w[i] = solution[i]*w[i] + w[i];
            t[i] = solution[i]*t[i] + t[i]; 
        } 
    }
}

int main()
{
    // TODO: Definir o domínio de integração
    flu a, b;
    cout << "Definir domínio de integração [a,b]: ";
    cin >> a >> b;

    // TODO: Definir o número de integração N (lembrando que N pontos = 2N incógnitas)
    cout << "Definir o número de pontos de integração N: ";
    cin >> N;

    // TODO: Definir condição inicial w_0 e t_0, com N entradas cada uma, para o método iterativo de Newton (relações descritas no roteiro)
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

    //! para debuging
    /* for(int i =0; i < N;i++)
    {
        cout<<w0[i]<<" ";
    }*/

    cout.precision(12);

    // TODO: criar vetor com integrais de x^i para criar um sistema linear bitelo
    flu *g = new flu[2 * N]; // g = vetor de integrais de x^i
    for (int i = 0; i < 2 * N; i++)
    {
        g[i] = simpson_3_8(a, b, i);
    }

    // TODO: Montar sistema linear de funções f(w,t)
    /* for (int j = 0; j < 2 * N; j++)
    {
        debug << f(w0, t0, g)[j] << endl;
    } */


    // TODO: Montar matriz Jacobiana formada pela derivada da função f(w,t) (descrito no roteiro)
    newtonSisNLin(w0, t0, g);

    // TODO: Escolher uma tolerancia na ordem de 10e-8 ou menos!
}