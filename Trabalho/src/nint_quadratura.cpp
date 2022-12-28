#include <iostream>
#include <cmath>

using namespace std;

#define MIN_SIMPSON_PART 1500
//Definindo um nome para a variável de ponto flutuante para
//facilitar possíveis trocas de precisão de ponto flutuante
typedef long double fl;

//FIXME: essa função, está com precisão de 10e-20, porém era para estar exata. É erro de ponto flutuante!
fl simpson_3_8(fl a, fl b, int pot)
{
    fl I = 0; 
    //TODO: calcular simpson para no mínimo 1000 partições
    fl h = fl((b-a)/MIN_SIMPSON_PART);
    for(int i=0; i<MIN_SIMPSON_PART; i+=3)
    {
        I += (fl(3*h)/fl(8.0)) * (powl(a + i*h, pot) + 3*powl(a + (i+1)*h, pot) + 3*powl(a + (i+2)*h, pot) + powl(a +(i+3)*h, pot));
    }
    return I;
}
int main()
{
    //TODO: Definir o domínio de integração
    fl a, b;
    cout<<"Definir domínio de integração [a,b]: ";
    cin>>a>>b;

    //TODO: Definir o número de integração N (lembrando que N pontos = 2N incógnitas)
    int N;
    cout<<"Definir o número de pontos de integração N: ";
    cin>>N;

    //TODO: Definir condição inicial w_0 e t_0, com N entradas cada uma, para o método iterativo de Newton (relações descritas no roteiro)
    fl *w0 = new fl[N];
    fl *t0 = new fl[N];

    for(int i =0; i <N/2; i++)
    {
        w0[i] = fl((b-a)*(i+1))/fl(2*N);
        t0[i] = fl(a + (i+1)*(w0[i]/2));
        w0[N-1-i] = w0[i];
        t0[N-1-i] = (a+b) - t0[i];
    }
    
    if(N%2 == 1)
    {
        t0[N/2] = (a+b)/fl(2); 
        w0[N/2] =  fl((b-a)*(N/fl(2)))/fl(2*N); //FIXME: esse valor está correto?
    }

    //!para debuging
    /* for(int i =0; i < N;i++)
    {
        cout<<w0[i]<<" ";
    }*/

    cout.precision(12);

    //TODO: criar vetor com integrais de x^i para criar um sistema linear bitelo
    fl *g = new fl[2*N]; //g = vetor de integrais de x^i
    for(int i=0; i<2*N; i++)
    {
        g[i] = simpson_3_8(a, b, i);
        cout<<g[i]<<" ";
    }

    //TODO: Montar matriz Jacobiana formada pela derivada da função f(w,t) (descrito no roteiro)
    
    //TODO: Escolher uma tolerancia na ordem de 10e-8 ou menos!
}