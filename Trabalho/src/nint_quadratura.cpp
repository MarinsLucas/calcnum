#include <iostream>
#include <cmath>

using namespace std;

//Definindo um nome para a variável de ponto flutuante para
//facilitar possíveis trocas de precisão de ponto flutuante
typedef float fl;

int main()
{
    //TODO: Definir o domínio de integração
    int a, b;
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

    //!para debuging
    /*for(int i =0; i < N;i++)
    {
        cout<<t0[i]<<" ";
    } */

    //TODO: Montar matriz Jacobiana formada pela derivada da função f(w,t) (descrito no roteiro)
    //TODO: Escolher uma tolerancia na ordem de 10e-8 ou menos!
    cout<<"hello world"<<endl;
}