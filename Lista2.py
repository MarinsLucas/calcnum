#Lista 1 - Parte 2
#Jonatas Dias Machado Costa - 202165077AC
#Lucas Marins Ramalho de Lima - 202165555C

from math import e
from cmath import sqrt
import matplotlib.pyplot as plt
import numpy as np

""" def TDMASolve(a, b, c, d): 
    #a,b,c são as diagonais 
    #d é onde vamos armazenar a resposta
    
    n = len(d)
    nmax = n #??

    # Modifica o primeiro coeficiente de cada linha
    c[0] /= b[0]
    d[0] /= b[0]
    
    # Modifica os coeficientes de cada linha
    for i in range(1, nmax):
        ptemp = b[i] - (a[i] * c[i-1]) 
        c[i] /= ptemp
        d[i] = (d[i] - a[i] * d[i-1])/ptemp

    #Substituição de volta
    x = [0 for i in range(nmax)]
    x[-1] = d[-1]
    for i in range(-2,- nmax-1,-1): 
        x[i] = d[i] - c[i] * x[i+1] 

    return x """

def solExata(epslon):
    t_ini = 0
    t_final = 1
    c2 =  (pow(e, (-1/sqrt(epslon))) - 1) / (pow(e, (1/sqrt(epslon))) - pow(e, (-1/sqrt(epslon))))
    c1 = - 1 - c2
    x = np.linspace(t_ini, t_final, 200)
    ux = c1* pow(e, (-x/sqrt(epslon)))  +  c2*pow(e, (x/sqrt(epslon))) + 1
    plt.plot( x, ux, label = "Solução exata")
    #plt.show()


def difFinita(epslon, h, npoints):
    U = np.zeros(npoints)
    U[0] = 0
    

    for i in range(1 , npoints-1):
        U[i+1] = -U[i-1] + 2*U[i] - (((1 - U[i])* (h * h))/epslon)  

    x = np.float128(np.linspace(0, 1, npoints))
    
    plt.plot(x, U, '-*')
    plt.show()

def main():
    solExata(0.1)   
    difFinita(0.1, 1/99, 100)


main()