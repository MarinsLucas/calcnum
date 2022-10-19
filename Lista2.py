#Lista 1 - Parte 2
#Jonatas Dias Machado Costa - 202165077AC
#Lucas Marins Ramalho de Lima - 202165555C

from math import e
from cmath import sqrt
import matplotlib.pyplot as plt
import numpy as np

def solExata(epslon):
    t_ini = 0
    t_final = 1
    c2 =  (pow(e, (-1/sqrt(epslon))) - 1) / (pow(e, (1/sqrt(epslon))) - pow(e, (-1/sqrt(epslon))))
    c1 = - 1 - c2
    x = np.linspace(t_ini, t_final, 200)
    ux = c1* pow(e, (-x/sqrt(epslon)))  +  c2*pow(e, (x/sqrt(epslon))) + 1
    plt.plot( x, ux, label = "Solução exata")
    plt.show()

def main():
    solExata(0.1)

main()