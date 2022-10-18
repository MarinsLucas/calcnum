#Lista 1 - Parte 2
#Jonatas Dias Machado Costa - 202165077AC
#Lucas Marins Ramalho de Lima - 202165555C

from math import e
from cmath import sqrt
import matplotlib.pyplot as plt
import numpy as np

def solExata(eu):
    t_ini = 0
    t_final = 1
    c2 =  (pow(e, -(1/sqrt(eu))) - 1) / (pow(e, (1/sqrt(eu))) - pow(e, -1/sqrt(eu)))
    c1 = 0 - 1 - c2
    x = np.linspace(t_ini, t_final, 200)
    plt.plot( x, c1* pow(e, -(x/sqrt(eu)))  +  c2*pow(e, (x/sqrt(eu))+ 1) )
    plt.show()

def main():
    solExata(0.1)


main()