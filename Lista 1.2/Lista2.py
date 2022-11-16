#Lista 1 - Parte 2
#Jonatas Dias Machado Costa - 202165077AC
#Lucas Marins Ramalho de Lima - 202165555C

from math import e
from cmath import sqrt
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

def eraseData():
    f = open("datafile.csv", "w")
    f.write('\0')
    f.close()

def writeData(x, exata, listaSolucoes, todosErros):
    f = open("float128.csv", "a")
    f.write("novoteste;novoteste;novoteste;novoteste;novoteste;novoteste\n")
    f.write("TEMPO;EXATA;0.1;0.01;0.001;0.0001\n")
    for i in range(0, len(listaSolucoes[0])):
        f.write(str(x[i]) + "; " + str(exata[i]) + "; " + str(listaSolucoes[0][i]) + "; " + str(listaSolucoes[1][i]) + "; " + str(listaSolucoes[2][i]) + "; " + str(listaSolucoes[3][i]) + '\n')
    f.close()

def normaMaximo(x):
    size = len(x)
    maximo = np.abs(x[0])   
    
    for i in range(size):
        temp = np.abs(x[i])        
        if(temp > maximo):
            maximo = temp
            
    return maximo

def distanciaMaximo(x1, x2):
    if(len(x1) != len(x2)):
        print("O tamanho dos vetores x1 e x2 precisa ser o mesmo")
        return 0
        
    size = len(x1)
    dist = abs(x1[0] - x2[0])  
    
    for i in range(size):
        temp = abs(x1[i] - x2[i])
        if(temp > dist):
            dist = temp
            
    return dist

def calculaErro(x_prox, x_atual):
    return distanciaMaximo(x_prox, x_atual)/normaMaximo(x_prox)


def Tomas(a, b, c, d):
    
    #Print para Debug
    #print('diagonais: ', a, '\n', b, '\n', c, '\n',d)
    
    n = len(d) 
    c_ = [ c[0] / b[0] ]
    d_ = [ d[0] / b[0] ]
    
    for i in range(1, n):
        aux = b[i] - c_[i-1]*a[i-1]
        if i < n-1:
            c_.append( c[i] / aux )
        d_.append( (d[i] - d_[i-1]*a[i-1])/aux )
    
    # Substituicao de volta
    x = [d_[-1]]
    for i in range(n-2, -1, -1):
        x = [ d_[i] - c_[i]*x[0] ] + x
    
    return x

def F(x, epslon):
    c2 =  (pow(e, (-1/sqrt(epslon))) - 1) / (pow(e, (1/sqrt(epslon))) - pow(e, (-1/sqrt(epslon))))
    c1 = - 1 - c2
    ux = c1* pow(e, (-x/sqrt(epslon)))  +  c2*pow(e, (x/sqrt(epslon))) + 1
    return ux

def difFinita(epslon, h, npart, it):
    
    ordem = npart -1
    diagonalSuperior = []
    diagonalPrincipal = []
    diagonalInferior = []
    termoIndependente = []
    solApprox = []
    tempoExata = []
    tempo = []
    solExata = []
    todosErros = []
    listaSolucoes = []
    
    for b in range(4):

        diagonalInferior.clear()
        diagonalSuperior.clear()
        diagonalPrincipal.clear()
        termoIndependente.clear()
        tempoExata.clear()
        solExata.clear()
        tempo.clear()

        for i in range(ordem - 1):
            diagonalSuperior.append(epslon[b] * -1.0)
            diagonalInferior.append(epslon[b]  * -1.0)
        
        for i in range(ordem):
            diagonalPrincipal.append(2 * epslon[b]  + (h*h))
            termoIndependente.append(h*h)
        
        particoes_exata = 100
        
        he = np.float128(1/float(particoes_exata-1))
        for i in range(particoes_exata):
            dhe = np.float128(he * i)
            tempoExata.append(dhe)
            solExata.append(F(dhe, epslon[b]))

        for i in range(npart + 1):
            dh = np.float128(h * i)
            tempo.append(dh)
        
        solApprox = Tomas(diagonalInferior, diagonalPrincipal, diagonalSuperior, termoIndependente)
        
        solApprox.insert(0, 0)
        solApprox.append(0)
        listaSolucoes.append(solApprox)

        partErro = len(solApprox)
        hErro = np.float128(1/float(partErro-1))
        solExataErro = []
        
        for i in range (partErro):
            dh = np.float128(hErro * i)
            solExataErro.append(F(dh, epslon[b]))
             
        erroNormaMax = calculaErro(solExataErro, solApprox)
        print("Erro Norma Max: ", repr(erroNormaMax))
        
        erro = []
        for j in range(len(solExataErro)):
            erro.append(abs(solExataErro[j] - solApprox[j]))
        todosErros.append(erro)

        plt.plot(tempoExata, solExata, '-',tempo, solApprox, '-.')
        
        plt.ylabel(u"Valor de u(h)")
        plt.xlabel(u"Valor de h, " + str(npart) + u" partições")
        
        se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução Exata', linestyle='-')
        ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Solução Aprox.', linestyle= '-.')
        
        plt.legend(handles=[se_line, ac_line])
        
        plt.title("Metodos de Resolucao " + "ε =" + str(epslon[b]))
        plt.show() 

    plt.plot(tempo, todosErros[0]) #erro 0.1
    plt.plot(tempo, todosErros[1]) #erro 0.01
    plt.plot(tempo, todosErros[2]) #erro 0.001
    plt.plot(tempo, todosErros[3]) #erro 0.0001
    
    plt.legend(['e = 0.1', 'e = 0.01', 'e = 0.001', 'e = 0.0001'])
    plt.ylabel(u"Valor de erro")
    plt.xlabel(u"Valor de h, " + str(npart) + u" partições")
    plt.title("i = " + str(it))

    plt.show()
    
    writeData(tempo, solExataErro, listaSolucoes, todosErros)

def main(): 
    eraseData()
    e = [np.float128(0.1) , np.float128(0.01), np.float128(0.001), np.float128(0.0001)]
    h = 4
    for i in range(1,6):
        difFinita(e, np.float128(1/h**i), h**i, i)

main()

    



