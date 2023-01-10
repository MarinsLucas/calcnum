#Jonatas Dias Machado Costa - 202165077AC
#Lucas Marins Ramalho de Lima - 202165555C

import matplotlib.pyplot as plt
import numpy as np
import pprint

def calcErro(exata, aprox):
    return max(abs(exata - aprox))

def taxaConvergencia(erro, pontos):
    taxa = []
    for i in range(0, len(pontos)-1):
        taxa.append((np.log(erro[i+1]) - np.log(erro[i])) / (np.log(pontos[i+1]) - np.log(pontos[i])))
    return taxa

def cooling(T0, k, T_m, dt, npoints, type):
    T = np.zeros(npoints)
    T[0] = T0
    for i in range(npoints-1):
        if type == 1:
            T[i+1] = -dt*k*(T[i] - T_m) + T[i]  # Euler Explícito

        elif type == 2:
            T[i+1] = (T[i] + k * dt * T_m)/(1+k*dt)  # Euler Implícito

        elif (type == 3):  
            T[i+1] = (-k/2 *dt*(T[i] - 2*T_m) + T[i])/(1 + (k*dt)/2)
    return T

def writeData(x, exata, eulere, euleri, cranckn):
    f = open("datafile.csv", "a")
    f.write("novoteste;novoteste;novoteste;novoteste;novoteste\n")
    f.write("TEMPO;EXATA;EULER EXPLÍCITO;EULER IMPLÍCITO;CRANCK NICOLSON\n")
    for i in range(0, len(exata)):
        f.write(str(x[i]) + "; " + str(exata[i]) + "; " + str(eulere[i]) + "; " + str(euleri[i]) + "; " + str(cranckn[i]) + '\n')
    f.close()

def eraseData():
    f = open("datafile.csv", "w")
    f.write('\0')
    f.close()

def main():
    #dados do enunciado 
    T0 = np.float16(99)   #temperatura inicial
    k = np.float16(0.0358)  #coeficiente de transferência de calor 
    T_m = np.float16(27)  #temperatura final
    t_ini = np.float16(0) #tempo inicial
    t_end = np.float16(50) #tempo final
    ref = 5
    vPoints = []
    erroExp = []
    erroImp = []
    erroCN = []
    
    erroEe = []
    erroEi = []
    erroCn = []
    dtt = []
    eraseData()

    for it in range(0, ref):
        npoints = 4**(it+1)
        vPoints.append(npoints)

        t = np.float128(np.linspace(t_ini, t_end, npoints))
        dt = (t_end-t_ini)/(npoints-1)
        dtt.append(dt)

        T = cooling(T0, k, T_m, dt, npoints, 1)
        print("Exata Eulere: " + str(calcErro(((T0 - T_m)*np.exp(-k*t)  + T_m), T)))
        plt.plot(t, T, '-*')
        erroEe.append(calcErro(((T0 - T_m)*np.exp(-k*t)  + T_m), T))
        
        T = cooling(T0, k, T_m, dt, npoints, 2)
        print("Exata Euleri: " + str(calcErro(((T0 - T_m)*np.exp(-k*t)  + T_m), T)))
        plt.plot(t, T, '-o')
        erroEi.append(calcErro(((T0 - T_m)*np.exp(-k*t)  + T_m), T))

        T = cooling(T0, k, T_m, dt, npoints, 3)
        print("Exata Cranck-Nicolson: " + str(calcErro(((T0 - T_m)*np.exp(-k*t)  + T_m), T)))
        plt.plot(t, T, '-^')
        erroCn.append(calcErro(((T0 - T_m)*np.exp(-k*t)  + T_m), T))


        tt = np.linspace(t_ini, t_end, 200)
        # plt.plot(tt, ((T0 - T_m)*np.exp(-k*tt)) + T_m, '--')
        
        #print("Delta Time " + str(dt))

        """ plt.legend(['Euler Explícito', 'Euler Implícito',
                'Crank-Nicolson', 'Exato']) """

        """ plt.xlabel('tempo (s)')
        plt.ylabel('temperatura (C)')
        plt.title('Resfriamento de Newton')
        plt.grid()
        plt.show() """

        writeData(t, exata, cooling(T0, k, T_m, dt, npoints, 1),\
             cooling(T0, k, T_m, dt, npoints, 2), cooling(T0, k, T_m, dt, npoints, 3))
    plt.plot(-np.log(dtt),np.log(erroEe), '-*')
    plt.plot(-np.log(dtt),np.log(erroEi),'-o')
    plt.plot(-np.log(dtt),np.log(erroCn),'-^')
    plt.grid()
    plt.legend(['Euler Explícito', 'Euler Implícito','Crank-Nicolson'])
    plt.xlabel('-log(h)')
    plt.ylabel('log(erro)')
    print("Taxas de convergência:")
    print("Euler Explícito " + str((np.log(erroEe[-1])-np.log(erroEe[0]))/(np.log(dtt[-1])-np.log(dtt[0]))))
    print("Euler Implícito " + str((np.log(erroEi[-1])-np.log(erroEi[0]))/(np.log(dtt[-1])-np.log(dtt[0]))))
    print("Cranck-Nicolson " + str((np.log(erroCn[-1])-np.log(erroCn[0]))/(np.log(dtt[-1])-np.log(dtt[0]))))

    plt.show()
    
if __name__ == "__main__":
    main()
