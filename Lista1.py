import matplotlib.pyplot as plt
import numpy as np

def cooling(T0, k, T_m, dt, npoints, type):
    k *= -1
    T = np.zeros(npoints)
    T[0] = T0
    for i in range(npoints-1):
        if type == 1:
            T[i+1] = dt*k*(T[i] - T_m) + T[i]  # Euler Explícito

        elif type == 2:
            T[i+1] = (T[i] - k * dt * T_m)/(1-k*dt)  # Euler Implícito

        elif (type == 3):  
            T[i+1] = T[i] + dt*k * (T[i] + (1/2)*dt*k*(T[i] - T_m) - T_m)  # Crank-Nicolson

    return T

def writeData(x, exata, eulere, euleri, cranckn):
    f = open("datafile.csv", "a")
    f.write("novoteste;novoteste;novoteste;novoteste;novoteste\n")
    f.write("TEMPO;EXATA;EULER EXPLÍCITO;EULER IMPLÍCITO;CRANCK NICOLSON\n")
    for i in range(0, len(exata)):
        f.write(str(x[i]) + "; " + str(exata[i]) + "; " + str(eulere[i]) + "; " + str(euleri[i]) + "; " + str(cranckn[i]) + '\n')

    f.close()


def main():
    #dados do enunciado b)
    T0 = 99   #temperatura inicial
    k = 0.0358  #coeficiente de transferência de calor 
    T_m = 27  #temperatura final
    t_ini = 0 #tempo inicial
    t_end = 100 #tempo final
    ref = 5
    
    for it in range(0, ref):
        npoints = 3**(it+1)+1

        t = np.linspace(t_ini, t_end, npoints)
        dt = (t_end-t_ini)/(npoints-1)

        T = cooling(T0, k, T_m, dt, npoints, 1)
        plt.plot(t, T, '-*')
        T = cooling(T0, k, T_m, dt, npoints, 2)
        plt.plot(t, T, '-o')
        T = cooling(T0, k, T_m, dt, npoints, 3)
        plt.plot(t, T, '-^')

        tt = np.linspace(t_ini, t_end, 200)
        plt.plot(tt, ((T0 - T_m)*np.exp(-k*tt)) + T_m, '--')

        plt.legend(['Euler Explícito', 'Euler Implícito',
                'Crank-Nicolson', 'Exato'])

        plt.xlabel('tempo (s)')
        plt.ylabel('temperatura (C)')
        plt.title('Resfriamento de Newton')
        plt.grid()
        plt.show()

        writeData(t, ((T0 - T_m)*np.exp(-k*t)  + T_m)  , cooling(T0, k, T_m, dt, npoints, 1),\
             cooling(T0, k, T_m, dt, npoints, 2), cooling(T0, k, T_m, dt, npoints, 3))


if __name__ == "__main__":
    main()
