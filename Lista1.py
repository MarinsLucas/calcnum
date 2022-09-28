from telnetlib import TM
import matplotlib.pyplot as plt
import numpy as np
import cmath as math


def resfriamentoNewtonExata(thetaInicial, thetaM, k, tt):
    return ((thetaInicial - thetaM)*np.exp(-k*tt)) + thetaM


refin = 10
erro = np.zeros(refin)
vetorDeltaT = np.zeros(refin)
condicaoInicial = 99
coeficienteDifusao = 0.358
tf = 100.0
ti = 0.0
T_s = 27

for k in range(0, refin):

    quantPontos = 4**(k+1)+1
    dt = (tf-ti)/(quantPontos - 1)

    t = np.linspace(ti, tf, quantPontos)
    f = np.zeros(quantPontos)
    
    f[0] = condicaoInicial

    for n in range(0, quantPontos-1):
        f[n+1] = (1.0-k*dt*(1.0-0.5))/(1.0+k*dt*-0.5) * \
        f[n] + (k*dt*(1.0-0.5)*T_s)/(1.0+k*dt*0.5)

        tt = np.linspace(ti, tf, 200)

        plt.plot(t, f, '-*', tt, resfriamentoNewtonExata(condicaoInicial, T_s, coeficienteDifusao, tt))
        plt.legend(['aproximada', 'exata'])
        plt.show()
