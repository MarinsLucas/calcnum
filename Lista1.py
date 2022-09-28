from telnetlib import TM
import matplotlib.pyplot as plt
import numpy as np
import cmath as math

def resfriamentoNewtonExata():
    thetaInicial = 99
    thetaM = 27
    dt = 0.01
    k = 0.035871952
    x = np.arange(0, 500, dt)
    y = ((thetaInicial - thetaM)*np.exp(-k*x)) + thetaM        
    return(x, y)

def resfriamentoNewtonEulerExplicito(x,y):
    refin = 5
    erro = np.zeros(refin)
    vetorDeltaT = np.zeros(refin)
    condicaoInicial = 99
    coeficienteDifusao = 0.1
    tf = 500
    ti = 0.0

    for k in range(0,refin):
        quantPontos = 4**(k+1)+1
        dt = (tf-ti)/(quantPontos - 1)
        t = np.linspace(ti,tf,quantPontos)
        f = np.zeros(quantPontos) 
        f[0] = condicaoInicial

        for i in range(0, quantPontos-1):
            f[i+1] = (1.0/(1.0+coeficienteDifusao*dt))*f[i] 
        plt.plot(t,f)
        plt.plot(x,y)
        plt.legend(['aproximada','exata'])
        plt.show()


x,y = resfriamentoNewtonExata()
resfriamentoNewtonEulerExplicito(x,y)

plt.plot(x, y)
plt.show()