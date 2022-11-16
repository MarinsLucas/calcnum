import numpy as np
import matplotlib.pyplot as plt

def uexata(a,b,t): #função exata
  return b*np.exp(-a*t) # Solução exata

ti = 0.0 # Tempo inicial
tf = 1.0 # Tempo final
refin = 5 # Refinamento

erro = np.zeros(refin) # Vetor para armazenar os erros
dtt = np.zeros(refin) # Vetor para armazenar os passos de tempo

for k in range(0,refin): # Loop para refinamento
    
  npoints = 4**(k+1)+1 # Número de pontos
  dt = (tf-ti)/(npoints-1) # Passo de tempo
  
  a = 1.0 # Coeficiente de difusão
  b = 1.0 # Condicao inicial
  
  # coordenadas
  t = np.linspace(ti,tf,npoints) # Vetor de tempo
  u = np.zeros(npoints) # Vetor para armazenar a solução
  u[0] = b # Condicao inicial

  for n in range(0,npoints-1): # Loop para resolução
    u[n+1] = (1.0/(1.0+a*dt))*u[n] # Solução numérica
    # u[n+1] = ((2.-a*dt)/(2.0+a*dt))*u[n]

  tt = np.linspace(ti,tf,200) # Vetor de tempo para a solução exata
  
  plt.plot(t,u,'-*',tt,uexata(a,b,tt))  # Plotando a solução numérica e exata
  plt.legend(['aproximada','exata']) # Legenda
  plt.show() # Mostrando o gráfico
  
  erro[k] = np.max(np.abs(u-uexata(a,b,t))) # Calculando o erro
  dtt[k] = dt # Armazenando o passo de tempo
  
  print(dt) # Mostrando o passo de tempo

plt.plot(-np.log(dtt),np.log(erro)) # Plotando o erro em função do passo de tempo
plt.show() # Mostrando o gráfico

print((np.log(erro[-1])-np.log(erro[0]))/(np.log(dtt[-1])-np.log(dtt[0]))) # Ordem de convergência