# Introdução
O método de integração numérica conhecido como Quadratura de Gauss-Legendre, ou simplesmente, Quadratura de Gauss busca obter uma aproximação para a uma integral através da criação de uma fórmula no seguinte formato:
    I = (integral de f(x) no intevalo a,b) = w_0 * f(x0) + ...  + w_n*f(x_n)
Onde os coeficientes w_i, também chamados nesse estudo de pesos, devem ser calculados de forma a obter a melhor precisão possível.

O objetivo do estudo a seguir é desenvolver um código que seja capaz de calcular pontos de integração da quadratura de Gauss para qualquer intervalo [a, b]. O programa desenvolvido deve retornar os pesos w_i e os pontos de integração x_i com i = 1,2...n

