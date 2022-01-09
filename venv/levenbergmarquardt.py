# Método de Levenberg-Marquardt

import math
import numpy as np
import matplotlib.pyplot as plt

# Parâmetros de entrada
x1 = 0.1
tol = 10 ** (-8)
h = 10 ** (-4)
lambida = 10 ** (-2)
analise = 4
analise1 = ['Sphere Function', 'Beale Function','Three Ramp Camel', 'Rastrigin', 'Default Problem' ]

print(f"Método de Levenberg Marquardt\n{30 * '-'}\nX_inicial {x1}\nTolerância: {tol}\nh: {h}\nLambda: {lambida}")

e = 50_000_000  # [N/cm^2]
i = 3_0000  # [cm^4]
le = 600  # [cm]
w = 2_500  # [N/cm^2]


def funcaox(x):
    if analise == 1:
        valor = sum([x ** 2])  # Sphere function
    elif analise == 2:
        valor = (1.5 - x) ** 2 + (2.25 - x) ** 2 + (2.625 - x ** 2)  # Beale Function
    elif analise == 3:
        valor = 2 * x ** 2 - 1.05 * x ** 4 + x ** 6 / 6  # Three ramp camel
    elif analise == 4:
        valor = 10 + sum([(x ** 2 - 10 * np.cos(2 * math.pi * x))])  # Rastrigin
    else:
        valor = (w / (120 * e * i * le)) * (- x ** 5 + 2 * le ** 2 * x ** 3 - le ** 4 * x)  # Função Default
    return valor


def derivadas(x):
    derivada1 = (funcaox(x + h) - funcaox(x - h)) / (2 * h)  # Derívada numérica - Diferenças finitas
    derivada2 = (funcaox(x + h) - 2 * funcaox(x) + funcaox(x - h)) / (h ** 2)  # Derivada numérica - D
    # derivada1 = (w / (120 * e * i * le)) * (-5 * x ** 4 + 6 * le ** 2 * x ** 2 - le ** 4)
    # derivada2 = (w / (120 * e * i * le)) * (-20 * x ** 3 + 12 * le ** 2 * x)
    return derivada1, derivada2


cp = 1
xk = [x1]
fun = [funcaox(x1)]
cont = 0
print(f'''Iteração xk{' '*10} f'(xk){' '*5} f''(xk){' '*5}Lambda{' '*5} x(k+1){' '*10} cp''')
lambida_velho = 0
beta = 2
poisson = 2
alfa = 1 / 2
while cp > tol:
    xk.append(xk[-1] - (derivadas(xk[-1])[0] / (derivadas(xk[-1])[1]) + lambida))
    cp = abs(xk[-1] - xk[-2])
    # lambida = abs(funcaox(xk[-1]) - funcaox(xk[-2])) / abs(funcaox(xk[0]))
    lambida_velho = lambida
    if lambida < lambida_velho:
        lambida = lambida * poisson
        poisson = 2 * poisson
    else:
        lambida = lambida * alfa
        poisson = beta
    lambida_velho = lambida
    cont += 1
    print(f"{cont}{' '*8}{xk[-2]:.3e}{' '*3}{derivadas(xk[-2])[0]:.3e}{' '*3}{derivadas(xk[-2])[1]:.3e}{' '*3}{lambida:.3e}{' '*3}{xk[-1]:.3e}{' '*3}{cp:.3e}")
    fun.append(funcaox(xk[-1]))

if xk[-1] < 10 ** -8:
    xk[-1] = 0

print(f'{"-" * 30}\nMínimo da função: {funcaox(xk[-1]):.3e}\nPonto correspondente: {xk[-1]:.3e}')

vetory = []
# yy = range(0, 600_000)  # Intervalo 1
yy = range(-5120, 5120)  # Intervalo 2
# yy = range(-600_000, 600_000)  # Intervalo 3
yyy = [abacate / 1000 for  abacate in yy]

for j in yyy:
    vetory.append(funcaox(j))

iteracoes = range(cont + 1)
plt.plot(yyy, vetory, xk, fun, 'o--')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'{analise1[analise - 1]}\n Levenberg Marquardt ')
plt.grid()
plt.savefig(f"{analise1[analise - 1]} Levenberg Marquardt.png",dpi=300)
plt.show()

