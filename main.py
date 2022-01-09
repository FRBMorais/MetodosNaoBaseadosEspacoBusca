# Método de Levenberg-Marquardt

import math
import numpy as np

# Parâmetros de entrada
x1 = 300
tol = 10 ** (-8)
h = 10 ** (-4)
lambida = 1
analise = 5

print(f"X_inicial {x1}\nTolerância: {tol}\nh: {h}\nLambda: {lambida}")

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
    derivada2 = (funcaox(x + h) - 2 * funcaox(x) + funcaox(x - h)) / (h ** 2)
    derivada1 = (funcaox(x + h) - funcaox(x - h)) / (2 * h)
    # derivada1 = (w / (120 * e * i * le)) * (-5 * x ** 4 + 6 * le ** 2 * x ** 2 - le ** 4)
    # derivada2 = (w / (120 * e * i * le)) * (-20 * x ** 3 + 12 * le ** 2 * x)
    return derivada1, derivada2


cp = 1
xk = [x1]
cont = 0
print(f'''Iteração xk{' '*10} f'(xk){' '*5} f''(xk){' '*5}Lambda{' '*5} x(k+1){' '*10} cp''')
while cp > tol:
    xk.append(xk[-1] - (derivadas(xk[-1])[0] / (derivadas(xk[-1])[1]) + lambida))
    cp = abs(xk[-1] - xk[-2])
    lambida = abs(funcaox(xk[-1]) - funcaox(xk[-2])) / abs(funcaox(xk[0]))
    cont += 1
    print(f"{cont}{' '*8}{xk[-2]:.3e}{' '*3}{derivadas(xk[-2])[0]:.3e}{' '*3}{derivadas(xk[-2])[1]:.3e}{' '*3}{lambida:.3e}{' '*3}{xk[-1]:.3e}{' '*3}{cp:.3e}")

print(f'Mínimo da função: {funcaox(xk[-1])}\nAvaliada no ponto: {xk[-1]}')