import matplotlib.pyplot as plt
import math
import numpy as np


# Parâmetros de entrada
x1 = 0.1
tol = 10 ** (-8)
h = 10 ** (-4)


print(f"Método de Newton Raphson\n{30 * '-'}\nX_inicial {x1}\nTolerância: {tol}\nh: {h}\n{30*'-'}")

e = 50_000_000  # [KN/cm^2]
i = 30_000  # [cm^4]
le = 600  # [cm]
w = 2_500  # [KN/cm]
analise = 1
analise1 = ['Sphere Function', 'Beale Function','Three Ramp Camel', 'Rastrigin', 'Default Problem' ]


def funcaox(x):
    if analise == 1:
        valor = sum([x ** 2])  # Sphere function
    elif analise == 2:
        valor = (1.5 - x) ** 2 + (2.25 - x) ** 2 + (2.625 - x ** 2)  # Beale Function
    elif analise == 3:
        valor = 2 * x ** 2 - 1.05 * x ** 4 + x ** 6 / 6  # Three ramp camel
    elif analise == 4:
        valor = 10 + sum([(x**2 - 10 * np.cos(2 * math.pi * x))])  # Rastrigin
    else:
        valor = (w / (120 * e * i * le)) * (- x**5 + 2 * le**2 * x**3 - le**4 * x)  # Função Default
    return valor


def derivadas(x):
    derivada1 = (funcaox(x + h) - funcaox(x - h)) / (2 * h)  # Derívada numérica - Diferenças finitas
    derivada2 = (funcaox(x + h) - 2 * funcaox(x) + funcaox(x - h)) / (h ** 2)  # Derivada numérica - Diferenças finitas
    # derivada1 = (w / (120 * e * i * le)) * (-5 * x ** 4 + 6 * le ** 2 * x ** 2 - le ** 4)  # Derivadas analíticas
    # derivada2 = (w / (120 * e * i * le)) * (-20 * x ** 3 + 12 * le ** 2 * x)  # Derivada analítica
    return derivada1, derivada2


cp = 1
xk = [x1]
fun = [funcaox(x1)]
cont = 0
print(f'''Iteração xk{' '*10} f'(xk){' '*5} f''(xk){' '*10} x(k+1){' '*10} cp''')


# Métod numérico de aproximação
while cp > tol:
    xk.append(xk[-1] - (derivadas(xk[-1])[0] / derivadas(xk[-1])[1]))
    cp = abs(xk[-1] - xk[-2])
    cont += 1
    print(f"{cont}{' '*8}{xk[-2]:.3e}{' '*3}{derivadas(xk[-2])[0]:.3e}{' '*3}{derivadas(xk[-2])[1]:.3e}{' '*7}{xk[-1]:.3e}{' '*7}{cp:.3e}")
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
plt.title(f'{analise1[analise - 1]}\n NR ')
plt.grid()
plt.savefig(f"{analise1[analise - 1]} NR.png",dpi=300)
plt.show()
