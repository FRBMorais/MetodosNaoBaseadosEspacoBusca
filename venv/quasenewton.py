# Quase-Newton
import matplotlib.pyplot as plt
import math
import numpy as np

# Parâmetros de entrada
x1 = 0.2
tol = 10 ** (-8)
h = 10 ** (-4)
xp = 0.1
analise = 2
analise1 = ['Sphere Function', 'Beale Function','Three Ramp Camel', 'Rastrigin', 'Default Problem' ]

print(f"Método Quase Newton\n{30 * '-'}\nX_inicial {x1}\nTolerância: {tol:3e}\nh: {h:3e}\nX_p: {xp}\n{30*'-'}")

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
    derivada1 = (funcaox(x + h) - funcaox(x - h)) / (2 * h)
    # derivada1 = (w / (120 * e * i * le)) * (-5 * x ** 4 + 6 * le ** 2 * x ** 2 - le ** 4)
    return derivada1


cp = 1
xk = [x1]
fun = [funcaox(x1)]
cont = 0
print(f'''Iteração xk{' '*15} f'(xk){' '*10} x(k+1){' '*12} cp''')
while cp > tol:
    xk.append(xk[-1] - ((derivadas(xk[-1])) / ((derivadas(xk[-1]) - derivadas(xp)) / (xk[-1] - xp))))
    cp = abs(xk[-1] - xk[-2])
    cont += 1
    print(f"{cont}{' '*8}{xk[-2]:.7e}{' ' * 3}{derivadas(xk[-1]):.3e}{' '*3}{xk[-1]:.7e}{' '*10}{cp:.3e}")
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
plt.title(f'{analise1[analise - 1]}\n Quase newton ')
plt.grid()
plt.savefig(f"{analise1[analise - 1]} Quase Newton.png",dpi=300)
plt.show()