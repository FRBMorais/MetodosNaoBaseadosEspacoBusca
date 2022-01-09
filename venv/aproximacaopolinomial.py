# Aproximação Polinomial
#
import matplotlib.pyplot as plt
import math
import numpy as np

# Parâmetros de entrada
#
#


xmin = 0
xmax = 600
gp = 5 # grau do polinomio
analise = 0  # Função a ser analisada

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
        e = 50_000_000  # [KN/cm^2]
        i = 30_000  # [cm^4]
        le = 600  # [cm]
        w = 2_500  # [KN/cm]

        valor = (w / (120 * e * i * le)) * (- x**5 + 2 * le**2 * x**3 - le**4 * x)  # Função Default
    return valor


# Particionamento do intervalo para a geração dos pontos
x = np.linspace(xmin, xmax, num=gp + 1)  # Vetor contendo os valores em x do particionamento
yy = []

for i in range(1, gp + 1):
    yy.append(funcaox(x[i - 1]))

# Gerando as matrizes x
dim = (gp + 1, gp + 1)
b = np.ones(dim)  # gera matriz de 1s a ser preenchida pelos xs

for i in range(len(x)):
    cont = gp
    for j in range(len(x)):
        b[i][j] = x[i] ** (cont)
        cont -= 1

dim = (gp + 1, 1)
y = np.zeros(dim)  # gera matriz de 1s a ser preenchida pelos xs

for i in range(len(yy)):
    y[i][0] = yy[i]

binv = np.linalg.inv(b)  # inversa de b
a = np.dot(binv, y)  # Matriz dos coeficientes do polinômio
print(a)

# Polinômio gerado com o método de aproximação polinomial
def funfun(xi):
    pol = 0
    cont = gp
    for i in range(gp + 1):
        pol += a[i][0] * xi ** cont
        cont -= 1

    return pol


vetorx = np.linspace(xmin, xmax, num=100)
vetoryestimado = [funfun(vetorx[i]) for i in range(len(vetorx))]
vetoryreal = [funcaox(vetorx[i]) for i in range(len(vetorx))]

plt.plot(vetorx, vetoryreal, vetorx ,vetoryestimado)
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Grau do polinômio {gp}')
plt.grid()
plt.savefig(f"AP_grau = {gp}{analise}.png",dpi=300)
plt.show()

