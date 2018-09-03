import numpy as np
import matplotlib.pyplot as plt
import math

N = 16

def poly_interp(x, coeffs):
    temp = 0
    for i in range(len(coeffs)):
        temp += coeffs[i]*(x**(N - i))

    return temp

x = np.linspace(0, 6, 1000)
f = np.zeros_like(x)
for i in range(len(x)):
    f[i] = x[i] - math.cos(1.5*x[i]) + math.sin(0.4*x[i])

cheb_x = np.zeros((N+1,))
for i in range(N+1):
    cheb_x[i] = (math.cos(i*math.pi/N) + 1)*3.0

samples = np.zeros_like(cheb_x)
for i in range(N+1):
    samples[i] = cheb_x[i] - math.cos(1.5*cheb_x[i]) + math.sin(0.4*cheb_x[i])

poly_params = np.polyfit(cheb_x, samples, N)

f_interp = np.zeros_like(f)
for i in range(len(x)):
    f_interp[i] = poly_interp(x[i], poly_params)

plt.plot(x, f, 'b', x, f_interp, 'g')
for i in cheb_x:
    plt.axvline(x=i, color = 'y', linestyle='dashed')
plt.show()

norm = np.sum(f - f_interp)/len(f)
print(norm)
