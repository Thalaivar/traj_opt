import numpy as np
import matplotlib.pyplot as plt
import math

N = 16

# cheb points
cheb_x = np.zeros((N+1,))
for i in range(N+1):
    cheb_x[i] = math.cos(i*math.pi/N)

cheb_t = 3.00 - 3.00*cheb_x

# values at cheb points
v = np.zeros((N+1,1))
for i in range(N+1):
    v[i][0] = cheb_t[i] - math.cos(1.5*cheb_t[i]) + math.sin(0.4*cheb_t[i])

# setup differentiation matrix, D
D = np.zeros((N+1, N+1))

for i in range(N+1):
    for j in range(N+1):
        if i == j:
            if i == 0 or i == N:
                D[0][0] = (2*(N**2)+ 1)/(6)
                D[N][N] = -1*D[0][0]
            else:
                D[i][j] = -1*cheb_x[j]/(2*(1 - (cheb_x[j])**2))

        else:
            if i == 0 or i == N: ci = 2
            else: ci = 1
            if j == 0 or j == N: cj = 2
            else: cj = 1

            D[i][j] = (ci/cj)*((-1)**(i+j))/(cheb_x[i] - cheb_x[j])

# discretised ode output
w = np.dot((-1/3)*D, v)
w = w.reshape((17,))

# actual ode output
xdot = np.zeros_like(cheb_x)
for i in range(N+1):
    xdot[i] = 1 + 1.5*math.sin(1.5*cheb_t[i]) + 0.4*math.cos(0.4*cheb_t[i])

err = xdot - w
print(cheb_x)
