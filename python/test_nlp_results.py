import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode


def model(t, X, p):
    x, y = X
    xdot = p[0]*x - p[1]*x*y
    ydot = p[2]*x*y - p[3]*y
    return [xdot, ydot]

def main(params):
    x0 = [0.9, 0.9]
    t0 = 0
    t1 = 100
    dt = 0.01

    x = np.zeros((int((t1 - t0)/dt), 2))
    t = np.zeros((int((t1 - t0)/dt), 1))
    r = ode(model).set_integrator('dopri5').set_initial_value(x0, t0).set_f_params(params)
    i = 0
    while r.successful() and r.t < t1:
        r.integrate(r.t + dt)
        x[i] = r.y
        t[i] = r.t
        i += 1

    plt.plot(t, x[:,0], 'b', t, x[:,1], 'r')
    plt.show()

if __name__ == '__main__':
    main([2/3, 4/3, 1, 1])
