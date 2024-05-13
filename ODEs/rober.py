import numpy as np


# Robertson

def f(t, u):
    deriv = [-0.04 * u[0] + 10000. * u[1] *u[2],\
             0.04 * u[0] - 10000. * u[1] * u[2] - 3.0e7 * u[1]**2, \
             3.0e7 * u[1]**2 ] 
    return np.array(deriv)

t = 0

dt = 1.0

tn = 1.0e+5 # 1.0e+11 # 1050. #

u = [1.0, 0.0, 0.0]

u = np.array(u)
