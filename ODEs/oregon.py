import numpy as np


# Oregonator Problem

def f(t, u):
    s = 77.27
    w = 0.161
    q = 8.375e-5
    deriv = [s*(u[1] + u[0]*(1. - q*u[0] - u[1])), \
             (u[2] - u[1]*(1. + u[0]))/s, \
             w*(u[0] - u[2]) ]
    return np.array(deriv)

t = 0
dt = 0.1
tn = 400.
u = [400.0, 1.0, 400.0]
#u = [1., 2., 3.]
u = np.array(u)
