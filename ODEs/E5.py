import numpy as np


# E5 Problem

def f(t, u):
    A = 7.89e-10
    B = 1.1e+7
    C = 1.13e+3
    CM = 1.13e+9
    deriv = [-A*u[0] - B*u[0]*u[2], \
             A*u[0] - CM*u[1]*u[2], \
             A*u[0] - B*u[0]*u[2] - CM*u[1]*u[2] + C*u[3], \
             B*u[0]*u[2] - C*u[3]]
    return np.array(deriv)

t = 0
dt = 1.0
tn = 1.e+13 #1000.# 
u = [1.76e-3, 0.0, 0.0, 0.0]
u = np.array(u)
