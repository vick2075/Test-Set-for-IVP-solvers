import numpy as np


# Chemical Akzo Nobel Problem
ks = 115.83
k1 = 18.7; k2 = 0.58; k3 = 0.09; k4 = 0.42
kbig = 34.4; kla = 3.3; po2 = 0.9; hen = 737.


def f(t, u):

    if u[1] < 0.0:
        return

    r1 = k1*(u[0]**4)*np.sqrt(u[1])
    r2 = k2*u[2]*u[3]
    r3 = k2/kbig*u[0]*u[4]
    r4 = k3*u[0]*(u[3]**2)
    r5 = k4*(u[5]**2)*np.sqrt(u[1])
    fin = kla*(po2/hen-u[1])

    deriv = [-2.*r1 +r2 -r3 -r4,
             -0.5*r1-r4-0.5*r5 + fin,
             r1 -r2 +r3,
             -r2 +r3 -2.*r4,
             r2 -r3 +r5,
             ks*u[0]*u[3]-u[5]]
    return np.array(deriv)

t = 0
dt = 0.01
tn = 180.

u = np.zeros(6)
u[0] = 0.444
u[1] = 0.00123
u[2] = 0.0
u[3] = 0.007
u[4] = 0.0
u[5] = ks*u[0]*u[3]

M = np.array([[1,0,0,0,0,0],
              [0,1,0,0,0,0],
              [0,0,1,0,0,0],
              [0,0,0,1,0,0],
              [0,0,0,0,1,0],
              [0,0,0,0,0,0]])

var_index = np.array([0,0,0,0,0,0])

