import numpy as np


# Medical Akzo Nobel Problem
        

def f(t, u):
    k = 100.; c = 4.; N = int(400/2); dzeta = 1./N
    dzeta2 = dzeta*dzeta
    dum = (dzeta-1.0)*(dzeta-1.0)/c
    alpha  = 2.0*(dzeta-1.0)*dum/c
    beta = dum*dum

    if t <= 5.:
        phi = 2.
    else:
        phi = 0.

    deriv = np.zeros(400)

    deriv[0] = (phi-2.0*u[0]+u[2])*beta/dzeta2 + \
               alpha*(u[2]-phi)/(2.0*dzeta)-k*u[0]*u[1]
    deriv[1] = -k*u[0]*u[1]

    for j in range(2,N):
        i = 2*j-1
        zeta = j*dzeta
        dum = (zeta-1.0)*(zeta-1.0)/c
        alpha = 2.0*(zeta-1.0)*dum/c
        beta  = dum*dum
        deriv[i-1] = (u[i-2-1]-2.0*u[i-1]+u[i+2-1])*beta/dzeta2+ \
                     alpha*(u[i+2-1]-u[i-2-1])/(2.0*dzeta)-k*u[i-1]*u[i+1-1]
        i = 2*j
        deriv[i-1] = -k*u[i-1]*u[i-1-1]

    deriv[2*N-1-1] = -k*u[2*N-1-1]*u[2*N-1]
    deriv[2*N-1] = -k*u[2*N-1-1]*u[2*N-1]
    
    return deriv

t = 0
dt = 0.1
tn = 20.

u=np.zeros(400)
for i in range(1, 200+1):
    u[2*i-1-1] = 0.
    u[2*i-1] = 1.


