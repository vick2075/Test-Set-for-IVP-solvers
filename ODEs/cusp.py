import numpy as np


# CUSP Problem
NNERV = 32; ANERV = NNERV*NNERV; DIFFUS = 1.*ANERV/144.


def f(t, u):
    deriv = np.zeros(96)
    for i in range(1, NNERV+1):
        X = u[3*i-3]
        A = u[3*i-2]
        B = u[3*i-1]

        if i == 1:
            XRIGHT = u[3*NNERV-3]
            ARIGHT = u[3*NNERV-2]
            BRIGHT = u[3*NNERV-1]
        else:
            XRIGHT = u[3*i-6]
            ARIGHT = u[3*i-5]
            BRIGHT = u[3*i-4]
        if i == NNERV:
            XLEFT = u[0]
            ALEFT = u[1]
            BLEFT = u[2]
        else:
            XLEFT = u[3*i]
            ALEFT = u[3*i+1]
            BLEFT = u[3*i+2]
        XDOT = -10000.*(B+X*(A+X*X))
        UU = (X-0.7)*(X-1.3)
        V = UU/(UU+0.1)
        ADOT = B+0.07*V
        BDOT = (1.*(1.-A*A)*B-A)-0.4*X+0.035*V
        deriv[3*i-3]=XDOT+DIFFUS*(XLEFT-2.*X+XRIGHT)
        deriv[3*i-2]=ADOT+DIFFUS*(ALEFT-2.*A+ARIGHT)
        deriv[3*i-1]  =BDOT+DIFFUS*(BLEFT-2.*B+BRIGHT)
            
    return deriv

t = 0
dt = 0.1
tn = 1.1

u=np.zeros(96)
DEL = 2.*np.pi/32
for i in range(1, NNERV+1):
    u[3*i-3] = 0.
    u[3*i-2] = -2.*np.cos(i*DEL)
    u[3*i-1] = 2.*np.sin(i*DEL)


