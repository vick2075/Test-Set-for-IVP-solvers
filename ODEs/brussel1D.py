import numpy as np


# Brusselator 1D Problem
N= 500; N2=2*N;USDELQ=(N+1)**2; GAMMA=0.02*USDELQ
GAMMA2 = 2.*GAMMA


def f(t, u):
    deriv = np.zeros(1000)
    i=1
    IU=2*i-1
    IV=2*i
    UI=u[IU-1]
    VI=u[IV-1]
    UIM=1.
    VIM=3.
    UIP=u[IU+2-1]
    VIP=u[IV+2-1]
    PROD=UI*UI*VI
    deriv[IU-1]=1.+PROD-4.*UI+GAMMA*(UIM-2.*UI+UIP)
    deriv[IV-1]=3.*UI-PROD+GAMMA*(VIM-2.*VI+VIP)
    for i in range(2,N):
        IU=2*i-1
        IV=2*i
        UI=u[IU-1]
        VI=u[IV-1]
        UIM=u[IU-2-1]
        VIM=u[IV-2-1]
        UIP=u[IU+2-1]
        VIP=u[IV+2-1]
        PROD=UI*UI*VI
        deriv[IU-1]=1.+PROD-4.*UI+GAMMA*(UIM-2.*UI+UIP)
        deriv[IV-1]=3.*UI-PROD+GAMMA*(VIM-2.*VI+VIP)

    i=N
    IU=2*i-1
    IV=2*i
    UI=u[IU-1]
    VI=u[IV-1]
    UIM=u[IU-2-1]
    VIM=u[IV-2-1]
    UIP=1.
    VIP=3.
    PROD=UI*UI*VI
    deriv[IU-1]=1.+PROD-4.*UI+GAMMA*(UIM-2.*UI+UIP)
    deriv[IV-1]=3.*UI-PROD+GAMMA*(VIM-2.*VI+VIP)
            
    return deriv

t = 0
dt = 0.1
tn = 10.

u=np.zeros(1000)
for i in range(1, N+1):
    ANP1 = (N+1)
    XI = i/ANP1
    u[2*i-1] = 3.
    u[2*i-2] = 1.+0.5*np.sin(2.*np.pi*XI)


