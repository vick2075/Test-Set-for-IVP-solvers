import numpy as np


# motion of an elastic, inextensible thin beam clamped at one end and subject to a force acting at the free end
neqn = 80
N = 40
NN = 2 * N
NSQ = N * N
NQUATR = NSQ * NSQ


def f(t, u):
    DF = np.empty(NN)
    #TH = np.zeros(150)
    TH = u
    
    U = np.zeros(150)
    V = np.zeros(150)
    W = np.zeros(150)
    ALPHA = np.zeros(150)
    BETA = np.zeros(150)
    STH = np.zeros(150)
    CTH = np.zeros(150)
    
    for i in range(1, N):
        THDIFF = TH[i] - TH[i-1]
        STH[i] = np.sin(THDIFF)
        CTH[i] = np.cos(THDIFF)
    
    if t > np.pi:
        TERM1 = (-3.0 * TH[0] + TH[1]) * NQUATR
        V[0] = TERM1
        for i in range(1, N-1):
            TERM1 = (TH[i-1] - 2.0 * TH[i] + TH[i+1]) * NQUATR
            V[i] = TERM1
        TERM1 = (TH[N-2] - TH[N-1]) * NQUATR
        V[N-1] = TERM1
    else:
        FABS = 1.5 * np.sin(t) * np.sin(t)
        FX = -FABS
        FY = FABS

        TERM1 = (-3.0 * TH[0] + TH[1]) * NQUATR
        TERM2 = NSQ * (FY * np.cos(TH[0]) - FX * np.sin(TH[0]))
        V[0]= TERM1 + TERM2

        for i in range(1, N-1):
            TERM1 = (TH[i-1] - 2.0 * TH[i] + TH[i+1]) * NQUATR
            TERM2 = NSQ * (FY * np.cos(TH[i]) - FX * np.sin(TH[i]))
            V[i] = TERM1 + TERM2

        TERM1=(TH[N-2] - TH[N-1]) * NQUATR
        TERM2 = NSQ * (FY * np.cos(TH[N-1]) - FX * np.sin(TH[N-1]))
        V[N-1] = TERM1 + TERM2
        

    W[0] = STH[1] * V[1]
    for i in range(1, N-1):
        W[i] = -STH[i] * V[i-1] + STH[i+1] * V[i+1]
    W[N-1] = -STH[N-1] * V[N-2]

    for i in range(0, N):
        W[i] = W[i] + TH[N+i] * TH[N+i]

    ALPHA[0] = 1.0
    for i in range(1, N):
        ALPHA[i] = 2.0
        BETA[i-1] = -CTH[i]
    ALPHA[N-1] = 3.0
    for i in range(N-2, 0-1, -1):
        Q = BETA[i]/ALPHA[i+1]
        W[i] = W[i] - W[i+1] * Q
        ALPHA[i] = ALPHA[i] - BETA[i] * Q
    W[0] = W[0]/ALPHA[0]
    for i in range(1, N):
        W[i] = (W[i] - BETA[i-1] * W[i-1]) / ALPHA[i]

    U[0] = V[0] - CTH[1] * V[1] + STH[1] * W[1]
    for i in range(1, N-1):
        U[i] = 2.0 * V[i] - CTH[i] * V[i-1] - CTH[i+1] * V[i+1] \
               - STH[i] * W[i-1] + STH[i+1] * W[i+1]
    U[N-1] = 3.0 * V[N-1] - CTH[N-1] * V[N-2] - STH[N-1] * W[N-2]

    for i in range(N):
        DF[i] = TH[N+i]
        DF[N+i] = U[i]

    return DF

t = 0
dt = 1./N
tn = 5.

u = np.empty(neqn)
for i in range(neqn):
    u[i] = 0.

