import numpy as np
from dprod import ddot

# Slider Crank problem

def PDOT(N, DX, INCX, DY, INCY):
    DTEMP = 0.
    if N <= 0:
        return 0.
    if INCX == 1 and INCY == 1:
        M = N % 5
        if M != 0:
            for i in range(M):
                DTEMP = DTEMP + DX[i]*DY[i]
            if N < 5:
                return DTEMP
        
        for i in range(M, N, 5):
            DTEMP = DTEMP + DX[i]*DY[i] + DX[i + 1]*DY[i + 1] + \
                    DX[i + 2]*DY[i + 2] + DX[i + 3]*DY[i + 3] + DX[i + 4]*DY[i + 4]
    else:
        IX = 0 
        IY = 0 
        if INCX < 0:
            IX = (-N+1) * INCX
        if INCY < 0:
            IY = (-N+1) * INCY
        for i in range(N):
            DTEMP = DTEMP + DX[IX]*DY[IY]
            IX = IX + INCX
            IY = IY + INCY
    return DTEMP

def resmbs(ITYP,IEQUA,ICALL,T,X,XD,DELTA,IRES,RPAR,IPAR):
    M1 = 0.36; M2 = 0.151104; M3 = 0.075552
    L1 = 0.15; L2 = 0.30
    J1 = 0.002727; J2 = 0.0045339259
    EE = 0.20e+12; NUE= 0.30; BB = 0.0080; HH = 0.0080
    RHO= 7870.0; GRAV= 0.0; OMEGA = 150.

    Q = np.zeros(20); QD = np.zeros(20); MQ = np.zeros((20,20))
    KQ = np.zeros((20,20)); BQ = np.zeros((20,20)); DQ = np.zeros((20,20))
    c1 = np.zeros(20); c2 = np.zeros(20); c12 = np.zeros(20); c21 = np.zeros(20)
    MQQ = np.zeros(20); KQQ = np.zeros(20); DQQD = np.zeros(20)
    QTBQ = np.zeros(20); BQQD = np.zeros(20); V = np.zeros(2); ALC = np.zeros(3)
    PLC = np.zeros(3); VLC = np.zeros(3); AM = np.zeros((23,23)); GP = np.zeros((3,23))
    F = np.zeros(23)

    IRES = 0; NQ = 4; NP = 7; NL = 3; NX = 3*NP+NL; KU = 4; KV = 0

    FIRST = True

    if FIRST:
        FACM = RHO*BB*HH*L2
        FACK = EE*BB*HH/L2
        FACB = BB*HH*L2

        for i in range(NQ):
            for j in range(NQ):
                MQ[j,i] = 0.
                KQ[j,i] = 0.
                BQ[j,i] = 0.
                DQ[j,i] = 0.
            c1[i] = 0.
            c2[i] = 0.
            c12[i] = 0.
            c21[i] = 0.

        MQ[0,0] = FACM*.5
        MQ[1,1] = FACM*.5
        MQ[2,2] = FACM*8.
        MQ[2,3] = FACM*1.
        MQ[3,2] = FACM*1.
        MQ[3,3] = FACM*2.

        KQ[0,0] = FACK*np.pi**4/24.*(HH/L2)**2
        KQ[1,1] = FACK*np.pi**4*2./3.*(HH/L2)**2
        KQ[2,2] = FACK*16./3.
        KQ[2,3] = -FACK*8./3.
        KQ[3,2] = -FACK*8./3.
        KQ[3,3] = FACK*7./3.

        BQ[0,2] = -FACB*16./np.pi**3
        BQ[0,3] = FACB*(8./np.pi**3-1./np.pi)
        BQ[1,3] = FACB*0.5/np.pi
        BQ[2,0] = FACB*16./np.pi**3
        BQ[3,0] = -FACB*(8./np.pi**3-1./np.pi)
        BQ[3,1] = -FACB*0.5/np.pi

        c1[2]  = FACB*2./3.
        c1[3]  = FACB*1./6.
        c2[0]  = FACB*2./np.pi
        c12[2] = L2*FACB*1./3.
        c12[3] = L2*FACB*1./6.
        c21[0] = L2*FACB*1./np.pi
        c21[1] = -L2*FACB*0.5/np.pi

        if IPAR[1] == 1:
            DQ[0,0] = 5.
            DQ[1,1] = 25.
            DQ[3,2] = 0.5*2.308375455264791e+02
            DQ[2,3] = -0.5*2.62688487992052e+02
            DQ[3,2] = -0.5*2.626884879920526e+02
            DQ[3,3] = 0.5*4.217421837156818e+02

        FIRST = False

    COSP1  = np.cos(X[0])
    COSP2  = np.cos(X[1])
    SINP1  = np.sin(X[0])
    SINP2  = np.sin(X[1])
    COSP12 = np.cos(X[0]-X[1])
    SINP12 = np.sin(X[0]-X[1])
    V[0]   = X[NP]
    V[1]   = X[NP+1]

    for i in range(NQ):
        Q[i] = X[3+i]
        QD[i] = X[NP+3+i]

    c1TQ  = PDOT(NQ,c1,1,Q,1)
    c1TQD = PDOT(NQ,c1,1,QD,1)
    c2TQ  = PDOT(NQ,c2,1,Q,1)
    c2TQD = PDOT(NQ,c2,1,QD,1)
    c12TQ = PDOT(NQ,c12,1,Q,1)
    c12TQD= PDOT(NQ,c12,1,QD,1)
    #'''
    for i in range(NQ):
        MQQ[i] = np.dot(MQ[i],Q)
        KQQ[i] = np.dot(KQ[i],Q)
        DQQD[i] = np.dot(DQ[i],QD)
        QTBQ[i] = np.dot(Q, BQ[i])
        BQQD[i] = np.dot(BQ[i],QD)
    #'''
    '''
    ss1 = PDOT(NQ,MQ,1,Q,1)
    ss2 = PDOT(NQ,KQ,1,Q,1)
    ss3 = PDOT(NQ,DQ,1,QD,1)
    ss4 = PDOT(NQ,Q,1,BQ,1)
    ss5 = PDOT(NQ,BQ,1,QD,1)
    '''
    '''
    for i in range(NQ):
        MQQ[i] = PDOT(NQ,MQ[0,i],1,Q,1)
        KQQ[i] = PDOT(NQ,KQ[0,i],1,Q,1)
        DQQD[i]= PDOT(NQ,DQ[0,i],1,QD,1)
        QTBQ[i]= PDOT(NQ,Q,1,BQ[0,i],1)
        BQQD[i]= PDOT(NQ,BQ[i,0],NQMAX,QD,1)
    '''
    QTMQQ   = PDOT(NQ,Q,1,MQQ,1)
    QDTMQQ  = PDOT(NQ,QD,1,MQQ,1)
    QDTBQQD = PDOT(NQ,QD,1,BQQD,1)

    
    for i in range(NP):
        DELTA[i]    = XD[i]    - X[NP+i]
        DELTA[NP+i] = XD[NP+i] - X[2*NP+i]

    AM[0,0] = J1 + M2*L1*L1
    AM[0,1] = .5*L1*L2*M2*COSP12
    AM[1,1] = J2
    AM[0,2] = 0.
    AM[1,2] = 0.
    AM[2,0] = 0.
    AM[2,1] = 0.
    AM[2,2] = M3
    AM[0,1] = AM[0,1] + RHO*L1*(SINP12*c2TQ+COSP12*c1TQ)
    AM[1,1] = AM[1,1] + QTMQQ + 2.0*RHO*c12TQ

    for i in range(NQ):
        AM[0,3+i] = RHO*L1*(-SINP12*c1[i] + COSP12*c2[i])
        AM[1,3+i] = RHO*c21[i] + RHO*QTBQ[i]
        AM[2,3+i] = 0.

    for i in range(NQ):
        for j in range(i+1):
            AM[3+j,3+i] = MQ[j,i]

    for i in range(NP):
        for j in range(NP):
            AM[j,i] = AM[i,j]
    
    if KU == 0:
        QKU = 0.
    else:
        QKU = Q[KU-1]

    if KV == 0:
        QKV = 0.
    else:
        QKV = Q[KV-1]

    GP[0,0] = L1*COSP1
    GP[0,1] = L2*COSP2 + QKU*COSP2 - QKV*SINP2
    GP[0,2] = 0.
    GP[1,0] = L1*SINP1
    GP[1,1] = L2*SINP2 + QKU*SINP2 + QKV*COSP2
    GP[1,2] = 1.
    GP[2,0] = 1.
    GP[2,1] = 0.
    GP[2,2] = 0.
    
    for i in range(NQ):
        GP[0,3+i] = 0.
        GP[1,3+i] = 0.
        GP[2,3+i] = 0.

    if KU != 0:
        GP[0,3+KU-1] = SINP2
        GP[1,3+KU-1] = -COSP2

    if KV != 0:
        GP[0,3+KV-1] = COSP2
        GP[1,3+KV-1] = SINP2

    F[0] = -.5*L1*GRAV*(M1+2.0*M2)*COSP1-.5*L1*L2*M2*V[1]*V[1]*SINP12
    F[1] = -.5*L2*GRAV*M2*COSP2+.5*L1*L2*M2*V[0]*V[0]*SINP12
    F[2] = 0.
    F[0] = F[0]+ RHO*L1*V[1]*V[1]*(-SINP12*c1TQ+COSP12*c2TQ)- 2.0*RHO*L1*V[1]*(COSP12*c1TQD+SINP12*c2TQD)
    F[1] = F[1]+ RHO*L1*V[0]*V[0]*(SINP12*c1TQ-COSP12*c2TQ) - 2.0*RHO*V[1]*c12TQD - 2.0*V[1]*QDTMQQ- RHO*QDTBQQD - RHO*GRAV*(COSP2*c1TQ-SINP2*c2TQ)

    for i in range(NQ):
        F[3+i] = V[1]*V[1]*MQQ[i]+ \
                 RHO*(V[1]*V[1]*c12[i]+ L1*V[0]*V[0]*(COSP12*c1[i]+SINP12*c2[i])+ \
                      2.0*V[1]*BQQD[i] )- RHO*GRAV*(SINP2*c1[i]+COSP2*c2[i])

    for i in range(NQ):
        F[3+i] = F[3+i] - KQQ[i] - DQQD[i]

    if IPAR[0] == 1:
        FACK = 0.5*EE*BB*HH/L2**2*np.pi**2
        FACB = 80./(np.pi**2*9.)
        F[3] = F[3] - FACK*(Q[0]*Q[3]-FACB*Q[1]*(-4*Q[2]+2*Q[3]))
        F[4] = F[4] - FACK*(4*Q[1]*Q[3]-FACB*Q[0]*(-4*Q[2]+2*Q[3]))
        F[5] = F[5] - FACK*4.*FACB*Q[0]*Q[1]
        F[6] = F[6] - FACK*(0.5*Q[0]**2+2*Q[1]**2-2*FACB*Q[0]*Q[1])

    
   
    ss = np.zeros(23)
    AM = np.vstack([AM,ss])
    sxx = ddot(AM, X, 2*NP)[0:NP]
    
    # PDOT(NP,AM[0,i],1,X[2*NP],1)
    for i in range(NP):
        DELTA[2*NP+i] = sxx[i]- \
                        F[i] + GP[0,i]*X[NX-3]+GP[1,i]*X[NX-2]+GP[2,i]*X[NX-1]
    #print(DELTA)
    if KU == 0:
        QDKU = 0.
    else:
        QDKU = QD[KU-1]
    if KV == 0:
        QDKV = 0.
    else:
        QDKV = QD[KV-1]

    ALC[0] = -L1*SINP1*V[0]*V[0] - (L2+QKU)*SINP2*V[1]*V[1] + \
             2.0*V[1]*(COSP2*QDKU-SINP2*QDKV) - COSP2*V[1]*V[1]*QKV
    ALC[1] =  L1*COSP1*V[0]*V[0] + (L2+QKU)*COSP2*V[1]*V[1] +\
             2.0*V[1]*(SINP2*QDKU+COSP2*QDKV) - SINP2*V[1]*V[1]*QKV
    ALC[2] = 0.0

    for i in range(NP):
        ALC[0] = ALC[0] + GP[0,i]*X[2*NP+i]
        ALC[1] = ALC[1] + GP[1,i]*X[2*NP+i]
        ALC[2] = ALC[2] + GP[2,i]*X[2*NP+i]

    PLC[0] = L1*SINP1 + L2*SINP2 + QKU*SINP2 + QKV*COSP2
    PLC[1] = X[2] - L1*COSP1 - L2*COSP2 -QKU*COSP2 + QKV*SINP2
    PLC[2] = X[0] - OMEGA*T

    VLC[0] = 0.0
    VLC[1] = 0.0
    VLC[2] = -OMEGA

    for i in range(NP):
        VLC[0] = VLC[0] + GP[0,i]*X[NP+i]
        VLC[1] = VLC[1] + GP[1,i]*X[NP+i]
        VLC[2] = VLC[2] + GP[2,i]*X[NP+i]

    if IEQUA == 2:
        DELTA[0] = PLC[0]
        DELTA[1] = PLC[1]
        DELTA[2] = PLC[2]
        DELTA[3] = VLC[0]
        DELTA[4] = VLC[1]
        DELTA[5] = VLC[2]
    else:
        if ITYP == 0:
            DELTA[NX-3] = PLC[0]
            DELTA[NX-2] = PLC[1]
            DELTA[NX-1] = PLC[2]
        elif ITYP == 1:
            DELTA[NX-3] = VLC[0]
            DELTA[NX-2] = VLC[1]
            DELTA[NX-1] = VLC[2]
        elif ITYP == 3:
            DELTA[NX-3] = ALC[0]
            DELTA[NX-2] = ALC[1]
            DELTA[NX-1] = ALC[2]
        
    return DELTA, IRES


                

def f(t, u):
    ipar = np.empty(2)
    fi = np.zeros(24)
    ipar[0] = 0; ipar[1] = 0; ityp = 1; iequa = 0

    deriv = []
    dy = []
    for i in range(24):
        dy.insert(i , 0.)

    fi, ires = resmbs(ityp,iequa,0,t,u,np.array(dy),fi,0,0,ipar)
    
    for i in range(24):
        #fi[i] = -fi[i]
        deriv.insert(i, -fi[i])
    
    return np.array(deriv)
    

t = 0
dt = 1.e-3
tn = 0.1

u = np.zeros(24)
# Position variables
# phi1, phi2, x3, q1, q2, q3, q4
u[0]= 0.000000000000000
u[1]= 0.000000000000000
u[2]= 4.500169330000000e-1
u[3]= 0.000000000000000
u[4]= 0.000000000000000
u[5] = 1.033398630000000e-5
u[6] = 1.693279690000000e-5
# Initial values velocity variables
u[7] = 1.500000000000000e+2
u[8] =-7.499576703969453e+1
u[9] =-2.689386719979040e-6
u[10] = 4.448961125815990e-1
u[11] = 4.634339319238670e-3
u[12] =-1.785910760000550e-6
u[13] =-2.689386719979040e-6
# Initial values acceleration variables
u[14] = 0.000000000000000
u[15] =-1.344541576709835e-3
u[16] =-5.062194924490193e+3
u[17] =-6.829725665986310e-5
u[18] = 1.813207639590617e-20
u[19] =-4.268463266810281
u[20] = 2.098339029337557e-1
# Lagrange multipliers
u[21] =-6.552727150584648e-8
u[22] = 3.824589509350831e+2
u[23] =-4.635908708561371e-9


M = np.zeros((24,24))
for i in range(14):
    M[i,i] = 1.
for i in range(14, 24):
    M[i,i] = 0.

var_index = np.empty(24)
for i in range(14):
    var_index[i] = 1.
for i in range(14, 24):
    var_index[i] = 2.


