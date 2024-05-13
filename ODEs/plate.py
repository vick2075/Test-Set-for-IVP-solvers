import numpy as np

# Plate Problem
MX = 8; MY = 5; MACHS1 = 2; MACHS2 = 4
NX=MX;NY=MY;NACHS1=MACHS1;NACHS2=MACHS2
NXM1 = NX-1;NYM1 = NY-1;NDEMI=NX*NY
OMEGA=1000.;STIFFN=100.;WEIGHT=200.
DENOM=NX+1;DELX=2./DENOM;USH4=1./(DELX**4)
FAC=STIFFN*USH4


def f(t, u):
    deriv = np.zeros(80)
    for i in range(1, NX+1):
        for j in range(1, NY+1):
            K = i+NX*(j-1)
            deriv[K-1] = u[K+NDEMI-1]
            UC = 16. *u[K-1]
            if i > 1:
                UC += u[K-1] - 8.*u[K-1-1]
            if i < NX:
                UC += u[K-1] - 8. *u[K+1-1]
            if j > 1:
                UC += u[K-1] - 8. * u[K-NX-1]
            if j < NY:
                UC += u[K-1] - 8. * u[K+NX-1]
            if i > 1 and j > 1:
                UC += 2. * u[K-NX-1-1]
            if i < NX and j > 1:
                UC += 2. * u[K-NX+1-1]
            if i > 1 and j < NY:
                UC += 2. * u[K+NX-1-1]
            if i < NX and j < NY:
                UC += 2. * u[K+NX+1-1]
            if i > 2:
                UC += u[K-2-1]
            if i < NXM1:
                UC += u[K+2-1]
            if j > 2:
                UC += u[K-2*NX-1]
            if j < NYM1:
                UC += u[K+2*NX-1]
            if j == NACHS1 or j == NACHS2:
                XI = i*DELX
                FORCE = np.exp(-5.*(t-XI-2.)**2)+np.exp(-5.*(t-XI-5.)**2)
            else:
                FORCE = 0.
            deriv[K+NDEMI-1] = -OMEGA*u[K+NDEMI-1]-FAC*UC+FORCE*WEIGHT
            
            
    return deriv

t = 0
dt = 0.1
tn = 7.

u=np.zeros(80)


