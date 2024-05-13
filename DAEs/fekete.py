import numpy as np


# Fekete problem
    

def f(t, u):
    alpha=0.5; nart = 20
    p=np.zeros((150,3))
    q=np.zeros((150,3))
    pp=np.zeros((150,3))
    qp=np.zeros((150,3))
    lam=np.zeros(150)
    mu=np.zeros(150)
    phi=np.zeros(150)
    gpq=np.zeros(150)
    ff = np.zeros((150,150,3))
    
    deriv = np.zeros(160)
    
    for i in range(nart):
        for k in range(3):
            p[i,k] = u[3*i+k]
            q[i,k] = u[3*nart+3*i+k]
        lam[i] = u[6*nart+i]
        mu[i] = u[7*nart+i]

    for i in range(nart):
        for j in range(nart):
            if i == j:
                for k in range(3):
                    ff[i,j,k] = 0.
            else:
                rn = 0.
                for k in range(3):
                    rn += (p[i,k]-p[j,k])**2
                for k in range(3):
                    ff[i,j,k] = (p[i,k]-p[j,k])/rn

    for i in range(nart):
        for k in range(3):
            pp[i,k] = q[i,k]+2.*mu[i]*p[i,k]
            qp[i,k] = -alpha*q[i,k]+2.*lam[i]*p[i,k]
            for j in range(nart):
                qp[i,k] = qp[i,k]+ff[i,j,k]

    for i in range(nart):
        phi[i] = -1.
        gpq[i] = 0.
        for k in range(3):
            phi[i] = phi[i]+p[i,k]**2
            gpq[i] = gpq[i]+2.*p[i,k]*q[i,k]
        
    for i in range(nart):
        for k in range(3):
            deriv[3*i+k] = pp[i,k]
            deriv[3*nart+3*i+k] = qp[i,k]
        deriv[6*nart+i] = phi[i]
        deriv[7*nart+i] = gpq[i]
        
    return deriv
    

t = 0
dt = 0.001
tn = 1000.

u = np.zeros(160)

for i in range(3):
    alpha=2.*np.pi*(i+1)/3.+np.pi/13.
    beta=3.*np.pi/8.
    u[3*i+0] = np.cos(alpha)*np.cos(beta)
    u[3*i+1] = np.sin(alpha)*np.cos(beta)
    u[3*i+2] = np.sin(beta)

for i in range(3,10):
    alpha=2.*np.pi*(i+1-3)/7.+np.pi/29.
    beta=np.pi/8.
    u[3*i+0]=np.cos(alpha)*np.cos(beta)
    u[3*i+1]=np.sin(alpha)*np.cos(beta)
    u[3*i+2]=np.sin(beta)

for i in range(10,16):
    alpha=2.*np.pi*(i+1-10)/6.+np.pi/7.
    beta=-2.*np.pi/15.
    u[3*i+0]=np.cos(alpha)*np.cos(beta)
    u[3*i+1]=np.sin(alpha)*np.cos(beta)
    u[3*i+2]=np.sin(beta)

for i in range(16,20):
    alpha=2.*np.pi*(i+1-17)/4.+np.pi/17.
    beta=-3.*np.pi/10.
    u[3*i+0]=np.cos(alpha)*np.cos(beta)
    u[3*i+1]=np.sin(alpha)*np.cos(beta)
    u[3*i+2]=np.sin(beta)

for i in range(60,120):
    u[i] = 0.

for i in range(120,160):
    u[i] = 0.


M = np.zeros((160,160))
for i in range(120):
    M[i,i] = 1.
for i in range(120,160):
    M[i,i] = 0.

var_index = np.empty(160)
for i in range(120):
    var_index[i] = 1
for i in range(120,160):
    var_index[i] = 2


