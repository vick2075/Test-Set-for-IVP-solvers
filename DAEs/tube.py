import numpy as np


# Water Tube System

nu = 1.31e-6; g = 9.8; rho = 1.0e+3
rcrit = 2.3e+3; length= 1.0e+3; k = 2.0e-4
d = 1.0; b = 2.0e+2; a = np.pi*d**2/4.


def f(t, u):
    
    nnodes = 13; mu = nu*rho
    phi = np.zeros((nnodes,nnodes)); lamda = np.zeros((nnodes,nnodes))
    p = np.zeros(nnodes); netflo = np.zeros(nnodes); ein = np.zeros(nnodes)
    eout = np.zeros(nnodes); rghres = np.zeros((nnodes,nnodes)); fdba = np.zeros((nnodes,nnodes))

    for j in range(nnodes):
        ein[j] = 0.
        eout[j] = 0.

    that=t/3600.
    that2=that*that

    ein[0] = (1.-np.cos(np.exp(-that)-1.))/200.
    ein[12] = (1.-np.cos(np.exp(-that)-1.))/80.
    eout[9] = that2*(3.*that2-92.*that+720.)/1e+6

    for j in range(nnodes):
        for i in range(nnodes):
            phi[i,j] = 0.
            lamda[i,j] = 1.

    phi[0, 1] = u[0]
    phi[1, 2] = u[1]
    phi[1, 5] = u[2]
    phi[2, 3] = u[3]
    phi[2, 4] = u[4]
    phi[3, 4] = u[5]
    phi[4,9] = u[6]
    phi[5, 4] = u[7]
    phi[6, 3] = u[8]
    phi[6, 7] = u[9]
    phi[7, 4] = u[10]
    phi[7,9] = u[11]
    phi[8, 7] = u[12]
    phi[10, 8] = u[13]
    phi[10,11] = u[14]
    phi[11, 6] = u[15]
    phi[11, 7] = u[16]
    phi[12,10] = u[17]

    lamda[0, 1] = u[18]
    lamda[1, 2] = u[19]
    lamda[1, 5] = u[20]
    lamda[2, 3] = u[21]
    lamda[2, 4] = u[22]
    lamda[3, 4] = u[23]
    lamda[4, 9] = u[24]
    lamda[5, 4] = u[25]
    lamda[6, 3] = u[26]
    lamda[6, 7] = u[27]
    lamda[7, 4] = u[28]
    lamda[7, 9] = u[29]
    lamda[8, 7] = u[30]
    lamda[10, 8] = u[31]
    lamda[10, 11] = u[32]
    lamda[11, 6] = u[33]
    lamda[11, 7] = u[34]
    lamda[12, 10] = u[35]

    p[4] = u[36]
    p[7] = u[37]
    p[0] = u[38]
    p[1] = u[39]
    p[2] = u[40]
    p[3] = u[41]
    p[5] = u[42]
    p[6] = u[43]
    p[8] = u[44]
    p[9] = u[45]
    p[10] = u[46]
    p[11] = u[47]
    p[12] = u[48]

    for j in range(nnodes):
        for i in range(nnodes):
            if lamda[i,j] < 0.:
                ierr = -1
                return
            rtla = np.sqrt(lamda[i,j])
            r = abs(phi[i,j]*d/(nu*a))
            if r > rcrit:
                rghres[i,j] = 1/rtla - 1.74 + \
                              2.*np.log10(2.*k/d + 18.7/(r*rtla))
                fdba[i,j] = p[i] - p[j] - \
                            lamda[i,j]*rho*length*phi[i,j]**2/(a**2*d)
            else:
                rghres[i,j] = 1/rtla - 1.74 + \
                              2.*np.log10(2.*k/d + 18.7/(rcrit*rtla))
                fdba[i,j] = p[i] - p[j] - \
                            32.*mu*length*phi[i,j]/(a*d**2)

    for n in range(nnodes):
        netflo[n] = ein[n]-eout[n]
        for i in range(nnodes):
            netflo[n] = netflo[n] + phi[i,n]
        for j in range(nnodes):
            netflo[n] = netflo[n] - phi[n,j]

    deriv = [fdba[0,1],
             fdba[1,2],
             fdba[1,5],
             fdba[2,3],
             fdba[2,4],
             fdba[3,4],
             fdba[4,9],
             fdba[5,4],
             fdba[6,3],
             fdba[6,7],
             fdba[7,4],
             fdba[7,9],
             fdba[8,7],
             fdba[10,8],
             fdba[10,11],
             fdba[11,6],
             fdba[11,7],
             fdba[12,10],

             rghres[0,1],
             rghres[1,2],
             rghres[1,5],
             rghres[2,3],
             rghres[2,4],
             rghres[3,4],
             rghres[4,9],
             rghres[5,4],
             rghres[6,3],
             rghres[6,7],
             rghres[7,4],
             rghres[7,9],
             rghres[8,7],
             rghres[10,8],
             rghres[10,11],
             rghres[11,6],
             rghres[11,7],
             rghres[12,10],

             netflo[4],
             netflo[7],
             netflo[0],
             netflo[1],
             netflo[2],
             netflo[3],
             netflo[5],
             netflo[6],
             netflo[8],
             netflo[9],
             netflo[10],
             netflo[11],
             netflo[12]]
    return np.array(deriv)
    

t = 0
dt = 0.01
tn = 17. * 3600.

u = np.empty(49)
for i in range(49):
    u[i] = 0.
for i in range(18, 36):
    u[i] = 0.47519404529185289807e-1
for i in range(36, 49):
    u[i] = 109800.


c = b/(rho*g)
v = rho*length/a

M = np.zeros((49,49))
for i in range(18):
    M[i,i] = v
for i in range(18, 36):
    M[i,i] = 0.
for i in range(36, 38):
    M[i,i] = c
for i in range(38, 49):
    M[i,i] = 0.


var_index = np.empty(49)
for i in range(38):
    var_index[i] = 1.
for i in range(38, 49):
    var_index[i] = 2.


