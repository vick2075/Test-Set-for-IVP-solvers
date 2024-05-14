import numpy as np


# Two Bit Adding Unit problem
        
tt = np.zeros(65)
for i in range(65):
    tt[i] = 5.0 * i

CTIME = 1.e+4; STIFF = 5.
RGS= .4e+2/(CTIME*STIFF)
RGD= .4e+2/(CTIME*STIFF)
RBS= .1e+3/(CTIME*STIFF)
RBD= .1e+3/(CTIME*STIFF)
CGS= .6e-4*CTIME
CGD= .6e-4*CTIME
CBD= 2.4e-5*CTIME
CBS= 2.4e-5*CTIME
DELTA= 0.2e-1
CURIS= 1.e-15*CTIME*STIFF
VTH= 25.85
VDD=5.
VBB=-2.5
CLOAD= 0.0
COUT= 2.e-4*CTIME-CLOAD

def FCN(n, x, y, ff, ierr):

    v1, v1d = PULSE(x, 0.0, 5.0, 0.0, 5.0, 5.0, 5.0, 20.0)
    v2, v2d = PULSE(x, 0.0, 5.0, 10.0, 5.0, 15.0, 5.0, 40.0)
    v3, v3d = PULSE(x, 0.0, 5.0, 30.0, 5.0, 35.0, 5.0, 80.0)
    v4, v4d = PULSE(x, 0.0, 5.0, 70.0, 5.0, 75.0, 5.0, 160.0)
    CIN, CIND = PULSE(x, 0.0, 5.0, 150.0, 5.0, 155.0, 5.0, 320.0)

    ff = NOR(n, 1, v1, v2, v1d, v2d, y, ff, ierr)
    ff = ANDOI(n, 14, y[4], v2, v1, 0.0, v2d, v1d, y, ff, ierr)
    ff = NOR(n, 32, y[17], CIN, 0.0, CIND, y, ff, ierr)
    ff = ANDOI(n, 45, y[35], CIN, y[17], 0.0, CIND, 0.0, y, ff, ierr)
    ff = ANDOI(n, 63, y[4], CIN, y[17], 0.0, CIND, 0.0, y, ff, ierr)
    ff = NOR(n, 81, v3, v4, v3d, v4d, y, ff, ierr)
    ff = ANDOI(n, 94, y[84], v4, v3, 0.0, v4d, v3d, y, ff, ierr)
    ff = NAND(n, 112, y[66], y[97], 0.0, 0.0, y, ff, ierr)
    ff = ORANI(n, 126, y[115], y[66], y[97], 0.0, 0.0, 0.0, y, ff, ierr)
    ff = ANDOIP(n, 144, y[84], y[4], y[97], 0.0, 0.0, 0.0, y, ff, ierr)

    ff[161] = -(y[161] - y[165]) / RGS - IDS(3, y[162] - y[161], y[97] - y[161], y[163] - y[165], y[97] - y[162], y[164] - y[147], ierr)
    ff[162] = -(y[162] - y[147]) / RGD + IDS(3, y[162] - y[161], y[97] - y[161], y[163] - y[165], y[97] - y[162], y[164] - y[147], ierr)
    ff[163] = -(y[163] - VBB) / RBS + IBS(y[163] - y[165])
    ff[164] = -(y[164] - VBB) / RBD + IBD(y[164] - y[147])
    ff[165] = -IBS(y[163] - y[165]) - (y[165] - y[161]) / RGS - IBD(y[169] - y[165]) - (y[165] - y[167]) / RGD
    ff[166] = -(y[166] - y[170]) / RGS - IDS(3, y[167] - y[166], y[17] - y[166], y[168] - y[170], y[17] - y[167], y[169] - y[165], ierr)
    ff[167] = -(y[167] - y[165]) / RGD + IDS(3, y[167] - y[166], y[17] - y[166], y[168] - y[170], y[17] - y[167], y[169] - y[165], ierr)
    ff[168] = -(y[168] - VBB) / RBS + IBS(y[168] - y[170])
    ff[169] = -(y[169] - VBB) / RBD + IBD(y[169] - y[165])
    ff[170] = -IBS(y[168] - y[170]) - (y[170] - y[166]) / RGS - IBD(y[174] - y[170]) - (y[170] - y[172]) / RGD
    ff[171] = CGS * CIND - y[171] / RGS - IDS(3, y[172] - y[171], CIN - y[171], y[173], CIN - y[172], y[174] - y[170], ierr)
    ff[172] = CGD * CIND - (y[172] - y[170]) / RGD + IDS(3, y[172] - y[171], CIN - y[171], y[173], CIN - y[172], y[174] - y[170], ierr)
    ff[173] = -(y[173] - VBB) / RBS + IBS(y[173])
    ff[174] = -(y[174] - VBB) / RBD + IBD(y[174] - y[170])

    if ierr == -1:
        return

    return ff

def PULSE(X, LOW, HIGH, DELAY, T1, T2, T3, PERIOD):
    TIME = X % PERIOD

    if TIME > (DELAY+T1+T2+T3):
        VIN = LOW
        VIND= 0.
    elif TIME > (DELAY+T1+T2):
        VIN = ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW
        VIND= -((HIGH-LOW)/T3)
    elif TIME > (DELAY+T1):
        VIN = HIGH
        VIND= 0.
    elif TIME > DELAY:
        VIN = ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW
        VIND= ((HIGH-LOW)/T1)
    else:
        VIN = LOW
        VIND=0.

    return VIN, VIND

def NOR(N,I,U1,U2,U1D,U2D,Y,FF,ierr):
    
    FF[I-1] = -(Y[I-1]-Y[I+4-1])/RGS -IDS(0,Y[I+1-1]-Y[I-1],Y[I+4-1]-Y[I-1],Y[I+2-1]-Y[I+4-1],Y[I+4-1]-Y[I+1-1],Y[I+3-1]-VDD,ierr)
    FF[I+1-1]= -(Y[I+1-1]-VDD)/RGD+IDS(0,Y[I+1-1]-Y[I-1],Y[I+4-1]-Y[I-1], \
                                  Y[I+2-1]-Y[I+4-1],Y[I+4-1]-Y[I+1-1],Y[I+3-1]-VDD,ierr)
    FF[I+2-1]= -(Y[I+2-1]-VBB)/RBS + IBS(Y[I+2-1]-Y[I+4-1])
    FF[I+3-1]= -(Y[I+3-1]-VBB)/RBD + IBD(Y[I+3-1]-VDD)

    FF[I+4-1] = -(Y[I+4-1] - Y[I-1])/RGS - IBS(Y[I+2-1] - Y[I+4-1]) - (Y[I+4-1] - Y[I+6-1])/RGD - IBD(Y[I+8-1] - Y[I+4-1]) - (Y[I+4-1] - Y[I+10-1])/RGD - IBD(Y[I+12-1] - Y[I+4-1])
    FF[I+5-1] = CGS*U1D - Y[I+5-1]/RGS - IDS(1, Y[I+6-1] - Y[I+5-1], U1 - Y[I+5-1], Y[I+7-1], U1 - Y[I+6-1], Y[I+8-1] - Y[I+4-1], ierr)
    FF[I+6-1] = CGD*U1D - (Y[I+6-1] - Y[I+4-1])/RGD + IDS(1, Y[I+6-1] - Y[I+5-1], U1 - Y[I+5-1], Y[I+7-1], U1 - Y[I+6-1], Y[I+8-1] - Y[I+4-1], ierr)
    FF[I+7-1] = -(Y[I+7-1] - VBB)/RBS + IBS(Y[I+7-1])
    FF[I+8-1] = -(Y[I+8-1] - VBB)/RBD + IBD(Y[I+8-1] - Y[I+4-1])
    FF[I+9-1] = CGS*U2D - Y[I+9-1]/RGS - IDS(1, Y[I+10-1] - Y[I+9-1], U2 - Y[I+9-1], Y[I+11-1], U2 - Y[I+10-1], Y[I+12-1] - Y[I+4-1], ierr)
    FF[I+10-1] = CGD*U2D - (Y[I+10-1] - Y[I+4-1])/RGD + IDS(1, Y[I+10-1] - Y[I+9-1], U2 - Y[I+9-1], Y[I+11-1], U2 - Y[I+10-1], Y[I+12-1] - Y[I+4-1], ierr)
    FF[I+11-1] = -(Y[I+11-1] - VBB)/RBS + IBS(Y[I+11-1])
    FF[I+12-1] = -(Y[I+12-1] - VBB)/RBD + IBD(Y[I+12-1] - Y[I+4-1])

    if ierr == -1:
        return

    return FF

def ANDOIP(N,I,U1,U2,U3,U1D,U2D,U3D,Y,FF,ierr):
    

    FF[I-1] = -(Y[I-1] - Y[I+4-1])/RGS - IDS(0, Y[I+1-1] - Y[I-1], Y[I+4-1] - Y[I-1], Y[I+2-1] - Y[I+4-1], Y[I+4-1] - Y[I+1-1], Y[I+3-1] - VDD, ierr)
    FF[I+1-1] = -(Y[I+1-1] - VDD)/RGD + IDS(0, Y[I+1-1] - Y[I-1], Y[I+4-1] - Y[I-1], Y[I+2-1] - Y[I+4-1], Y[I+4-1] - Y[I+1-1], Y[I+3-1] - VDD, ierr)
    FF[I+2-1] = -(Y[I+2-1] - VBB)/RBS + IBS(Y[I+2-1] - Y[I+4-1])
    FF[I+3-1] = -(Y[I+3-1] - VBB)/RBD + IBD(Y[I+3-1] - VDD)

    FF[I+4-1] = -(Y[I+4-1] - Y[I-1])/RGS - IBS(Y[I+2-1] - Y[I+4-1]) - (Y[I+4-1] - Y[I+6-1])/RGD - IBD(Y[I+8-1] - Y[I+4-1]) - (Y[I+4-1] - Y[I+10-1])/RGD - IBD(Y[I+12-1] - Y[I+4-1]) - (Y[I+4-1] - Y[163-1])/RGD - IBD(Y[165-1] - Y[I+4-1])
    FF[I+5-1] = CGS*U1D - Y[I+5-1]/RGS - IDS(1, Y[I+6-1] - Y[I+5-1], U1 - Y[I+5-1], Y[I+7-1], U1 - Y[I+6-1], Y[I+8-1] - Y[I+4-1], ierr)
    FF[I+6-1] = CGD*U1D - (Y[I+6-1] - Y[I+4-1])/RGD + IDS(1, Y[I+6-1] - Y[I+5-1], U1 - Y[I+5-1], Y[I+7-1], U1 - Y[I+6-1], Y[I+8-1] - Y[I+4-1], ierr)
    FF[I+7-1] = -(Y[I+7-1] - VBB)/RBS + IBS(Y[I+7-1])
    FF[I+8-1] = -(Y[I+8-1] - VBB)/RBD + IBD(Y[I+8-1] - Y[I+4-1])
    FF[I+9-1] = CGS*U2D - (Y[I+9-1] - Y[I+13-1])/RGS - IDS(2, Y[I+10-1] - Y[I+9-1], U2 - Y[I+9-1], Y[I+11-1] - Y[I+13-1], U2 - Y[I+10-1], Y[I+12-1] - Y[I+4-1], ierr)
    FF[I+10-1] = CGD*U2D - (Y[I+10-1] - Y[I+4-1])/RGD + IDS(2, Y[I+10-1] - Y[I+9-1], U2 - Y[I+9-1], Y[I+11-1] - Y[I+13-1], U2 - Y[I+10-1], Y[I+12-1] - Y[I+4-1], ierr)
    FF[I+11-1] = -(Y[I+11-1] - VBB)/RBS + IBS(Y[I+11-1] - Y[I+13-1])
    FF[I+12-1] = -(Y[I+12-1] - VBB)/RBD + IBD(Y[I+12-1] - Y[I+4-1])
    FF[I+13-1] = -(Y[I+13-1] - Y[I+9-1])/RGS - IBS(Y[I+11-1] - Y[I+13-1]) - (Y[I+13-1] - Y[I+15-1])/RGD - IBD(Y[I+17-1] - Y[I+13-1])
    FF[I+14-1] = CGS*U3D - Y[I+14-1]/RGS - IDS(2, Y[I+15-1] - Y[I+14-1], U3 - Y[I+14-1], Y[I+16-1], U3 - Y[I+15-1], Y[I+17-1] - Y[I+13-1], ierr)
    FF[I+15-1] = CGD*U3D - (Y[I+15-1] - Y[I+13-1])/RGD + IDS(2, Y[I+15-1] - Y[I+14-1], U3 - Y[I+14-1], Y[I+16-1], U3 - Y[I+15-1], Y[I+17-1] - Y[I+13-1], ierr)
    FF[I+16-1] = -(Y[I+16-1] - VBB)/RBS + IBS(Y[I+16-1])
    FF[I+17-1] = -(Y[I+17-1] - VBB)/RBD + IBD(Y[I+17-1] - Y[I+13-1])

    if ierr == -1:
        return

    return FF


def ANDOI(N,I,U1,U2,U3,U1D,U2D,U3D,Y,FF,ierr):
    

    FF[I-1] = -(Y[I-1] - Y[I+4-1])/RGS - IDS(0, Y[I+1-1] - Y[I-1], Y[I+4-1] - Y[I-1], Y[I+2-1] - Y[I+4-1], Y[I+4-1] - Y[I+1-1], Y[I+3-1] - VDD, ierr)
    FF[I+1-1] = -(Y[I+1-1] - VDD)/RGD + IDS(0, Y[I+1-1] - Y[I-1], Y[I+4-1] - Y[I-1], Y[I+2-1] - Y[I+4-1], Y[I+4-1] - Y[I+1-1], Y[I+3-1] - VDD, ierr)
    FF[I+2-1] = -(Y[I+2-1] - VBB)/RBS + IBS(Y[I+2-1] - Y[I+4-1])
    FF[I+3-1] = -(Y[I+3-1] - VBB)/RBD + IBD(Y[I+3-1] - VDD)

    FF[I+4-1] = -(Y[I+4-1]-Y[I-1])/RGS - IBS(Y[I+2-1]-Y[I+4-1]) - (Y[I+4-1]-Y[I+6-1])/RGD - IBD(Y[I+8-1]-Y[I+4-1]) - (Y[I+4-1]-Y[I+10-1])/RGD - IBD(Y[I+12-1]-Y[I+4-1])
    FF[I+5-1] = CGS*U1D - Y[I+5-1]/RGS - IDS(1, Y[I+6-1]-Y[I+5-1], U1-Y[I+5-1], Y[I+7-1], U1-Y[I+6-1], Y[I+8-1]-Y[I+4-1], ierr)
    FF[I+6-1] = CGD*U1D - (Y[I+6-1]-Y[I+4-1])/RGD + IDS(1, Y[I+6-1]-Y[I+5-1], U1-Y[I+5-1], Y[I+7-1], U1-Y[I+6-1], Y[I+8-1]-Y[I+4-1], ierr)
    FF[I+7-1] = -(Y[I+7-1]-VBB)/RBS + IBS(Y[I+7-1])
    FF[I+8-1] = -(Y[I+8-1]-VBB)/RBD + IBD(Y[I+8-1]-Y[I+4-1])
    FF[I+9-1] = CGS*U2D - (Y[I+9-1]-Y[I+13-1])/RGS - IDS(2, Y[I+10-1]-Y[I+9-1], U2-Y[I+9-1], Y[I+11-1]-Y[I+13-1], U2-Y[I+10-1], Y[I+12-1]-Y[I+4-1], ierr)
    FF[I+10-1] = CGD*U2D - (Y[I+10-1]-Y[I+4-1])/RGD + IDS(2, Y[I+10-1]-Y[I+9-1], U2-Y[I+9-1], Y[I+11-1]-Y[I+13-1], U2-Y[I+10-1], Y[I+12-1]-Y[I+4-1], ierr)
    FF[I+11-1] = -(Y[I+11-1]-VBB)/RBS + IBS(Y[I+11-1]-Y[I+13-1])
    FF[I+12-1] = -(Y[I+12-1]-VBB)/RBD + IBD(Y[I+12-1]-Y[I+4-1])
    FF[I+13-1] = -(Y[I+13-1]-Y[I+9-1])/RGS - IBS(Y[I+11-1]-Y[I+13-1]) - (Y[I+13-1]-Y[I+15-1])/RGD - IBD(Y[I+17-1]-Y[I+13-1])
    FF[I+14-1] = CGS*U3D - Y[I+14-1]/RGS - IDS(2, Y[I+15-1]-Y[I+14-1], U3-Y[I+14-1], Y[I+16-1], U3-Y[I+15-1], Y[I+17-1]-Y[I+13-1], ierr)
    FF[I+15-1] = CGD * U3D - (Y[I+15-1] - Y[I+13-1]) / RGD + IDS(2, Y[I+15-1] - Y[I+14-1], U3 - Y[I+14-1], Y[I+16-1], U3 - Y[I+15-1], Y[I+17-1] - Y[I+13-1], ierr)
    FF[I+16-1] = -(Y[I+16-1] - VBB) / RBS + IBS(Y[I+16-1])
    FF[I+17-1] = -(Y[I+17-1] - VBB) / RBD + IBD(Y[I+17-1] - Y[I+13-1])

    if ierr == -1:
        return

    return FF

def NAND(N,I,U1,U2,U1D,U2D,Y,FF,ierr):
    

    FF[I-1] = -(Y[I-1] - Y[I+4-1]) / RGS - IDS(0, Y[I+1-1] - Y[I-1], Y[I+4-1] - Y[I-1], Y[I+2-1] - Y[I+4-1], Y[I+4-1] - Y[I+1-1], Y[I+3-1] - VDD, ierr)
    FF[I+1-1] = -(Y[I+1-1] - VDD) / RGD + IDS(0, Y[I+1-1] - Y[I-1], Y[I+4-1] - Y[I-1], Y[I+2-1] - Y[I+4-1], Y[I+4-1] - Y[I+1-1], Y[I+3-1] - VDD, ierr)
    FF[I+2-1] = -(Y[I+2-1] - VBB) / RBS + IBS(Y[I+2-1] - Y[I+4-1])
    FF[I+3-1] = -(Y[I+3-1] - VBB) / RBD + IBD(Y[I+3-1] - VDD)

    FF[I+4-1]= -(Y[I+4-1]-Y[I-1])/RGS-IBS(Y[I+2-1]-Y[I+4-1])-(Y[I+4-1]-Y[I+6-1])/RGD-IBD(Y[I+8-1]-Y[I+4-1])
    FF[I+5-1]=CGS*U1D-(Y[I+5-1]-Y[I+9-1])/RGS-IDS(2,Y[I+6-1]-Y[I+5-1],U1-Y[I+5-1],Y[I+7-1]-Y[I+9-1],U1-Y[I+6-1],Y[I+8-1]-Y[I+4-1],ierr)
    FF[I+6-1]=CGD*U1D-(Y[I+6-1]-Y[I+4-1])/RGD+IDS(2,Y[I+6-1]-Y[I+5-1],U1-Y[I+5-1],Y[I+7-1]-Y[I+9-1],U1-Y[I+6-1],Y[I+8-1]-Y[I+4-1],ierr)
    FF[I+7-1]=-(Y[I+7-1]-VBB)/RBS + IBS(Y[I+7-1]-Y[I+9-1])
    FF[I+8-1]=-(Y[I+8-1]-VBB)/RBD + IBD(Y[I+8-1]-Y[I+4-1])

    FF[I+9-1]=-(Y[I+9-1]-Y[I+5-1])/RGS-IBS(Y[I+7-1]-Y[I+9-1])-(Y[I+9-1]-Y[I+11-1])/RGD-IBD(Y[I+13-1]-Y[I+9-1])
    FF[I+10-1]=CGS*U2D-Y[I+10-1]/RGS-IDS(2,Y[I+11-1]-Y[I+10-1],U2-Y[I+10-1],Y[I+12-1],U2-Y[I+11-1],Y[I+13-1]-Y[I+9-1],ierr)
    FF[I+11-1]=CGD*U2D-(Y[I+11-1]-Y[I+9-1])/RGD+IDS(2,Y[I+11-1]-Y[I+10-1],U2-Y[I+10-1],Y[I+12-1],U2-Y[I+11-1],Y[I+13-1]-Y[I+9-1],ierr)
    FF[I+12-1]=-(Y[I+12-1]-VBB)/RBS + IBS(Y[I+12-1])
    FF[I+13-1]=-(Y[I+13-1]-VBB)/RBD + IBD(Y[I+13-1]-Y[I+9-1])

    if ierr == -1:
        return

    return FF


def ORANI(N,I,U1,U2,U3,U1D,U2D,U3D,Y,FF,ierr):
    

    FF[I-1] = -(Y[I-1] - Y[I + 4-1]) / RGS - IDS(0, Y[I + 1-1] - Y[I-1], Y[I + 4-1] - Y[I-1], Y[I + 2-1] - Y[I + 4-1], Y[I + 4-1] - Y[I + 1-1], Y[I + 3-1] - VDD,ierr)
    FF[I + 1-1] = -(Y[I + 1-1] - VDD) / RGD + IDS(0, Y[I + 1-1] - Y[I-1], Y[I + 4-1] - Y[I-1], Y[I + 2-1] - Y[I + 4-1], Y[I + 4-1] - Y[I + 1-1], Y[I + 3-1] - VDD,ierr)
    FF[I + 2-1] = -(Y[I + 2-1] - VBB) / RBS + IBS(Y[I + 2-1] - Y[I + 4-1])
    FF[I + 3-1] = -(Y[I + 3-1] - VBB) / RBD + IBD(Y[I + 3-1] - VDD)

    FF[I + 4-1] = -(Y[I + 4-1] - Y[I-1]) / RGS - IBS(Y[I + 2-1] - Y[I + 4-1]) - (Y[I + 4-1] - Y[I + 6-1]) / RGD - IBD(Y[I + 8-1] - Y[I + 4-1])
    FF[I + 5-1] = CGS * U1D - (Y[I + 5-1] - Y[I + 9-1]) / RGS - IDS(2, Y[I + 6-1] - Y[I + 5-1], U1 - Y[I + 5-1], Y[I + 7-1] - Y[I + 9-1],U1 - Y[I + 6-1], Y[I + 8-1] - Y[I + 4-1],ierr)
    FF[I + 6-1] = CGD * U1D - (Y[I + 6-1] - Y[I + 4-1]) / RGD + IDS(2, Y[I + 6-1] - Y[I + 5-1], U1 - Y[I + 5-1], Y[I + 7-1] - Y[I + 9-1],U1 - Y[I + 6-1], Y[I + 8-1] - Y[I + 4-1],ierr)
    FF[I + 7-1] = -(Y[I + 7-1] - VBB) / RBS + IBS(Y[I + 7-1] - Y[I + 9-1])
    FF[I + 8-1] = -(Y[I + 8-1] - VBB) / RBD + IBD(Y[I + 8-1] - Y[I + 4-1])

    FF[I + 9-1] = -(Y[I + 9-1] - Y[I + 5-1]) / RGS - IBS(Y[I + 7-1] - Y[I + 9-1]) - (Y[I + 9-1] - Y[I + 11-1]) / RGD - IBD(Y[I + 13-1] - Y[I + 9-1]) - (Y[I + 9-1] - Y[I + 15-1]) / RGD - IBD(Y[I + 17-1] - Y[I + 9-1])
    FF[I + 10-1] = CGS * U2D - Y[I + 10-1] / RGS - IDS(2, Y[I + 11-1] - Y[I + 10-1], U2 - Y[I + 10-1], Y[I + 12-1], U2 - Y[I + 11-1],Y[I + 13-1] - Y[I + 9-1],ierr)
    FF[I + 11-1] = CGD * U2D - (Y[I + 11-1] - Y[I + 9-1]) / RGD + IDS(2, Y[I + 11-1] - Y[I + 10-1], U2 - Y[I + 10-1], Y[I + 12-1], U2 - Y[I + 11-1], Y[I + 13-1] - Y[I + 9-1],ierr)
    FF[I + 12-1] = -(Y[I + 12-1] - VBB) / RBS + IBS(Y[I + 12-1])
    FF[I + 13-1] = -(Y[I + 13-1] - VBB) / RBD + IBD(Y[I + 13-1] - Y[I + 9-1])
    FF[I + 14-1] = CGS * U3D - Y[I + 14-1] / RGS - IDS(2, Y[I + 15-1] - Y[I + 14-1], U3 - Y[I + 14-1], Y[I + 16-1], U3 - Y[I + 15-1],Y[I + 17-1] - Y[I + 9-1],ierr)
    FF[I + 15-1] = CGD * U3D - (Y[I + 15-1] - Y[I + 9-1]) / RGD + IDS(2, Y[I + 15-1] - Y[I + 14-1], U3 - Y[I + 14-1], Y[I + 16-1],U3 - Y[I + 15-1], Y[I + 17-1] - Y[I + 9-1],ierr)
    FF[I + 16-1] = -(Y[I + 16-1] - VBB) / RBS + IBS(Y[I + 16-1])
    FF[I + 17-1] = -(Y[I + 17-1] - VBB) / RBD + IBD(Y[I + 17-1] - Y[I + 9-1])

    if ierr == -1:
        return

    return FF


def GCN(N,U,G):

    for i in range(N):
        G[i] = 0.


    G = DNOR(N,U,1,G)
    G = DANDOI(N,U,14,G)
    G = DNOR(N,U,32,G)
    G = DANDOI(N,U,45,G)
    G = DANDOI(N,U,63,G)
    G = DNOR(N,U,81,G)
    G = DANDOI(N,U,94,G)
    G = DNAND(N,U,112,G)
    G = DORANI(N,U,126,G)
    G = DANDOI(N,U,144,G)

    # Capacitive coupling result node nor-gate 1
    G[5-1] += CGS * (U[5-1] - U[19-1]) + CGD * (U[5-1] - U[20-1]) + CGS * (U[5-1] - U[68-1]) + CGD * (U[5-1] - U[69-1]) + CGS * (U[5-1] - U[153-1]) + CGD * (U[5-1] - U[154-1])
    G[19-1] -= CGS * U[5-1]
    G[20-1] -= CGD * U[5-1]
    G[68-1] -= CGS * U[5-1]
    G[69-1] -= CGD * U[5-1]
    G[153-1] -= CGS * U[5-1]
    G[154-1] -= CGD * U[5-1]

    # Capacitive coupling result node andoi-gate 1
    G[18-1] += CGS * (U[18-1] - U[37-1]) + CGD * (U[18-1] - U[38-1]) + CGS * (U[18-1] - U[59-1]) + CGD * (U[18-1] - U[60-1]) + CGS * (U[18-1] - U[77-1]) + CGD * (U[18-1] - U[78-1]) + CGS * (U[18-1] - U[167-1]) + CGD * (U[18-1] - U[168-1])
    G[37-1] -= CGS * U[18-1]
    G[38-1] -= CGD * U[18-1]
    G[59-1] -= CGS * U[18-1]
    G[60-1] -= CGD * U[18-1]
    G[77-1] -= CGS * U[18-1]
    G[78-1] -= CGD * U[18-1]

    # Capacitive coupling result node nor-gate 2
    G[36-1] += CGS * (U[36-1] - U[50-1]) + CGD * (U[36-1] - U[51-1])
    G[50-1] -= CGS * U[36-1]
    G[51-1] -= CGD * U[36-1]

    # Capacitive coupling result node andoi-gate 2 === s0
    G[49-1] += COUT * U[49-1]

    # Capacitive coupling result node andoi-gate 3
    G[67-1] += CGS * (U[67-1] - U[117-1]) + CGD * (U[67-1] - U[118-1]) + CGS * (U[67-1] - U[136-1]) + CGD * (U[67-1] - U[137-1])
    G[117-1] -= CGS * U[67-1]
    G[118-1] -= CGD * U[67-1]
    G[136-1] -= CGS * U[67-1]
    G[137-1] -= CGD * U[67-1]

    # Capacitive coupling result node nor-gate 3
    G[85-1] += CGS * (U[85-1] - U[99-1]) + CGD * (U[85-1] - U[100-1]) + CGS * (U[85-1] - U[149-1]) + CGD * (U[85-1] - U[150-1])
    G[99-1] -= CGS * U[85-1]
    G[100-1] -= CGD * U[85-1]
    G[149-1] -= CGS * U[85-1]
    G[150-1] -= CGD * U[85-1]

    # Capacitive coupling result node andoi-gate 4
    G[98-1] += CGS * (U[98-1] - U[122-1]) + CGD * (U[98-1] - U[123-1]) + CGS * (U[98-1] - U[140-1]) + CGD * (U[98-1] - U[141-1]) + CGS * (
                U[98-1] - U[158-1]) + CGD * (U[98-1] - U[159-1]) + CGS * (U[98-1] - U[162-1]) + CGD * (U[98-1] - U[163-1])
    G[122-1] -= CGS * U[98-1]
    G[123-1] -= CGD * U[98-1]
    G[140-1] -= CGS * U[98-1]
    G[141-1] -= CGD * U[98-1]
    G[158-1] -= CGS * U[98-1]
    G[159-1] -= CGD * U[98-1]

    # Capacitive coupling result nand-gate
    G[116-1] += CGS * (U[116-1] - U[131-1]) + CGD * (U[116-1] - U[132-1])
    G[131-1] -= CGS * U[116-1]
    G[132-1] -= CGD * U[116-1]

    # Capacitive coupling result node orani-gate === s1
    G[130-1] += COUT * U[130-1]

    # Capacitive coupling result andoi-gate 5 === Cinvers
    G[148-1] += CBDBS(U[165-1] - U[148-1]) * (U[148-1] - U[165-1]) + COUT * U[148-1]

    # Charge-function of three additional transistors
    G[162-1] += CGS * (U[162-1] - U[98-1])
    G[163-1] += CGD * (U[163-1] - U[98-1])
    G[164-1] += CBDBS(U[164-1] - U[166-1]) * (U[164-1] - U[166-1])
    G[165-1] += CBDBS(U[165-1] - U[148-1]) * (U[165-1] - U[148-1])
    G[166-1] += (CBDBS(U[164-1] - U[166-1]) * (U[166-1] - U[164-1]) + CBDBS(U[170-1] - U[166-1]) * (U[166-1] - U[170-1]) + CLOAD * U[
        166-1])
    G[167-1] += CGS * (U[167-1] - U[18-1])
    G[168-1] += CGD * (U[168-1] - U[18-1])
    G[169-1] += CBDBS(U[169-1] - U[171-1]) * (U[169-1] - U[171-1])
    G[170-1] += CBDBS(U[170-1] - U[166-1]) * (U[170-1] - U[166-1])
    G[171-1] += (CBDBS(U[169-1] - U[171-1]) * (U[171-1] - U[169-1]) + CBDBS(U[175-1] - U[171-1]) * (U[171-1] - U[175-1]) + CLOAD * U[
        171-1])
    G[172-1] += CGS * U[172-1]
    G[173-1] += CGD * U[173-1]
    G[174-1] += CBDBS(U[174-1]) * U[174-1]
    G[175-1] += CBDBS(U[175-1] - U[171-1]) * (U[175-1] - U[171-1])

    return G


def DNOR(N,U,I,G):
    

    G[I-1] += CGS * (U[I-1] - U[I + 4-1])
    G[I + 1-1] += CGD * (U[I + 1-1] - U[I + 4-1])
    G[I + 2-1] += CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 2-1] - U[I + 4-1])
    G[I + 3-1] += CBDBS(U[I + 3-1] - VDD) * U[I + 3-1]
    G[I + 4-1] += (CGS * (U[I + 4-1] - U[I-1]) + CGD * (U[I + 4-1] - U[I + 1-1]) + CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 2-1]) + CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 8-1]) + CBDBS(U[I + 12-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 12-1]) + CLOAD * U[I + 4-1])
    G[I + 5-1] += CGS * U[I + 5-1]
    G[I + 6-1] += CGD * U[I + 6-1]
    G[I + 7-1] += CBDBS(U[I + 7-1]) * U[I + 7-1]
    G[I + 8-1] += CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 8-1] - U[I + 4-1])
    G[I + 9-1] += CGS * U[I + 9-1]
    G[I + 10-1] += CGD * U[I + 10-1]
    G[I + 11-1] += CBDBS(U[I + 11-1]) * U[I + 11-1]
    G[I + 12-1] += CBDBS(U[I + 12-1] - U[I + 4-1]) * (U[I + 12-1] - U[I + 14-1])
    
    return G

def DANDOI(N,U,I,G):
    

    G[I-1] += CGS * (U[I-1] - U[I + 4-1])
    G[I + 1-1] += CGD * (U[I + 1-1] - U[I + 4-1])
    G[I + 2-1] += CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 2-1] - U[I + 4-1])
    G[I + 3-1] += CBDBS(U[I + 3-1] - VDD) * U[I + 3-1]
    G[I + 4-1] += (CGS * (U[I + 4-1] - U[I-1]) + CGD * (U[I + 4-1] - U[I + 1-1]) + \
                 CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 2-1]) + \
                 CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 8-1]) + \
                 CBDBS(U[I + 12-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 12-1]) + \
                 CLOAD * U[I + 4-1])
    G[I + 5-1] += CGS * U[I + 5-1]
    G[I + 6-1] += CGD * U[I + 6-1]
    G[I + 7-1] += CBDBS(U[I + 7-1]) * U[I + 7-1]
    G[I + 8-1] += CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 8-1] - U[I + 4-1])
    G[I + 9-1] += CGS * U[I + 9-1]
    G[I + 10-1] += CGD * U[I + 10-1]
    G[I + 11-1] += CBDBS(U[I + 11-1] - U[I + 13-1]) * (U[I + 11-1] - U[I + 13-1])
    G[I + 12-1] += CBDBS(U[I + 12-1] - U[I + 4-1]) * (U[I + 12-1] - U[I + 4-1])
    G[I + 13-1] += (CBDBS(U[I + 11-1] - U[I + 13-1]) * (U[I + 13-1] - U[I + 11-1]) + \
                  CBDBS(U[I + 17-1] - U[I + 13-1]) * (U[I + 13-1] - U[I + 17-1]) + \
                  CLOAD * U[I + 13-1])
    G[I + 14-1] += CGS * U[I + 14-1]
    G[I + 15-1] += CGD * U[I + 15-1]
    G[I + 16-1] += CBDBS(U[I + 16-1]) * U[I + 16-1]
    G[I + 17-1] += CBDBS(U[I + 17-1] - U[I + 13-1]) * (U[I + 17-1] - U[I + 13-1])

    return G


def DNAND(N,U,I,G):
   

    G[I-1] += CGS * (U[I-1] - U[I + 4-1])
    G[I + 1-1] += CGD * (U[I + 1-1] - U[I + 4-1])
    G[I + 2-1] += CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 2-1] - U[I + 4-1])
    G[I + 3-1] += CBDBS(U[I + 3-1] - VDD) * U[I + 3-1]
    G[I + 4-1] += (CGS * (U[I + 4-1] - U[I-1]) + CGD * (U[I + 4-1] - U[I + 1-1]) + \
                 CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 2-1]) + \
                 CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 8-1]) + \
                 CLOAD * U[I + 4-1])
    G[I + 5-1] += CGS * U[I + 5-1]
    G[I + 6-1] += CGD * U[I + 6-1]
    G[I + 7-1] += CBDBS(U[I + 7-1] - U[I + 9-1]) * (U[I + 7-1] - U[I + 9-1])
    G[I + 8-1] += CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 8-1] - U[I + 4-1])
    G[I + 9-1] += (CBDBS(U[I + 7-1] - U[I + 9-1]) * (U[I + 9-1] - U[I + 7-1]) + \
                 CBDBS(U[I + 13-1] - U[I + 9-1]) * (U[I + 9-1] - U[I + 13-1]) + \
                 CLOAD * U[I + 9-1])
    G[I + 10-1] += CGS * U[I + 10-1]
    G[I + 11-1] += CGD * U[I + 11-1]
    G[I + 12-1] += CBDBS(U[I + 12-1]) * U[I + 12-1]
    G[I + 13-1] += CBDBS(U[I + 13-1] - U[I + 9-1]) * (U[I + 13-1] - U[I + 9-1])

    return G


def DORANI(N,U,I,G):
   

    G[I-1] += CGS * (U[I-1] - U[I + 4-1])
    G[I + 1-1] += CGD * (U[I + 1-1] - U[I + 4-1])
    G[I + 2-1] += CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 2-1] - U[I + 4-1])
    G[I + 3-1] += CBDBS(U[I + 3-1] - VDD) * U[I + 3-1]
    G[I + 4-1] += (CGS * (U[I + 4-1] - U[I-1]) + CGD * (U[I + 4-1] - U[I + 1-1]) + \
                 CBDBS(U[I + 2-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 2-1]) + \
                 CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 4-1] - U[I + 8-1]) + \
                 CLOAD * U[I + 4-1])
    G[I + 5-1] += CGS * U[I + 5-1]
    G[I + 6-1] += CGD * U[I + 6-1]
    G[I + 7-1] += CBDBS(U[I + 7-1] - U[I + 9-1]) * (U[I + 7-1] - U[I + 9-1]) 
    G[I + 8-1] += CBDBS(U[I + 8-1] - U[I + 4-1]) * (U[I + 8-1] - U[I + 4-1])
    G[I + 9-1] += (CBDBS(U[I + 7-1] - U[I + 9-1]) * (U[I + 9-1] - U[I + 7-1]) + \
                 CBDBS(U[I + 13-1] - U[I + 9-1]) * (U[I + 9-1] - U[I + 13-1]) + \
                 CBDBS(U[I + 17-1] - U[I + 9-1]) * (U[I + 9-1] - U[I + 17-1]) + \
                 CLOAD * U[I + 9-1])
    G[I + 10-1] += CGS * U[I + 10-1]
    G[I + 11-1] += CGD * U[I + 11-1]
    G[I + 12-1] += CBDBS(U[I + 12-1]) * U[I + 12-1]
    G[I + 13-1] += CBDBS(U[I + 13-1] - U[I + 9-1]) * (U[I + 13-1] - U[I + 9-1])
    G[I + 14-1] += CGS * U[I + 14-1]
    G[I + 15-1] += CGD * U[I + 15-1]
    G[I + 16-1] += CBDBS(U[I + 16-1]) * U[I + 16-1]
    G[I + 17-1] += CBDBS(U[I + 17-1] - U[I + 9-1]) * (U[I + 17-1] - U[I + 9-1])

    return G

def IDS (NED,VDS, VGS, VBS, VGD, VBD, ierr):
    if VDS > 0.:
        return GDSP (NED,VDS, VGS, VBS,ierr)
    elif VDS == 0.:
        return 0.
    elif VDS < 0.:
        return GDSM (NED,VDS, VGD, VBD,ierr)

    if ierr == -1:
        return


def GDSP (NED,VDS, VGS, VBS,ierr):
    if NED == 0:
        VT0=-2.43
        CGAMMA=.2
        PHI=1.28
        BETA=53.5e-6*CTIME*STIFF
    elif NED == 1:
        VT0=.2
        CGAMMA=0.035
        PHI=1.01
        BETA=4*43.7e-6*CTIME*STIFF
    elif NED == 2:
        VT0=.2
        CGAMMA=0.035
        PHI=1.01
        BETA=8*43.7e-6*CTIME*STIFF
    else:
        VT0=.2
        CGAMMA=0.035
        PHI=1.01
        BETA=12*43.7e-6*CTIME*STIFF

    if(PHI-VBS) < 0. or (PHI < 0.):
        ierr = -1
        return 0.

    #if np.isnan(VBS) == True:
    #    VBS = 0.

    VTE = VT0 + CGAMMA * ( np.sqrt(PHI-VBS) - np.sqrt(PHI) )

    if ( VGS-VTE) <= 0.:
        return 0.
    elif ( VGS-VTE) <= VDS:
        return - BETA * (VGS - VTE)**2. * (1. + DELTA*VDS)
    else:
        return - BETA * VDS * (2.*(VGS - VTE) - VDS) * (1. + DELTA*VDS)


def GDSM (NED,VDS, VGD, VBD, ierr):
    if NED == 0:
        VT0=-2.43
        CGAMMA=.2
        PHI=1.28
        BETA=53.5e-6*CTIME*STIFF
    elif NED == 1:
        VT0=.2
        CGAMMA=0.035
        PHI=1.01
        BETA=4*43.7e-6*CTIME*STIFF
    elif NED == 2:
        VT0=.2
        CGAMMA=0.035
        PHI=1.01
        BETA=8*43.7e-6*CTIME*STIFF
    else:
        VT0=.2
        CGAMMA=0.035
        PHI=1.01
        BETA=12*43.7e-6*CTIME*STIFF


    if(PHI-VBD) < 0. or PHI < 0.:
        ierr = -1
        return 0.

    #if np.isnan(VBD) == True:
    #    VBD = 0.
    
    VTE = VT0 + CGAMMA * ( np.sqrt(PHI-VBD) - np.sqrt(PHI) )

    if ( VGD-VTE) <= 0.:
        return 0.
    elif ( VGD-VTE) <= -VDS:
        return BETA * (VGD - VTE)*(VGD - VTE) * (1. - DELTA*VDS)
    else:
        return - BETA * VDS * (2.0*(VGD - VTE) + VDS) * (1. - DELTA*VDS)



def IBS (VBS):
    if VBS <= 0.:
        return - CURIS * ( np.exp( VBS/VTH ) - 1. )
    else:
        return 0.


def IBD(VBD):
    if VBD <= 0.:
        return - CURIS * ( np.exp( VBD/VTH ) - 1. )
    else:
        return 0.


def CBDBS(V):
    PHIB=0.87

    if V <= 0.:
        return CBD/np.sqrt(1.-V/PHIB)
    else:
        return CBD*(1.+V/(2.*PHIB))

           

def f(t, u):
    ierr = 1.
    x = np.zeros(175); res = np.zeros(175)
    deriv = np.zeros(350)
    
    for i in range(175):
        x[i] = u[i+175]

    res = FCN(175,t,x,res,ierr)

    if ierr == -1:
        return

    for i in range(175):
        deriv[i] = res[i]

    res = GCN(175,x,res)

    for i in range(175):
        deriv[i+175] = u[i] - res[i]


    return deriv
    

t = tt[0]
dt = 1.e-5
tn = tt[-1]

y = np.zeros(175)
u = np.zeros(350)
y[0] = 4.999999999996544
y[1] = 4.999999999999970
y[2] = -2.499999999999975
y[3] = -2.499999999999975
y[4] = 4.999999999996514
y[5] = 0.000000000000000
y[6] = 4.999999999996514
y[7] = -2.499999999999991
y[8] = -2.499999999999975
y[9] = 0.000000000000000
y[10] = 4.999999999996514
y[11] =	-2.499999999999991
y[12] =	-2.499999999999975
y[13] =	0.215858486765796
y[14] =	4.988182208251953
y[15] =	-2.499999999999990
y[16] =	-2.499999999999975
y[17] =	0.204040695017748
y[18] =	0.011817791748026
y[19] =	0.192222903269723
y[20] =	-2.499999999999991
y[21] =	-2.499999999999990
y[22] =	-0.228160951881239
y[23] =	0.204040695017748
y[24] =	-2.499999999999992
y[25] =	-2.499999999999990
y[26] =	-0.228160951881241
y[27] =	0.000000000000000
y[28] =	-0.228160951881239
y[29] =	-2.499999999999991
y[30] =	-2.499999999999992
y[31] =	4.999999999996547
y[32] =	4.999999999999970
y[33] =	-2.499999999999975
y[34] =	-2.499999999999975
y[35] =	4.999999999996517
y[36] =	0.000000000000000
y[37] =	4.999999999996517
y[38] =	-2.499999999999991
y[39] =	-2.499999999999975
y[40] =	0.000000000000000
y[41] =	4.999999999996517
y[42] =	-2.499999999999991
y[43] =	-2.499999999999975
y[44] =	0.215858484247529
y[45] =	4.988182208251953
y[46] =	-2.499999999999990
y[47] =	-2.499999999999975
y[48] =	0.204040692499482
y[49] =	0.011817791748035
y[50] =	0.192222900751447
y[51] =	-2.499999999999991
y[52] =	-2.499999999999990
y[53] =	-0.026041071738432
y[54] =	0.204040692499482
y[55] =	-2.499999999999992
y[56] =	-2.499999999999990
y[57] =	-0.026041071738434
y[58] =	0.000000000000000
y[59] =	-0.026041071738432
y[60] =	-2.499999999999991
y[61] =	-2.499999999999992
y[62] =	0.215858484880918
y[63] =	4.988182208251953
y[64] =	-2.499999999999990
y[65] =	-2.499999999999975
y[66] =	0.204040693132870
y[67] =	0.011817791748026
y[68] =	0.192222901384845
y[69] =	-2.499999999999991
y[70] =	-2.499999999999990
y[71] =	-0.026041071737961
y[72] =	0.204040693132870
y[73] =	-2.499999999999992
y[74] =	-2.499999999999990
y[75] =	-0.026041071737963
y[76] =	0.000000000000000
y[77] =	-0.026041071737961
y[78] =	-2.499999999999991
y[79] =	-2.499999999999992
y[80] =	4.999999999996546
y[81] =	4.999999999999970
y[82] =	-2.499999999999975
y[83] =	-2.499999999999975
y[84] =	4.999999999996516
y[85] =	0.000000000000000
y[86] =	4.999999999996516
y[87] =	-2.499999999999991
y[88] =	-2.499999999999975
y[89] =	0.000000000000000
y[90] =	4.999999999996516
y[91] =	-2.499999999999991
y[92] =	-2.499999999999975
y[93] =	0.215858481060569
y[94] =	4.988182208251953
y[95] =	-2.499999999999990
y[96] =	-2.499999999999975
y[97] =	0.204040689312522
y[98] =	0.011817791748023
y[99] =	0.192222897564498
y[100] = -2.499999999999991
y[101] = -2.499999999999990
y[102] = 4.734672533390068
y[103] = 0.204040689312522
y[104] = -2.499999999999977
y[105] = -2.499999999999990
y[106] = 4.734672533390062
y[107] = 0.000000000000000
y[108] = 4.734672533390068
y[109] = -2.499999999999991
y[110] = -2.499999999999977
y[111] = 4.999999999996870
y[112] = 4.999999999999972
y[113] = -2.499999999999975
y[114] = -2.499999999999975
y[115] = 4.999999999996843
y[116] = -0.025968303070038
y[117] = 4.999999999996843
y[118] = -2.499999999999992
y[119] = -2.499999999999975
y[120] = -0.025968303070040
y[121] = 0.000000000000000
y[122] = -0.025968303070038
y[123] = -2.499999999999991
y[124] = -2.499999999999992
y[125] = 4.999999999997699
y[126] = 4.999999999999980
y[127] = -2.499999999999975
y[128] = -2.499999999999975
y[129] = 4.999999999997678
y[130] = 4.744923533081106
y[131] = 4.999999999997678
y[132] = -2.499999999999977
y[133] = -2.499999999999975
y[134] = 4.744923533081098
y[135] = 0.000000000000000
y[136] = 4.744923533081106
y[137] = -2.499999999999991
y[138] = -2.499999999999977
y[139] = 0.000000000000000
y[140] = 4.744923533081106
y[141] = -2.499999999999991
y[142] = -2.499999999999977
y[143] = 0.215858484844162
y[144] = 4.988182208251953
y[145] = -2.499999999999990
y[146] = -2.499999999999975
y[147] = 0.204040693096114
y[148] = 0.011817791748023
y[149] = 0.192222901348091
y[150] = -2.499999999999991
y[151] = -2.499999999999990
y[152] = 0.204040693096045
y[153] = 0.204040693096107
y[154] = -2.499999999999990
y[155] = -2.499999999999990
y[156] = 0.204040693096037
y[157] = 0.000000000000000
y[158] = 0.204040693096037
y[159] = -2.499999999999991
y[160] = -2.499999999999990
y[161] = -0.026017361873565
y[162] = 0.204040693096114
y[163] = -2.499999999999992
y[164] = -2.499999999999990
y[165] = -0.026017361873568
y[166] = -0.026017590106916
y[167] = -0.026017361873565
y[168] = -2.499999999999992
y[169] = -2.499999999999992
y[170] = -0.026017590106918
y[171] = 0.000000000000000
y[172] = -0.026017590106916
y[173] = -2.499999999999991
y[174] = -2.499999999999992

u = GCN(175,y,u)

for i in range(175):
    u[i+175] = y[i]

M = np.zeros((350,350))
for i in range(175):
    M[i,i] = 1.
for i in range(175, 350):
    M[i,i] = 0.

var_index = np.zeros(350)
for i in range(350):
    var_index[i] = 1.


