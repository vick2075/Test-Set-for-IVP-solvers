import numpy as np


# Pollution

def f(t, u):
    k =[0.35, 26.6, 12300.0, 0.00086, 0.00082, 15000.0, 0.00013, 24000.0, 16500.0, 9000.0, 0.022, 12000.0, 1.88, \
        16300.0, 4.8e+6, 0.00035, 0.0175, 1.0e+8, 4.44e+11, 1240.0, 2.1, 5.78, 0.0474, 1780.0, 3.12]
    
    deriv = [-k[0]*u[0] - k[9]*u[10]*u[0] - k[13]*u[0]*u[5] - k[22]*u[0]*u[3] - \
             k[23]*u[18]*u[0] + k[1]*u[1]*u[3] + k[2]*u[4]*u[1] + k[8]*u[10]*u[1] + \
             k[10]*u[12] + k[11]*u[9]*u[1] + k[21]*u[18] + k[24]*u[19],

             -k[1]*u[1]*u[3] - k[2]*u[4]*u[1] - k[8]*u[10]*u[1] - k[11]*u[9]*u[1] + \
             k[0]*u[0] + k[20]*u[18],

             -k[14]*u[2] + k[0]*u[0] + k[16]*u[3] + k[18]*u[15] + k[21]*u[18],

             -k[1]*u[1]*u[3] - k[15]*u[3] - k[16]*u[3] -k[22]*u[0]*u[3] + k[14]*u[2],

             -k[2]*u[4]*u[1] + 2*(k[3]*u[6]) + k[5]*u[6]*u[5] + k[6]*u[8] + \
             k[12]*u[13] + k[19]*u[16]*u[5],

             -k[5]*u[6]*u[5] - k[7]*u[8]*u[5] - k[13]*u[0]*u[5] - k[19]*u[16]*u[5] + \
             k[2]*u[4]*u[1] + 2*(k[17]*u[15]),

             -k[3]*u[6] - k[4]*u[6] - k[5]*u[6]*u[5] + k[12]*u[13],

             k[3]*u[6] + k[4]*u[6] + k[5]*u[6]*u[5] + k[6]*u[8],

             -k[6]*u[8] - k[7]*u[8]*u[5],

             -k[11]*u[9]*u[1] + k[6]*u[8] + k[8]*u[10]*u[1],

             -k[8]*u[10]*u[1] - k[9]*u[10]*u[0] + k[7]*u[8]*u[5] + k[10]*u[12],

             k[8]*u[10]*u[1],

             -k[10]*u[12] + k[9]*u[10]*u[0],

             -k[12]*u[13] + k[11]*u[9]*u[1],

             k[13]*u[0]*u[5],

             -k[17]*u[15] - k[18]*u[15] + k[15]*u[3],

             -k[19]*u[16]*u[5],

             k[19]*u[16]*u[5],

             -k[20]*u[18] - k[21]*u[18] - k[23]*u[18]*u[0] + k[22]*u[0]*u[3] + k[24]*u[19],

             -k[24]*u[19] + k[23]*u[18]*u[0]]
    return np.array(deriv)


t = 0
dt = 0.1
tn = 60.
u = [0.0, 0.2, 0.0, 0.04, 0.0, 0.0, 0.1, 0.3, 0.017, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007, 0.0, 0.0, 0.0]
u = np.array(u)
