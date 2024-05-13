import numpy as np

# Charge Pump

tt = np.empty(41)
tt[0] = 0.
tt[1] = 50.e-9
tt[2] = 60.e-9
tt[3] = 110.e-9
for i in range(4,40):
    tt[i] = tt[i%4] + int(i/4) * 120.e-9
tt[40] = 1200.e-9


def qgate(vgb, vgs, vgd):
    vt0 = 0.20; gamma = 0.350e-1; phi = 1.010; cox = 4.0e-12
    if (vgs-vgd) <= 0.:
        ugs = vgd
        ugd = vgs
    else:
        ugs = vgs
        ugd = vgd

    ugb = vgb
    ubs = ugs - ugb
    
    if np.isnan(ubs) == True:
        ubs = 0.
    
    vfb = vt0 - gamma * np.sqrt(phi) - phi
    vte = vt0 + gamma *( np.sqrt(phi - ubs) - np.sqrt(phi) )

    if ugb <= vfb:
        return cox * (ugb - vfb)
    elif ugb > vfb and ugs <= vte:
        return cox * gamma * ( np.sqrt((gamma/2.)**2. + ugb - vfb) - gamma/2. )
    elif ugb  > vfb and ugs > vte:
        ugst = ugs - vte
        if ugd > vte:
            ugdt = ugd - vte
        else:
            ugdt = 0.
        return cox * ( (2./3.) * (ugdt + ugst - ((ugdt * ugst)/(ugdt+ugst)) ) + gamma * np.sqrt(phi-ubs) )
    #else:
    #    return 0.


def qbulk(vgb,vgs,vgd):
    vt0 = 0.20; gamma = 0.350e-1; phi = 1.010; cox = 4.0e-12

    if (vgs - vgd) <= 0.:
        ugs = vgd
        ugd = vgs
    else:
        ugs = vgs
        ugd = vgd

    ugb = vgb
    ubs = ugs - ugb

    vfb = vt0 - gamma * np.sqrt(phi) - phi
    vte = vt0 + gamma *( np.sqrt(phi - ubs) - np.sqrt(phi) )

    if ugb <= vfb:
        return - cox * (ugb - vfb)
    elif (ugb > vfb) and (ugs <= vte):
        return - cox * gamma * ( np.sqrt((gamma/2.)**2. + ugb - vfb) - gamma/2. )
    elif ugb > vfb and ugs > vte:
        return - cox * gamma * np.sqrt(phi-ubs)
    #else:
    #    return 0.
        

def qsrc(vgb, vgs, vgd):
    vt0 = 0.20; gamma = 0.350e-1; phi = 1.010; cox = 4.0e-12

    if (vgs - vgd) <= 0.:
        ugs = vgd
        ugd = vgs
    else:
        ugs = vgs
        ugd = vgd

    ugb = vgb
    ubs = ugs - ugb

    if np.isnan(ubs) == True:
        ubs = 0.
        
    vfb = vt0 - gamma * np.sqrt(phi) - phi
    vte = vt0 + gamma *( np.sqrt(phi - ubs) - np.sqrt(phi) )

    if ugb <= vfb:
        return 0.
    elif (ugb > vfb) and (ugs <= vte):
        return 0.
    elif ( ugb > vfb) and (ugs > vte):
        ugst = ugs - vte
        if ugd >= vte:
            ugdt = ugd - vte
        else:
            ugdt = 0.
        return - cox * (1./3.) * (ugdt + ugst - ((ugdt * ugst)/(ugdt+ugst)) )
    #else:
    #    return 0.


def qdrain(vgb,vgs,vgd):
    vt0 = 0.20; gamma = 0.350e-1; phi = 1.010; cox = 4.0e-12

    if (vgs - vgd) <= 0.:
        ugs = vgd
        ugd = vgs
    else:
        ugs = vgs
        ugd = vgd

    ugb = vgb
    ubs = ugs - ugb

    if np.isnan(ubs) == True:
        ubs = 0.

    vfb = vt0 - gamma * np.sqrt(phi) - phi
    vte = vt0 + gamma *( np.sqrt(phi - ubs) - np.sqrt(phi) )

    if ugb <= vfb:
        return 0.
    elif (ugb > vfb) and (ugs <= vte):
        return 0.
    elif ( ugb > vfb) and (ugs > vte):
        ugst = ugs - vte
        if ugd >= vte:
            ugdt = ugd - vte
        else:
            ugdt = 0.
        return - cox * (1./3.) * (ugdt + ugst - ((ugdt * ugst)/(ugdt+ugst)) )
    #else:
    #    return 0.


def vin(t):
    vhigh = 20.

    deltat = 120.0e-9
    t1 = 50.0e-9
    t2 = 60.0e-9
    t3 = 110.0e-9

    dummy = t % deltat
    if dummy < t1:
        # 0 <= t' < 50ns
        return 0.0
    else:
        if dummy < t2:
            # 50ns <= t' < 60ns
            return (dummy - t1) * 0.10e+9 * vhigh
        else:
            if dummy < t3:
                # 60ns <= t' < 110ns
                return vhigh
            else:
                # 110ns <= t' < 120ns
                return (deltat - dummy) * 0.10e+9 * vhigh
            

def f(t, u):
            
    capd = 0.40e-12; caps = 1.60e-12
    
    deriv = [-u[8],
             0.,
             0.,
             -u[5] + vin(t),
             u[0] - qgate(u[5], u[5]-u[6], u[5]-u[7]),
             u[1] - caps*u[6],
             u[2] - qsrc(u[5], u[5]-u[6], u[5]-u[7]),
             u[3] - capd*u[7],
             u[4] - qdrain(u[5], u[5]-u[6], u[5]-u[7])]
    return np.array(deriv)
    

t = tt[0]
dt = 4.4e-13
tn = tt[-1]
u = [qgate(0.,0.,0.),
     0.,
     qsrc(0.,0.,0.),
     0.,
     qdrain(0.,0.,0.),
     0.,
     0.,
     0.,
     0.]

M = np.zeros((9,9))

M[0,0] = 1.
M[1,1] = 1.
M[1,2] = 1.
M[2,3] = 1.
M[2,4] = 1.


var_index = np.array([1,1,1,1,1,1,1,1,2])

u = np.array(u)
