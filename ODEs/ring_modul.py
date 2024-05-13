import numpy as np


# Ring Modulator Problem

def f(t, u):
    c=1.6e-8; cs=2.e-12; cp=1.e-8
    r=25.e+3; rp=50.
    lh=4.45; ls1=2.e-3; ls2=5.e-4; ls3=5.e-4
    rg1=36.3; rg2=17.3; rg3=17.3; ri=5.e+1; rc=6.e+2
    gamma=40.67286402e-9; delta=17.7493332

    uin1   = 0.5*np.sin(2.e+3*np.pi*t)
    uin2   = 2.*np.sin(2.e+4*np.pi*t)
    ud1    = +u[2]-u[4]-u[6]-uin2
    ud2    = -u[3]+u[5]-u[6]-uin2
    ud3    = +u[3]+u[4]+u[6]+uin2
    ud4    = -u[2]-u[5]+u[6]+uin2
    
    if (delta*max(ud1,ud2,ud3,ud4))> 300.0:
        return
        
    qud1   = gamma*(np.exp(delta*ud1)-1.)
    qud2   = gamma*(np.exp(delta*ud2)-1.)
    qud3   = gamma*(np.exp(delta*ud3)-1.)
    qud4   = gamma*(np.exp(delta*ud4)-1.)
    
    deriv = [(u[7]-0.5*u[9]+0.5*u[10]+u[13]-u[0]/r)/c, \
             (u[8]-0.5*u[11]+0.5*u[12]+u[14]-u[1]/r)/c, \
             (u[9]-qud1+qud4)/cs, \
             (-u[10]+qud2-qud3)/cs, \
             (u[11]+qud1-qud3)/cs, \
             (-u[12]-qud2+qud4)/cs, \
             (-u[6]/rp+qud1+qud2-qud3-qud4)/cp, \
             -u[0]/lh, \
             -u[1]/lh, \
             (0.5*u[0]-u[2]-rg2*u[9])/ls2, \
             (-0.5*u[0]+u[3]-rg3*u[10])/ls3, \
             (0.5*u[1]-u[4]-rg2*u[11])/ls2, \
             (-0.5*u[1]+u[5]-rg3*u[12])/ls3, \
              (-u[0]+uin1-(ri+rg1)*u[13])/ls1, \
             (-u[1]-(rc+rg1)*u[14])/ls1]
    return np.array(deriv)

t = 0.
dt = 1.e-9
tn = 1.e-3  
u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
u = np.array(u)
