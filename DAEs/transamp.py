import numpy as np

# Transistor Amplifier
C1=1.e-6; C2=2.e-6; C3=3.e-6; C4=4.e-6; C5=5.e-6
R0=1.e3
R1=R2=R3=R4=R5=R6=R7=R8=R9=9.e3
Ub=6.
uf = 0.026; alpha=0.99; beta=1.e-6


def f(t, u):
    if (u[1]-u[2])/uf > 300 or (u[4]-u[5])/uf > 300:
        return


    fac1 = beta*(np.exp((u[1]-u[2])/uf)-1)
    fac2 = beta*(np.exp((u[4]-u[5])/uf)-1)
    
    deriv = [(u[0] - 0.1*np.sin(200*np.pi*t))/R0,
             u[1]/R1+(u[1]-Ub)/R2+(1-alpha)*fac1,
             u[2]/R3-fac1,
             (u[3]-Ub)/R4+alpha*fac1,
             u[4]/R5+(u[4]-Ub)/R6+(1-alpha)*fac2,
             u[5]/R7-fac2,
             (u[6]-Ub)/R8+alpha*fac2,
             u[7]/R9]
    return np.array(deriv)

t = 0
dt = 0.01
tn = 0.2
u = [0., Ub/(R2/R1+1.), Ub/(R2/R1+1.), Ub, Ub/(R6/R5+1.), Ub/(R6/R5+1.), Ub, 0.]

M = np.array([[-C1,C1,0,0,0,0,0,0],
              [C1,-C1,0,0,0,0,0,0],
              [0,0,-C2,0,0,0,0,0],
              [0,0,0,-C3,C3,0,0,0],
              [0,0,0,C3,-C3,0,0,0],
              [0,0,0,0,0,-C4,0,0],
              [0,0,0,0,0,0,-C5,C5],
              [0,0,0,0,0,0,C5,-C5]])

var_index = np.array([0,0,0,0,0,0,0,0])

u = np.array(u)
