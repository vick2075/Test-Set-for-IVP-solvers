import numpy as np


# Stiff Pendulum problem
        

def f(t, u):
    hold = np.sqrt(u[0]**2+u[1]**2)
    deriv = [u[2],
             u[3],
             -u[0]*u[4],
             -u[1]*u[4]-1.0,
             1.e-6 * u[4]-(hold-1.0)/hold]
    return np.array(deriv)
    

t = 0
dt = 0.0001
tn = 10.
u = [1., 0., 0., 0., 0.]

M = np.array([[1.,0,0,0,0],
              [0,1.,0,0,0],
              [0,0,1.,0,0],
              [0,0,0,1.,0],
              [0,0,0,0,0]])

var_index = np.array([1,1,2,2,3])

u = np.array(u)
