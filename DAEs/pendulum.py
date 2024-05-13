import numpy as np

# Pendulum problem
        

def f(t, u):
    deriv = [u[2],
             u[3],
             -u[0]*u[4],
             -u[1]*u[4]-1.0,
             u[0]*u[0]+u[1]*u[1]-1.0]
    return np.array(deriv)
    

t = 0
dt = 0.0001
tn = 1.
u = [1., 0., 0., 1., 1.]

M = np.array([[1.,0,0,0,0],
              [0,1.,0,0,0],
              [0,0,1.,0,0],
              [0,0,0,1.,0],
              [0,0,0,0,0]])

var_index = np.array([1,1,2,2,3])

u = np.array(u)
