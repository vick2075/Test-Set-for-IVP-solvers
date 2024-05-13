import numpy as np

# van der Pol

def f(t, u):
    y = 1000.
    #y = 1.0e+6
    
    deriv = [u[1] ,\
             y *(1.0 - u[0]**2) * u[1] - u[0] ] 
    return np.array(deriv)

t = 0

dt = 0.1

tn = 3000.
#tn = 6.3

u = [2.0, 0.0]
#u = [0.0, np.sqrt(3.)] #[1.0, 1.0]#

u = np.array(u)
