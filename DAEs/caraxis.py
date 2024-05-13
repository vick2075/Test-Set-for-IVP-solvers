import numpy as np

# Car Axis problem
eps = 1.e-2; mm = 10.

def f(t, u):
    L = 1.; L0 = 0.5; r = 0.1
    w = 10.; g = 1.; yb = r*np.sin(w*t); xb = np.sqrt(L*L-yb*yb)
    xl = u[0]; yl = u[1]; xr = u[2]; yr = u[3]; lam1 = u[8]; lam2 = u[9]
    Ll = np.sqrt(xl**2+yl**2)
    Lr = np.sqrt((xr-xb)**2+(yr-yb)**2)
    
    deriv = [u[4],
             u[5],
             u[6],
             u[7],
             (L0-Ll)*xl/Ll + lam1*xb+2.*lam2*(xl-xr),
             (L0-Ll)*yl/Ll + lam1*yb+2.*lam2*(yl-yr)-mm*eps*eps*g/2.,
             (L0-Lr)*(xr-xb)/Lr -2.* lam2*(xl-xr),
             (L0-Lr)*(yr-yb)/Lr -2.*lam2*(yl-yr)-mm*eps*eps*g/2.,
             xb*xl+yb*yl,
             (xl-xr)**2+(yl-yr)**2-L*L]
    return np.array(deriv)
    

t = 0
dt = 0.0001
tn = 3.
u = [0., 0.5, 1., 0.5, -0.5/1., 0., -0.5/1., 0., 0., 0.]

ss = mm*eps*eps/2.
M = np.array([[1.,0,0,0,0,0,0,0,0,0],
              [0,1.,0,0,0,0,0,0,0,0],
              [0,0,1.,0,0,0,0,0,0,0],
              [0,0,0,1.,0,0,0,0,0,0],
              [0,0,0,0,ss,0,0,0,0,0],
              [0,0,0,0,0,ss,0,0,0,0],
              [0,0,0,0,0,0,ss,0,0,0],
              [0,0,0,0,0,0,0,ss,0,0],
              [0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0]])

var_index = np.array([1,1,1,1,2,2,2,2,3,3])

u = np.array(u)
