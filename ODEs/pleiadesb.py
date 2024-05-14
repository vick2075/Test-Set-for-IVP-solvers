import numpy as np


# Pleiades problem
G = 1.0
bod = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
ma = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

### Please note that code for plotting for the below
### two types of code will be different!

## As per Fortran codes:
''' 
def f(t, u):
    deriv=np.zeros(28)
    for i in range(1, 7+1):
        sumx = 0.; sumy = 0.
        for j in range(1, 7+1):
            mj = j
            rij = (u[i-1]-u[j-1])**2 + (u[i+7-1]-u[j+7-1])**2
            rij32 = rij**(3./2.)
            if j != i:
                sumx += mj*(u[j-1]-u[i-1])/rij32
                sumy += mj*(u[j+7-1]-u[i+7-1])/rij32
        deriv[i+14-1] = sumx
        deriv[i+21-1] = sumy
    for i in range(1, 14+1):
        deriv[i-1] = u[i+14-1]
    
    return deriv

t = 0
dt = 0.075#0.01
tn = 15. #3.

u = [3.0, 3.0, -1.0, -3.0, 2.0, -2.0, 2.0, 3.0, -3.0, 2.0, 0.0, 0.0, -4.0, 4.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 1.75, -1.5, 0.0, 0.0, 0.0, -1.25, 1.0, 0.0, 0.0]

u = np.array(u)
'''


## Normal n-body system of ODEs
#'''        
def f(t, u):
    deriv = []

    m = ma
    N = len(m)

    axx = ayy = azz = 0
    
    for i in range(N):
        ioffset = i * 6

        ax = ay = az = 0

        for j in range(N):
            joffset = j * 6
            
            if i != j:
                
                dx =  u[joffset+0] - u[ioffset+0]
                dy =  u[joffset+1] - u[ioffset+1]
                dz =  u[joffset+2] - u[ioffset+2]
                r = (dx**2 + dy**2 + dz**2) ** 0.5
                ax += (G * m[j] / (r**3)) * dx
                ay += (G * m[j] / (r**3)) * dy
                az += (G * m[j] / (r**3)) * dz
            
        axx = ax; ayy = ay; azz = az
        
        deriv.insert(ioffset+0, u[ioffset+3])
        deriv.insert(ioffset+1, u[ioffset+4])
        deriv.insert(ioffset+2, u[ioffset+5])
        deriv.insert(ioffset+3, axx)
        deriv.insert(ioffset+4, ayy)
        deriv.insert(ioffset+5, azz)
    return np.array(deriv)


t = 0
dt = 0.075#0.01
tn = 15. #3.

plp = [[3.0, 3.0, 0.0],[3.0, -3.0, 0.0],[-1.0, 2.0, 0.0],[-3.0, 0.0, 0.0],[2.0, 0.0, 0.0],[-2.0, -4.0, 0.0],[2.0, 4.0, 0.0]] # pos of bodies
plv = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, -1.25, 0.0],[0.0, 1.0, 0.0],[1.75, 0.0, 0.0],[-1.5, 0.0, 0.0]] # vel of bodies


n,m = np.shape(plp)
u = []
for i in range(n):
    for j in range(m):
        u.append(plp[i][j])
    for k in range(m):
        u.append(plv[i][k])

u = np.array(u)
#'''
