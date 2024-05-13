import numpy as np


# Pleiades problem
        
au = 149597870700.
day = 3600. * 24.
G1 = 6.67430e-11


G = 1.0
bod = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
ma = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]


    
ma = np.array(ma)
N = len(bod)
mass = np.vstack(np.split(ma, N))
     
def f(t, u):
    [pos, vel] = np.hsplit(u, 2)
    xx = pos[:,0:1]
    yy = pos[:,1:2]
    zz = pos[:,2:3]
    dx = xx.T - xx
    dy = yy.T - yy
    dz = zz.T - zz
    inv_r3 = (dx**2 + dy**2 + dz**2 )
    inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)
    ax = G * (dx * inv_r3) @ mass
    ay = G * (dy * inv_r3) @ mass
    az = G * (dz * inv_r3) @ mass
    a = np.hstack((vel, np.hstack((ax,ay,az))))
    return a



t = 0
dt = 0.075#0.01
tn= 15. # 3. # 

plp = [[3.0, 3.0, 0.0],[3.0, -3.0, 0.0],[-1.0, 2.0, 0.0],[-3.0, 0.0, 0.0],[2.0, 0.0, 0.0],[-2.0, -4.0, 0.0],[2.0, 4.0, 0.0]] # pos of bodies
plv = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, -1.25, 0.0],[0.0, 1.0, 0.0],[1.75, 0.0, 0.0],[-1.5, 0.0, 0.0]] # vel of bodies



n,m = np.shape(plp)
pos = np.ones((n, 3))
vel = np.ones((n, 3))
for i in range(n):
    for j in range(m):
        pos[i,j] = plp[i][j]
        vel[i,j] = plv[i][j]
        
u = np.hstack((pos, vel))

