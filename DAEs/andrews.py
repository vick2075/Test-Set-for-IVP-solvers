import numpy as np


# Andrews'Squeezing Mechanism
        

def f(t, u):
    m1=.04325;m2=.00365;m3=.02373;m4=.00706
    m5=.07050;m6=.00706;m7=.05498
    xa=-.06934;ya=-.00227
    xb=-0.03635;yb=.03273
    xc=.014;yc=.072;c0=4530
    i1=2.194e-6;i2=4.410e-7;i3=5.255e-6;i4=5.667e-7
    i5=1.169e-5;i6=5.667e-7;i7=1.912e-5
    d=28.e-3;da=115.e-4;e=2.e-2;ea=1421.e-5
    rr=7.e-3;ra=92.e-5;l0=7785.e-5
    ss=35.e-3;sa=1874.e-5;sb=1043.e-5;sc=18.e-3;sd=2.e-2
    ta=2308.e-5;tb=916.e-5;uu=4.e-2;ua=1228.e-5;ub=449.e-5
    zf=2.e-2;zt=4.e-2;fa=1421.e-5;mom=33.e-3

    m = np.zeros((7,7))
    ff = np.zeros(7)
    gp = np.zeros((6,7))
    g = np.zeros(6)

    deriv = np.zeros(27)

    sibe = np.sin(u[0])
    sith = np.sin(u[1])
    siga = np.sin(u[2])
    siph = np.sin(u[3])
    side = np.sin(u[4])
    siom = np.sin(u[5])
    siep = np.sin(u[6])

    cobe = np.cos(u[0])
    coth = np.cos(u[1])
    coga = np.cos(u[2])
    coph = np.cos(u[3])
    code = np.cos(u[4])
    coom = np.cos(u[5])
    coep = np.cos(u[6])

    sibeth = np.sin(u[0]+u[1])
    siphde = np.sin(u[3]+u[4])
    siomep = np.sin(u[5]+u[6])

    cobeth = np.cos(u[0]+u[1])
    cophde = np.cos(u[3]+u[4])
    coomep = np.cos(u[5]+u[6])

    bep = u[7]
    thp = u[8]
    php = u[10]
    dep = u[11]
    omp = u[12]
    epp = u[13]

    for j in range(7):
        for i in range(7):
            m[i,j] = 0.

    m[0,0] = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 + i2
    m[1,0] = m2*(da**2-da*rr*coth) + i2
    m[1,1] = m2*da**2 + i2
    m[2,2] = m3*(sa**2+sb**2) + i3
    m[3,3] = m4*(e-ea)**2 + i4
    m[4,3] = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
    m[4,4] = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2) + i4 + i5
    m[5,5] = m6*(zf-fa)**2 + i6
    m[6,5] = m6*((zf-fa)**2-uu*(zf-fa)*siom) + i6
    m[6,6] = m6*((zf-fa)**2-2*uu*(zf-fa)*siom+uu**2) + m7*(ua**2+ub**2) + i6 + i7

    for j in range(1,7):
        for i in range(j):
            m[i,j] = m[j,i]

    xd = sd*coga + sc*siga + xb
    yd = sd*siga - sc*coga + yb
    lang  = np.sqrt ((xd-xc)**2 + (yd-yc)**2)
    force = - c0 * (lang - l0)/lang
    fx = force * (xd-xc)
    fy = force * (yd-yc)
    ff[0] = mom - m2*da*rr*thp*(thp+2*bep)*sith
    ff[1] = m2*da*rr*bep**2*sith
    ff[2] = fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga)
    ff[3] = m4*zt*(e-ea)*dep**2*coph
    ff[4] = - m4*zt*(e-ea)*php*(php+2*dep)*coph
    ff[5] = - m6*uu*(zf-fa)*epp**2*coom
    ff[6] = m6*uu*(zf-fa)*omp*(omp+2*epp)*coom

    for j in range(7):
        for i in range(6):
            gp[i,j] = 0.

    gp[0,0] = - rr*sibe + d*sibeth
    gp[0,1] = d*sibeth
    gp[0,2] = - ss*coga
    gp[1,0] = rr*cobe - d*cobeth
    gp[1,1] = - d*cobeth
    gp[1,2] = - ss*siga
    gp[2,0] = - rr*sibe + d*sibeth
    gp[2,1] = d*sibeth
    gp[2,3] = - e*cophde
    gp[2,4] = - e*cophde + zt*side
    gp[3,0] = rr*cobe - d*cobeth
    gp[3,1] = - d*cobeth
    gp[3,3] = - e*siphde
    gp[3,4] = - e*siphde - zt*code
    gp[4,0] = - rr*sibe + d*sibeth
    gp[4,1] = d*sibeth
    gp[4,5] = zf*siomep
    gp[4,6] = zf*siomep - uu*coep
    gp[5,0] = rr*cobe - d*cobeth
    gp[5,1] = - d*cobeth
    gp[5,5] = - zf*coomep
    gp[5,6] = - zf*coomep - uu*siep

    g[0] = rr*cobe - d*cobeth - ss*siga - xb
    g[1] = rr*sibe - d*sibeth + ss*coga - yb
    g[2] = rr*cobe - d*cobeth - e*siphde - zt*code - xa
    g[3] = rr*sibe - d*sibeth + e*cophde - zt*side - ya
    g[4] = rr*cobe - d*cobeth - zf*coomep - uu*siep - xa
    g[5] = rr*sibe - d*sibeth - zf*siomep + uu*coep - ya

    for i in range(14):
        deriv[i] = u[i+7]

    for i in range(14,21):
        deriv[i] = -ff[i-14]
        for j in range(7):
            deriv[i] = deriv[i]+m[i-14,j]*u[j+14]
        for j in range(6):
            deriv[i] = deriv[i]+gp[j,i-14]*u[j+21]

    for i in range(21, 27):
        deriv[i] = g[i-21]
    
    
    return deriv
    

t = 0
dt = 1.e-06
tn = 3.e-2
u = [-0.0617138900142764496358948458001,
     0.,
     0.455279819163070380255912382449,
     0.222668390165885884674473185609,
     0.487364979543842550225598953530,
     -0.222668390165885884674473185609,
     1.23054744454982119249735015568,
     0.,
     0.,
     0.,
     0.,
     0.,
     0.,
     0.,
     14222.4439199541138705911625887,
     -10666.8329399655854029433719415,
     0.,
     0.,
     0.,
     0.,
     0.,
     98.5668703962410896057654982170,
     -6.12268834425566265503114393122,
     0.,
     0.,
     0.,
     0.]

M = np.zeros((27,27))
for i in range(27):
    M[i,i] = 0.
for i in range(14):
    M[i,i] = 1.

var_index = np.empty(27)
for i in range(7):
    var_index[i] = 1
for i in range(7,14):
    var_index[i] = 2
for i in range(14,27):
    var_index[i] = 3

u = np.array(u)
