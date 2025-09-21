# 2D Brusselator
        
NS = 128
NSSQ = NS * NS
NSNSM1 = NS * (NS - 1)
N = 2 * NSSQ
ALF = 1.0e-1
ALPH = ALF

PI = np.pi
NDIM = 2
NF = [NS, NS]
TCOS = np.zeros(512)

for i in range(NS):
    TCOS[i - 1] = 2 * np.cos(PI * (i - 1) * 2.0 / NS) * ALF * NSSQ

NSF = NS
NSSQF = NSSQ
ANS = NS

@njit
def f(t, u):
    """Compute du/dt for the 2D Brusselator with diffusion."""
    deriv = np.zeros_like(u)

    RADSQ = 0.1**2
    BET = 5.0 if t >= 1.1 else 0.0

    for i in range(NSSQ):
        # Compute 2D indices
        row = i // NS  # IY-1
        col = i % NS   # IX-1

        # Periodic boundary indices
        left  = row * NS + (col - 1) % NS
        right = row * NS + (col + 1) % NS
        up    = ((row - 1) % NS) * NS + col
        down  = ((row + 1) % NS) * NS + col

        # Variables
        UIJ = u[i]
        VIJ = u[i + NSSQ]

        ULEFT = u[left]
        URIGHT = u[right]
        UUP = u[up]
        UDOWN = u[down]

        VLEFT = u[left + NSSQ]
        VRIGHT = u[right + NSSQ]
        VUP = u[up + NSSQ]
        VDOWN = u[down + NSSQ]

        # Diffusion term
        lap_u = ULEFT + URIGHT + UUP + UDOWN - 4.0 * UIJ
        lap_v = VLEFT + VRIGHT + VUP + VDOWN - 4.0 * VIJ

        deriv[i] = 1.0 + UIJ**2 * VIJ - 4.4 * UIJ + ALF * NSSQ * lap_u
        deriv[i + NSSQ] = 3.4 * UIJ - UIJ**2 * VIJ + ALF * NSSQ * lap_v

        # Inhomogeneity
        XX = (col + 1) / ANS
        YY = (row + 1) / ANS
        if (XX - 0.3)**2 + (YY - 0.6)**2 <= RADSQ:
            deriv[i] += BET 
    return deriv

t = 0
dt = 0.01
tn = 6.1 # 11.5 # 1.5 # 

u = np.zeros(N)

for j in range(NS):
    YY = (j) / ANS
    for i in range(NS):
        u[j * NS + i] = 22.0 * YY * (1.0 - YY)**1.5
        
for i in range(NS):
    XX = (i ) / ANS
    for j in range(NS):
        u[j * NS + i + NSSQ] = 27.0 * XX * (1.0 - XX)**1.5
