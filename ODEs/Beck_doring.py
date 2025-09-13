import numpy as np


# Becker–Döring ODE system

# Parameters
ND = 5000  # Number of cluster sizes
RHO = 7.5

# Initialize detachment rates B(k), for k = 2 to ND
B = np.zeros(ND + 1)  # B[0] and B[1] are unused for 1-based indexing
R0 = 0.0
for i in range(1, ND):
    R1 = i**(2.0 / 3.0)
    B[i+1] = np.exp((2 * i - 1.0) / (R1**2 + R1 * R0 + R0**2))
    R0 = R1


def f(t, u):
    deriv = np.zeros_like(u)
    y1 = u[0]  # c1
    ND = len(u)
    
    # Compute last flux J_{ND-1}
    ajo = u[ND - 2] * y1 - B[ND] * u[ND - 1]
    deriv[ND - 1] = ajo
    sumj = ajo
    
    # Loop from ND-1 down to 2
    for i in range(ND - 2, 0, -1):
        aj = u[i - 1] * y1 - B[i + 1] * u[i]
        deriv[i] = aj - ajo
        sumj += aj
        ajo = aj

    # Monomer equation: dc1/dt = -sum of all J_k
    deriv[0] = -sumj - aj
    return deriv


t = 0
dt = 6.250e+13
tn = 1.e15

# Initial condition: all monomers, no clusters
u = np.zeros(ND)
u[0] = RHO  # All mass in c1


