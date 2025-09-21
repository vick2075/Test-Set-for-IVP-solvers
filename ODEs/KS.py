# Kuramoto-Sivashinsky Equation

# ---- Parameters
MMM = 9                 # N=1024
NH = 2 ** MMM           # 512
N = 2 * NH              # 1024
ND = 2 * NH - 2         # 1022
QQ = 0.025              # spectral scale
L = (2.0 * np.pi) / QQ  # domain length ≈251.327
UZERO = 0.0

# Precompute diagonals (linear operators, constant)
j_arr = np.arange(1, NH)
diag_arr = (QQ * j_arr)**2 * (1.0 - (QQ * j_arr)**2) # Linear term: j^2 q^2 - j^4 q^4
k_arr = QQ * j_arr  # For nonlinear term: i q j

#@njit
def f(t, u):
    # Build complex Fourier coeffs (averaged convention)
    c = np.zeros(NH, dtype=np.complex128)
    c[0] = UZERO
    c[1:] = u[::2] + 1j * u[1::2]

    # Include Nyquist=0
    c_nyq = np.zeros(NH + 1, dtype=np.complex128)
    c_nyq[:NH] = c

    # Inverse FFT to real space (adjust for sign convention and scaling)
    u_real = np.fft.irfft(np.conj(c_nyq) * N, n=N)

    # Nonlinear term in real space: u^2 / 2
    w_real = u_real**2 / 2.0

    # Forward FFT (averaged, adjust sign to match Hairer)
    w_cv = np.fft.rfft(w_real) / N
    w_hairer = np.conj(w_cv)

    # Extract to packed U [Re0, 0, Re1, Im1, ..., Re_{NH-1}, Im_{NH-1}]
    U = np.zeros(N)
    U[0] = np.real(w_hairer[0])
    U[1] = 0.0
    U[2::2] = np.real(w_hairer[1:NH])
    U[3::2] = np.imag(w_hairer[1:NH])

    # Compute derivatives (vectorized)
    f_Re = diag_arr * u[::2] + k_arr * U[3::2]
    f_Im = diag_arr * u[1::2] - k_arr * U[2::2]

    deriv = np.zeros(ND)
    deriv[::2] = f_Re
    deriv[1::2] = f_Im
    return deriv

t = 0.
dt = 1.e-6
tn = 100.

# Initial real-space u on [0, L] (scale X from [0,1] to [0,L])
AN = N
delx = 1.0 / AN  # adjusted spacing
u = np.zeros(N, dtype=np.float64)
for i in range(1, N + 1):
    X = delx * (i - 1)
    U1 = min(X - 0.0, 0.1 - X)
    U2 = 20.0 * (X - 0.2) * (0.3 - X)
    U3 = min(X - 0.6, 0.7 - X)
    U4 = min(X - 0.9, 1.0 - X)
    m = max(0.0, U1, U2, U3, U4)
    u[i - 1] = 16.0 * m

# Fourier transform initial (using new convention)
uu_cv = np.fft.rfft(u) / N
uu_hairer = np.conj(uu_cv)
UZERO = np.real(uu_hairer[0])
Y = np.zeros(ND, dtype=np.float64)
Y[::2] = np.real(uu_hairer[1:NH])
Y[1::2] = np.imag(uu_hairer[1:NH])
