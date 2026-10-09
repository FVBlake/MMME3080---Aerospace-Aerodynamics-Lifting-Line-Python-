import numpy as np
# wing geometry definitions

def wing_area(b, chordR, chordT):
    a1 = b*chordT
    tri1 = (chordR - chordT)*(b/2)/2
    total = a1 + 2*tri1
    return total


def aspect_ratio(b, area):
    AR = (b**2)/area
    return AR


def theta(y,b):
    """
    theta
    calculates the transformed coordinate 0
    y: spanwise coordinate
    b: wing span
    """
    return np.arccos(2*y/b)


def linear_function(y, val_root, val_tip, b):
    """
    linear funciton
    return the local zero lift angle for a wing with linear twist
    y: spanwise coordinate
    val_root: root value
    val_tip: tip value
    b: wing span
    """
    m = (val_root - val_tip)/(b/2)
    val_y = val_root - abs(y)*m
    return val_y

def uniform_grid(b, nPoints):
    half_span = b/2
    delta = half_span/nPoints
    grid = (np.arange(nPoints) + 0.5) * delta
    return grid

# lifting line funcitons
def C_yn(y, n, b, theta, c, m):
    """
    C_yn(y, n, theta, c, m)

    Returns the entries for the coefficient matrix `C` corresponding to a given `y` coordinate and Fourier coefficient `n` (odd coefficients only)
    - `y`: local spanwise coordinate
    - `n`: Fourier coefficient (converted to odd integer internally) 
    - `b`: wing span
    - `theta`: function for coordinate transform w.r.t `y`
    - `c`: function returning local chord w.r.t `y`
    - `m`: function returning local aerofoil lift slope w.r.t `y`
    """
    n_odd = (2*n)-1 # only 
    t1 = (4*b)/(m(y)*c(y))
    t2 = n_odd/np.sin(theta(y))
    t3 = t1 + t2
    return t3*np.sin(n_odd*theta(y))



def build_linear_system(y, b, theta, c, m, a, a0):
    """
    build_linear_system(y, b, theta, c, m, a, a0)

    Builds the linear system for solving Prandtl's Lifting Line Theory on a rectangular wing (can include wing taper, but no sweep). Returns a matrix of coefficients and a vector for the right-hand-side of the linear system. Input:
        - `y`: Vector of chord locations
        - `b`: Wing span 
        - `theta`: function for coordinate transform w.r.t `y`
        - `c`: function returning local chord w.r.t `y`
        - `m`: function returning local aerofoil lift slope w.r.t `y`
        - `a`: function returning local angle of attack w.r.t `y`
        - `a0`: function returning local zero-lift angle w.r.t `y`
    """
    nPoints = len(y)
    C = np.zeros((nPoints, nPoints))
    D = np.zeros(nPoints)
    for j in range(nPoints):
        for i in range(nPoints):
            C[i,j] = C_yn(y[i], j, b, theta, c, m)
        D[j] = np.deg2rad(a(y[j])- a0(y[j]))
    return C, D

        

def gamma0(y, A, theta):
    """
    Gamma0(y, A, nPoints::Int, theta)

    Return local circulation at spanwise location `y`.
    - `y`: local spanwise coordinate
    - `A`: wing Fourier coefficients (vector) 
    - `theta`: function to perform coordinate transformation w.r.t. `y`
    """
    nPoints = len(y)
    G = np.zeros(nPoints)
    for i in range(nPoints):
        for n in range(nPoints):
            G[i] += A[n]*np.sin((2*n-1)*theta(y[i]))
    return G

# Aero Coefficients
def delta(A):
    sum = 0.0
    for n in range(1, len(A)):
        sum += (2 * n + 1) * (A[n] / A[0])**2
    return sum

def efficiency_factor(delta):
    return 1/(1 + delta)
def drag_coefficient(Cl, AR, E):
    return Cl**2/(np.pi*E*AR)
def lift_coefficient(A, AR):
    return np.pi*A[1]*AR
