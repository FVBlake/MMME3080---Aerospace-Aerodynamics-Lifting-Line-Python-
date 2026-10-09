import function_definitions as fd
import matplotlib.pyplot as plt
import numpy as np

# Freestream properties
U = 50

# Define wing properties
span = 7.5
chord_root = 1
chord_tip = 0.85
alpha_root = 5
alpha_tip = 5
alpha0_const = 0
m_const = 6.5

#redefine geo functions
def Theta(y):
    return fd.theta(y,span)

def c(y):
    return fd.linear_function(y, chord_root, chord_tip, span)

def a(y):
    return fd.linear_function(y, alpha_root, alpha_tip, span)

def a0(y):
    return alpha0_const

def m(y):
    return m_const

#define grid
y = fd.uniform_grid(span, 6)

#build matrix system
C, D = fd.build_linear_system(y, span, Theta, c, m, a, a0)

#solve system
A = np.linalg.solve(C, D)

G0 = fd.gamma0(y, A, Theta)
area = fd.wing_area(span, chord_root, chord_tip)
AR = fd.aspect_ratio(span, area)
Delta = fd.delta(A)
E = fd.efficiency_factor(Delta)
Cl = fd.lift_coefficient(A, AR)
Cd = fd.drag_coefficient(Cl, AR, E)

plt.plot(y, (G0/max(G0)))
plt.title("Circulation")
plt.xlabel("Spanwise Coordinate [m]")
plt.ylabel(r"$\Gamma/\Gamma_0$")
plt.grid("on")
plt.show()

plt.plot(
    [-span/2, -span/2, 0, span/2, span/2, 0, -span/2],
    [-chord_tip/2, chord_tip/2, chord_root/2, chord_tip/2, -chord_tip/2, -chord_root/2, -chord_tip/2],
    label="Wing Geometry" 
)
plt.title(f"Cl = {round(Cl, 4)}, Cd = {round(Cd, 4)}")
plt.xlabel("Spanwise coordinate [m]")
plt.ylabel("Streamwise coordinate [m]")
plt.gca().set_aspect("equal")
plt.grid("on")
plt.show()
