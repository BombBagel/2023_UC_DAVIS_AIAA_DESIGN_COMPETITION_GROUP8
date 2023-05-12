import numpy as np
from scipy.optimize import fsolve

# L = length of fuselage
# diam = diameter of fuselage
# X_fcg = forward cg
# X_acg = aft cg
# b = wing span
# Wn = (assumed) nose gear position
# theta_TB = longitudinal tip over angle (must be >15)
# theta_LOF = longitudinal ground clearance angle (must be >15)
# phi = lateral tip over (overturn) angle (must be <55 ideally or <63 at max)
# psi = lateral ground clearange angle (is actually an output, just make sure this is >5)

def solve_equations(L, diam, X_fcg, X_acg, b, Wn, theta_TB, theta_LOF, phi, psi):
    # Function to solve for Wfm, H, and T
    def equations(vars):
        Wfm, H, T = vars
        Y = diam/2
        theta_TB_deg = np.arctan((Wfm - (X_acg - X_fcg)) / H) * (180/np.pi)
        theta_LOF_deg = np.arctan((H+Y) / (L - (X_fcg + Wfm))) * (180/np.pi)
        #phi_deg = np.arctan((H + diam) / ((b/2) - (T/2))) * (180/np.pi)
        A1 = np.arctan((T/2) / ((X_fcg - Wn) + Wfm))
        D = np.sin(A1) * (((X_fcg - Wn) + Wfm))
        psi_deg = np.arctan(H / D) * (180/np.pi)
        eq1 = theta_TB_deg - theta_TB
        eq2 = theta_LOF_deg - theta_LOF
        #eq3 = phi_deg - phi
        eq4 = psi_deg - psi
        return [eq1, eq2, eq4]

    # Initial guess for Wfm, H, and T
    initial_guess = [2, 8, 5]

    # Solve for Wfm, H, and T
    Wfm, H, T  = fsolve(equations, initial_guess)

    return Wfm, H, T
    # Wfm with respect to cg
    # H is height from bottom of landing gear to cg
    # T is width between each landing gear


# Example:
'''L = 83.0 #ft
diam = 9 #ft
X_fcg = 37.72 #ft
X_acg = 38.05 #ft
b = 111.20 #ft

Wn = 5.5 #ft

theta_TB = 16 #deg
theta_LOF = 16
phi = 6
psi = 54

Wfm, H, T = solve_equations(L, diam, X_fcg, X_acg, b, Wn, theta_TB, theta_LOF, phi, psi)

print("W_fm (w.r.t cg):", Wfm)
print("H:", H)
print("T:", T)
print("Ground Clearance Angle (Requires >5 deg):", np.arctan((H + diam) / ((b/2) - (T/2))) * (180/np.pi))'''
