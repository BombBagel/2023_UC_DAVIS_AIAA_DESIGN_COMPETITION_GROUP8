import numpy as np

# Inputs (airplane parameters)
# c_root [ft] chord of the root of the wing
# c_tip [ft] chord of the tip of the wing
# b_w [ft] wing span
# S_w [ft^2] wing area
# L_ht [ft] distance from wing aerodynamic center to horizontal tail aerodynamic center (momment arm)
# L_vt [ft] distance from wing aerodynamic center to vertical tail aerodynamic center (momment arm)

# Outputs
# taper ratio (lam)
# horizontal tail and vertical tail areas (S_ht, S_vt)

def getempennage(c_root,c_tip,b_w,S_w,L_ht,L_vt):
    # typical tail volume coefficients (table 6.4, Raymer)
    t_tail = 0.95 # reduction to vertical tail coefficient b/c end-plate effect and clean-air effect
    c_ht = 0.9*t_tail
    c_vt = 0.08*t_tail

    # rudder and elevator sizing (table 6.5, Raymer)
    #rudder_per = 0.32
    #elevator_per = 0.25

    lam = c_tip/c_root # taper ratio angle of the wing
    cbar_w = (2/3)*c_root*(1+lam+lam**2)/(1+lam)
    S_ht = c_ht*cbar_w*S_w/L_ht # horizontal tail area
    S_vt = c_vt*b_w*S_w/L_vt # vertical tail area
    return lam,S_ht,S_vt,cbar_w

'''# airplane parameters
c_root = 18.39660 # [ft] chord of the root of the wing
c_tip = 3.17520
b_w = 99.27966 # [ft] wing span
S_w = 922.97191 # [ft^2] wing area
L_ht = 46.916 # [ft] distance from wing aerodynamic center to horizontal tail aerodynamic center
L_vt = 39.899 # [ft] distance from wing aerodynamic center to vertical tail aerodynamic center

lam,S_ht,S_vt,cbar_w = getempennage(c_root,c_tip,b_w,S_w,L_ht,L_vt)
print('taper ratio = {}'.format(lam))
print('S_ht = {} [ft^2]'.format(S_ht))
print('S_vt = {} [ft^2]'.format(S_vt))
print('cbar_w = {} [ft]'.format(cbar_w))'''