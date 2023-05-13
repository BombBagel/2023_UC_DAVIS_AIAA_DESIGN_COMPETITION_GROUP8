import pandas as pd
import numpy as np
import sys

# W_static is the weight of the airplane in lbf
# V_stall is the stall speed in ft/s
# nose_num and main_num are the number of wheels for each section
# H, x_nose, x_main, FWD_cg, AFT_cg all must be in the same units of distance
# aircraft_type (1 = general aviation and 2 = transport/bomber)
# margin (1 = FAR 25 requirement, 2 = common practice and is stricter than FAR)
def calc_tire(W_static,V_stall,nose_num,main_num,H,x_nose,x_main,FWD_cg,AFT_cg,aircraft_type,margin):
    # gravitational acceleration
    g = 32.17 # ft/s^2
    
    # calculate needed distance values (Raymer Figure 11.6)
    B = x_main-x_nose
    N_a = AFT_cg-x_nose
    N_f = FWD_cg-x_nose # not used
    M_a = x_main-AFT_cg
    M_f = x_main-FWD_cg
    
    # check if values are within load limits for nose and main landing gears
    if M_a/B <= 0.05:
        M_A = 0.05*B
        print("Tire Sizing: M_a/B > 0.05, not satisfied - M_a/B = {} &(distance between AFT cg and main landing gear is too small relative to the distance between nose and main landing gears)".format(M_a/B))
        print('Requirement: M_a > {}'.format(M_A))
        sys.exit()
    elif M_f/B >= 0.20:
        M_F = 0.20*B
        B_F = M_f/0.20
        print("Tire Sizing: M_f/B < 0.20, not satisfied - M_f/B = {} (distance between FWD cg and main landing gear is too large relative to the distance between nose and main landing gears)".format(M_f/B))
        print('Suggestion 1: M_f < {} [ft], distance between FWD_cg and main landing gear must be less than {} [ft] without moving main landing gear'.format(M_F,M_F))
        print('Suggestion 2: B > {} [ft], distance between the nose and main landing gear without moving main gear'.format(B_F))
        sys.exit()
    
    # tire sizing parameters (Raymer table 11.1)
    if aircraft_type == 1: # general aviation
        A_d = 1.51
        B_d = 0.349
        A_w = 0.7150
        B_w = 0.312
    elif aircraft_type == 2: # transport/bomber
        A_d = 1.63
        B_d = 0.315
        A_w = 0.1043
        B_w = 0.480

    # wheel load margins
    if margin == 1: # FAR 25 (7% margin to wheel loads)
        r = 1.07
    elif margin == 2: # Common (25% margin to wheel loads)
        r = 1.25

    # calculate loads (Raymer equations 11.1 - 11.4)
    max_static_load_main = W_static*N_a/B*r
    max_static_load_nose = W_static*M_f/B*r
    min_static_load_nose = W_static*M_a/B*r
    dynamic_braking_load_nose = 0.31*H*W_static/B*r

    # calculate wheel diameter and width (Raymer table 11.1)
    dia_main = A_d*(max_static_load_main/main_num)**B_d # wheel diameter [in]
    width_main = A_w*(max_static_load_main/main_num)**B_w # wheel width [in]
    dia_nose = dia_main*0.80 # from statistical tire sizing (from 60% to 100% of main tires for tricycles)
    width_nose = width_main*0.80

    # check the tire sizing by kinetic energy from braking
    W_landing = 0.9*W_static*r # landing weight is between 80% to 100% takeoff weight
    KE_braking = 0.5*(W_landing/main_num)/g*V_stall**2 # [ft-lb/s]? (Raymer equation 11.7)
    KE_data = pd.read_csv('KE_tiresize.csv',skiprows=[0]) # kinetic energy per braked wheel vs wheel diameter plot (Raymer figure 11.8)
    large_transport = KE_data.iloc[:,:2] # assume it is comparable to a large transport
    #typical_fighter = KE_data.iloc[:,2:4]
    #expensive_complicated_brakes = KE_data.iloc[:,4:6]
    #general_aviation_and_small_gets = KE_data.iloc[:,6:8] # the curve goes back on itself 

    # format the data
    large_transport = large_transport.dropna()
    large_transport = large_transport.sort_values(by=['X'])
    large_transport = large_transport.values
    dia_KE_main = np.interp(KE_braking,large_transport[:,0],large_transport[:,1])*2 # [in] (the x2 is from the plot being about rim wheel diameter)

    # makes the larger size the saved one
    if dia_main >= dia_KE_main:
        dia_main = dia_main
    else:
        dia_main = dia_KE_main
        dia_nose = 0.80*dia_main # from statistical tire sizing (from 60% to 100% of main tires for tricycles)
        
    return dia_main,dia_nose,width_main,width_nose,KE_braking


# example
W_static = 79110 # lbf
V_stall = 238 # ft/s
nose_num = 2
main_num = 4
FWD_cg = 37.72
AFT_cg = 38.05
x_nose = 5.5
x_main = 3.01+FWD_cg # 1.4 ft offset to satisfy the M_f/B < 0.20 requirement
H = 9.27
aircraft_type = 2
margin = 2

dia_main,dia_nose,width_main,width_nose,KE_braking = calc_tire(W_static,V_stall,nose_num,main_num,H,x_nose,x_main,FWD_cg,AFT_cg,aircraft_type,margin)
print('dia_main = {} [in]'.format(dia_main))
print('width_main = {} [in]'.format(width_main))
print('dia_nose = {} [in]'.format(dia_nose))
print('width_nose = {} [in]'.format(width_nose))
print('KE_braking = {} x 10^6 [ft-lbf]'.format(KE_braking/10**6))