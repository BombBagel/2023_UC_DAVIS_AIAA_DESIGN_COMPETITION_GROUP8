import numpy as np

# Function to solve stall constraint
#
#
def stallConstraint(rho,C_L_maxL,graphPoints):
    V_stall = 238/1.3 # Stall speed [ft/s]
    W_S_Stall = (0.5*rho*(V_stall**2)*C_L_maxL)*np.ones(graphPoints)
    return W_S_Stall

# Function to solve for takeoff constraint
#
#
def takeoffConstraint(C_D_o_ToGDo,rho_SL,N_Engine,D_P,W_SArray,C_L_maxTO):
    P_TO = 5000 # Takeoff power [hp]
    BFL = 4500 # Balanced field length
    S_TOG = BFL*0.6 # Takeoff groundrun
    eta_p = 0.8 # Propulsive Efficiency
    mu_G = 0.025 # Roskam Table 3.2 ground coefficient
    l_p = 5.75 # Roskam pg 101
    k_1 = 0.0376 # Roskam pg 101
    C_D_o = C_D_o_ToGDo
    k_2 = l_p*((rho_SL/rho_SL)*(N_Engine*D_P**2/P_TO))**(1/3)
    W_P_Takeoff_Main = ((1/(k_2))*(((k_1*W_SArray/(S_TOG*rho_SL) + 0.72*C_D_o)/C_L_maxTO) + mu_G))**(-1)
    return W_P_Takeoff_Main

# Function to solve for landing constraint
#
#
def landingConstraint(C_L_maxL,rho,rho_SL,WL_WTO,graphPoints):
    S_FL = 4500 # Landing field length [ft]
    S_a = 1000 # Obstacle clearing distance [ft]
    W_SLanding = ((S_FL - S_a)/80)*(C_L_maxL*(rho/rho_SL) / WL_WTO)*np.ones(graphPoints)
    return W_SLanding

# Function to solve climb constraints
#
#
def climbConstraints(C_D_o_ToGUp,C_L_maxTO,W_SArray,rho_28000,e_Takeoff,AR,N_Engine,eta_p,graphPoints,C_D_o_ToGDo,C_D_o_Clean,C_L_maxCr,
                     e_Clean,C_D_o_LaGDo,C_L_maxL,e_Landing,rho,WL_WTO,Single_Run = False,**kwargs):
    file2 = kwargs.get('file', None)
    # Takeoff climb
    k_s = 1.2
    G = 0.012
    C_D_o = C_D_o_ToGUp
    C_L_maxCL = C_L_maxTO
    V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
    W_P_TakeoffClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Takeoff*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)) 
                              / (550*eta_p/V))**(-1))*np.ones(graphPoints)
    # Transition climb
    k_s = 1.19
    G = 0
    C_D_o = C_D_o_ToGDo
    C_L_maxCL = C_L_maxTO
    V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
    W_P_TransitionClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Takeoff*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)) 
                                 / (550*eta_p/V))**(-1))*np.ones(graphPoints)
    # Second segment climb
    k_s = 1.2
    G = 0.024
    C_D_o = C_D_o_ToGUp
    C_L_maxCL = C_L_maxTO
    V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
    W_P_SSClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Takeoff*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)) 
                         / (550*eta_p/V))**(-1))*np.ones(graphPoints)
    # En-route climb
    k_s = 1.25
    G = 0.012
    C_D_o = C_D_o_Clean
    C_L_maxCL = C_L_maxCr
    V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
    W_P_EnRouteClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Clean*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)/0.94) 
                              / (550*eta_p/V))**(-1))*np.ones(graphPoints)
    # Balked landing climb (AEO)
    k_s = 1.3
    G = 0.032
    C_D_o = C_D_o_LaGDo
    C_L_maxCL = C_L_maxL
    V = np.sqrt(2*W_SArray/(rho*C_L_maxCL))
    W_P_BalkedAEO_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Landing*AR*k_s**2)) + G)*(WL_WTO/0.8) 
                           / (550*eta_p/V))**(-1))*np.ones(graphPoints)
    # Balked landing climb (OEI)
    k_s = 1.5
    G = 0.021
    C_D_o = (C_D_o_ToGDo + C_D_o_LaGDo)/2
    if Single_Run == True:
        print('Average parasite drag for landing and takeoff gear down:', C_D_o, file = file2)
    C_L_maxCL = C_L_maxL
    V = np.sqrt(2*W_SArray/(rho*C_L_maxCL))
    W_P_BalkedOEI_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Landing*AR*k_s**2)) + G)*(N_Engine*WL_WTO/0.8/(N_Engine-1)) 
                           / (550*eta_p/V))**(-1))*np.ones(graphPoints)
    return W_P_TakeoffClimb_Main,W_P_TransitionClimb_Main,W_P_SSClimb_Main,W_P_EnRouteClimb_Main,W_P_BalkedAEO_Main,W_P_BalkedOEI_Main

# Function to solve for cruise constraint
#
#
def cruiseConstraint(V_cruise,C_D_o_Clean,rho_28000,W_SArray,WCruise_WTakeoff,AR,e_Clean,eta_p):
    C_D_o = C_D_o_Clean # Minimum drag from drag polar for clean configuration
    q = 0.5*V_cruise**2*rho_28000
    W_P_Cruise_Main = ((q*V_cruise*(C_D_o + (((W_SArray**2)*(WCruise_WTakeoff**2)) / ((q**2)*np.pi*AR*e_Clean)))) 
                       / (550*eta_p*W_SArray) * (1.125))**(-1)
    return W_P_Cruise_Main

# Function to solve ceiling constraint
#
#
def ceilingConstraint(C_D_o_Clean,W_SArray,rho_28000,e_Clean,AR,eta_p,graphPoints):
    C_D_o = C_D_o_Clean
    CL_maxLD = 0.75
    V = np.sqrt(2*W_SArray/(rho_28000*CL_maxLD))
    W_P_Ceiling_Main = ((2*np.sqrt(C_D_o/(np.pi*e_Clean*AR))*(V/(550*eta_p)))**(-1))*np.ones(graphPoints)
    return W_P_Ceiling_Main
