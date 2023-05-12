import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from powertrain_matrix import getpowertrain
from matplotlib.backends.backend_pdf import PdfPages
import os
from fpdf import FPDF
from pypdf import PdfMerger


# Function for estimating weight
## 
###
def Calc_weight(i, pilots, attendants, passengers, R, L_D, E, pe, v, D_P, N_Prop, N_Engine):
    C = (0.5*v) / (550*pe*3600) # Specific fuel consumption in lb/(lbf hr)

    # Establish Raymer's method coefficients
    A = df1["Raymer \nCoefficient A"][i]
    B = df1["Raymer \nCoefficient B"][i]

    # Calculate weight of crew and payload
    W_crew = (190+30)*(pilots+attendants)
    W_passengers = (200+40)*(passengers)

    # Establish fuel fractions at different flight segments
    W1_W0 = 0.970 # Takeoff fuel fraction
    W2_W1 = 0.985 # Climb fuel fraction
    W3_W2 = np.exp((-R*C) / (v*L_D)) # Cruise fuel fraction
    W4_W3 = 0.990 # Descent fuel fraction
    W5_W4 = np.exp((-E*C) / (L_D)) # Loiter fuel fraction
    W6_W5 = 0.995 # Landing fuel fraction

    # Calculate fuel fraction
    Wf_W0 = (1 - W1_W0*W2_W1*W3_W2*W4_W3*W5_W4*W6_W5)*1.06 # Multiplied by 1.06 to account for trapped and reserved fuel

    # Set tolerance
    sigma = 10**-6
    delta = 2*sigma

    # Guess MTOW
    W_guess = 60000
    W0 = W_guess

    while delta > sigma:
        # Calculate empty weight fraction 
        WE_W0 = A*W0**B
        W0_new = (W_crew + W_passengers) / (1-Wf_W0-WE_W0) 
        delta = abs(W0_new-W0) / abs(W0_new)
        W0 = W0_new
        #print(W0)

    M_Batt = (5443.1084*W0*Wf_W0*0.45*0.2*0.28/453) # Estimated mass of battery
    W_batt_W0 = M_Batt/(W0) # Battery mass fraction
    R = R*0.8 # Reduce fuel range to account for battery taking over during cruise
    
    # Calculate landing fuel fraction
    W1_W0 = 0.970 # Takeoff fuel fraction
    W2_W1 = 0.985 # Climb fuel fraction
    W3_W2 = 1 # No cruise segment
    ELanding = 0.25*3600 # 15 minute loiter time
    W4_W3 = 0.990 # Descent fuel fraction
    W5_W4 = np.exp((-ELanding*C) / (L_D)) # Loiter fuel fraction
    WL_WTO = (W1_W0*W2_W1*W3_W2*W4_W3*W5_W4)*1.06 # Landing weight fraction
    print('Landing Weight Ratio is:', WL_WTO, file = file1)

    # Recalculate fuel fraction since range has changed
    W1_W0 = 0.970 # Takeoff fuel fraction
    W2_W1 = 0.985 # Climb fuel fraction
    W3_W2 = np.exp((-R*C) / (v*L_D)) # Cruise fuel fraction
    W4_W3 = 0.990 # Descent fuel fraction
    W5_W4 = np.exp((-E*C) / (L_D)) # Loiter fuel fraction
    W6_W5 = 0.995 # Landing fuel fraction

    # Calculate fuel fraction
    Wf_W0 = (1 - W1_W0*W2_W1*W3_W2*W4_W3*W5_W4*W6_W5)*1.06 # Multiplied by 1.06 to account for trapped and reserved fuel
    WCruise_WTakeoff = W1_W0*W2_W1*W3_W2*1.06 # Cruise weight fraction
    print('Cruise Weight Ratio is:', WCruise_WTakeoff, file = file1)
    
    # Set tolerance
    sigma = 10**-6
    delta = 2*sigma

    # Guess MTOW
    W_guess = 30000
    W0 = W_guess

    while delta > sigma:
        # Calculate empty weight fraction 
        WE_W0 = A*W0**B
        W0_new = (W_crew + W_passengers) / (1-Wf_W0-WE_W0-W_batt_W0) 
        delta = abs(W0_new-W0) / abs(W0_new)
        W0 = W0_new
    print('Max weight of hybrid electric concept', (i+1), 'is:', W0, 'lbs', file = file1)
    Unit_Cost = Calc_cost(i, W0*WE_W0, D_P, N_Prop, N_Engine)/50
    print('Total cost of hybrid electric concept', (i+1), 'is:', Unit_Cost, '$', file = file1)
    return W0, WL_WTO, WCruise_WTakeoff

# Function for finding cost of aircraft
## Takes in the row number, i and WE which is the empty weight
### Returns the total cost and prints out the individual production costs
def Calc_cost(i, WE, D_P, N_Prop, N_Engine):
    W_Airframe = WE*0.45 # Airframe weight (lbs) (We assume is 45% of empty weight)
    V_H = df1["Speed (knots)"][i] # Maximum Level Airspeed (KTAS)
    Q = 50 # Number of aircraft to be produced over a 5 year period
    Q_M = 5 # Number of aircraft produced in 1 month
    Q_Proto = 2 # Number of prototype aircraft to be produced 
    CPI = 1.31 # Consumer price index for consideration of inflation since 2012
    # Cost of engineering
    R_Eng = 92 # Rate of engineering labor in US$ per hour (92$/h is recommended) 
    C_Eng = 0.083 * (W_Airframe**0.791) * (V_H**1.521) * (Q**0.183) * 0.67 * 2.00 * 1.03 * 1.03 * 1.66 * R_Eng * CPI
    print('Cost of engineering for concept', i+1, 'is', C_Eng/Q, '$', file = file1)
    # Cost of tooling
    R_Tool = 61 # Rate of tooling labor in US$ per hour (61$/h is recommended) 
    C_Tool = 2.1036 * (W_Airframe**0.764) * (V_H**0.899) * (Q**0.178) * (Q_M**0.066) * 2.00 * 0.95 * 1.02 * 1.01 * 1.10 * R_Tool * CPI
    print('Cost of tooling for concept', i+1, 'is', C_Tool/Q, '$', file = file1)
    # Cost of manufacturing
    R_Mfg = 53 # Rate of manufacturing labor in US$ per hour (53$/h is recommended) 
    C_Mfg = 20.2588 * (W_Airframe**0.74) * (V_H**0.543) * (Q**0.524) * 0.75 * 1.25 * 1.01 * 1.10 * R_Mfg * CPI
    print('Cost of manufacturing for concept', i+1, 'is', C_Mfg/Q, '$', file = file1)
    # Cost of development support
    C_Dev = 0.06458 * (W_Airframe**0.873) * (V_H**1.89) * (Q_Proto**0.346) * CPI * 0.5 * 1.5 * 1.01 * 1.03 * 1.05
    print('Cost of development for concept', i+1, 'is', C_Dev/Q, '$', file = file1)
    # Cost of flight test operations
    C_Ft = 0.009646 * (W_Airframe**1.16) * (V_H**1.3718) * (Q_Proto**1.281) * CPI * 0.5 * 1.5
    print('Cost of flight test operations for concept', i+1, 'is', C_Ft/Q, '$', file = file1)
    # Cost of quality control
    C_Qc = 0.13 * (C_Mfg) * 0.5 * 1.5 * 1.5
    print('Cost of quality control for concept', i+1, 'is', C_Qc/Q, '$', file = file1)
    # Cost of materials
    C_Mat = 24.896 * (W_Airframe**0.689) * (V_H**0.624) * (Q*0.792) * CPI * 0.75 * 1.02 * 1.01 * 1.05 
    print('Cost of materials for concept', i+1, 'is', C_Mat/Q, '$', file = file1)

    # Cost of electric motors
    P_Em = 2500 # Rated power P (hp)
    C_Em = 174 * N_Engine * P_Em * CPI
    print('Cost of electric motors for concept', i+1, 'is', C_Em/Q, '$', file = file1)
    # Cost of power management system 
    P_EmTotal = P_Em * N_Engine
    C_Pms = 150 * P_EmTotal * CPI
    print('Cost of power management system for concept', i+1, 'is', C_Pms/Q, '$', file = file1)
    # Cost of battery (assuming 200$/kWh)
    E_Bat = 918.400 # Total energy in battery (kWh)
    C_Bat = 200 * E_Bat * CPI
    print('Cost of battery for concept', i+1, 'is', C_Bat/Q, '$', file = file1)
    # Cost of Propeller
    P_SHP = 2500 # Shaft horsepower
    C_Prop = 210 * N_Prop * CPI * D_P**2 * ((P_SHP/D_P)**0.12)
    print('Cost of propeller for concept', i+1, 'is', C_Prop/Q, '$', file = file1)

    TotalCost = (C_Eng + C_Tool + C_Mfg + C_Dev + C_Ft + C_Mat + C_Em + C_Pms + C_Bat + C_Prop)*1.15
    return TotalCost

# Matrix solving function
##
###
def ComponentPowerLoading(W_P, ComponentMatrix):
    P_W = W_P**-1
    PowerLoadingMatrix =[]
    for i in range(len(P_W)):
        b = np.zeros([len(ComponentMatrix),1])
        b[len(b)-1] = P_W[i]
        x = np.linalg.solve(ComponentMatrix,b)
        PowerLoadingMatrix = np.append(PowerLoadingMatrix, x)
    PowerLoadingMatrix = np.reshape(PowerLoadingMatrix,[len(ComponentMatrix),len(P_W)], order='F')
    PowerLoadingMatrix = PowerLoadingMatrix**-1
    return PowerLoadingMatrix

# Function to solve for all constraints
##
###
def ConstraintSolve():
    global W_P_Takeoff_Matrix
    global W_P_TakeoffClimb_Matrix
    global W_P_TransitionClimb_Matrix
    global W_P_SSClimb_Matrix
    global W_P_EnRouteClimb_Matrix
    global W_P_BalkedOEI_Matrix
    global W_P_BalkedAEO_Matrix
    global W_P_Cruise_Matrix
    global W_P_Ceiling_Matrix
    W_P_Takeoff_Matrix = ComponentPowerLoading(W_P_Takeoff_Main, A)
    W_P_TakeoffClimb_Matrix = ComponentPowerLoading(W_P_TakeoffClimb_Main, A)
    W_P_TransitionClimb_Matrix = ComponentPowerLoading(W_P_TransitionClimb_Main, A)
    W_P_SSClimb_Matrix = ComponentPowerLoading(W_P_SSClimb_Main, A)
    W_P_EnRouteClimb_Matrix = ComponentPowerLoading(W_P_EnRouteClimb_Main, A)
    W_P_BalkedOEI_Matrix = ComponentPowerLoading(W_P_BalkedOEI_Main, A)
    W_P_BalkedAEO_Matrix = ComponentPowerLoading(W_P_BalkedAEO_Main, A)
    W_P_Cruise_Matrix = ComponentPowerLoading(W_P_Cruise_Main, B)
    W_P_Ceiling_Matrix = ComponentPowerLoading(W_P_Ceiling_Main, B)

# Plotting function
##
###
def PlotPowerLoadingMatrix(n): 
    global W_P_Takeoff 
    global W_P_TakeoffClimb 
    global W_P_TransitionClimb 
    global W_P_SSClimb
    global W_P_EnRouteClimb 
    global W_P_BalkedAEO 
    global W_P_BalkedOEI 
    global W_P_Cruise 
    global W_P_Ceiling 
    W_P_Takeoff = W_P_Takeoff_Matrix[n,:]
    W_P_TakeoffClimb = W_P_TakeoffClimb_Matrix[n,:]
    W_P_TransitionClimb = W_P_TransitionClimb_Matrix[n,:]
    W_P_SSClimb = W_P_SSClimb_Matrix[n,:]
    W_P_EnRouteClimb = W_P_EnRouteClimb_Matrix[n,:]
    W_P_BalkedAEO = W_P_BalkedAEO_Matrix[n,:]
    W_P_BalkedOEI = W_P_BalkedOEI_Matrix[n,:]
    W_P_Cruise = W_P_Cruise_Matrix[n,:]
    W_P_Ceiling = W_P_Ceiling_Matrix[n,:]
    n = n + 1
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_S_Stall, W_PArray, label = 'Stall', color = 'C8')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_Takeoff, label = 'Takeoff', color = 'C1')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SLanding, W_PArray, label = 'Landing', color = 'C0')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_TakeoffClimb, label = 'Takeoff climb', color = 'red')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_TransitionClimb, label = 'Transition climb', color = 'brown', linestyle = 'dashed')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_SSClimb, label = 'Second segment climb', color = 'purple')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_EnRouteClimb, label = 'En-route climb', color = 'darkgreen', linestyle = 'dotted')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_BalkedAEO, label = 'Balked landing climb (AEO)', color = 'black')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_BalkedOEI, label = 'Balked landing climb (OEI)', color = 'C2')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_Cruise, label = 'Cruise', color = 'cyan')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(W_SArray, W_P_Ceiling, label = 'Ceiling', color = 'grey')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_ylim([0,35])
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_xlim([0,100])
    DesignW_S, DesignW_P = ShadeSubplots(n)
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(DesignW_S, DesignW_P, 'ro', markersize = '4')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_title('Graph of Power Loading against Wing Loading for ' + Components[i])
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_xlabel(r'$W_{0}/S_{ref}$' '[lb/ft^2]')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_ylabel(r'$W_{0}/P_{0}$' '[hp/ft^2]')
    #axs[int(np.ceil(n/3)-1),int(n-3*(np.floor(n/3))-1)]..ylabel(r"$\frac{{W}}{{{}}}$".format(Components[i]) + "[lb/hp]")
    if Components[i] == 'P_e1':
        global EM_W_P
        EM_W_P = DesignW_P

# Shading function for main graph
## Also returns optimal point
###
def Shade():
    min = np.minimum(W_P_BalkedOEI_Main, W_P_Cruise_Main)
    min = np.minimum(min, W_P_Takeoff_Main)
    min = np.minimum(min, W_P_Ceiling_Main)
    min = np.minimum(min, W_P_BalkedAEO_Main)
    min = np.minimum(min, W_P_EnRouteClimb_Main)
    min = np.minimum(min, W_P_SSClimb_Main)
    min = np.minimum(min, W_P_TakeoffClimb_Main)
    min = np.minimum(min, W_P_TransitionClimb_Main)
    plt.fill_between(W_SArray, min, 0, color='grey')
    plt.fill_betweenx(W_PArray, np.minimum(W_SLanding, W_S_Stall), 100*np.ones(graphPoints), color = 'white')
    DesignW_S, DesignW_P = PlotPoint(min)
    return DesignW_S, DesignW_P

# Shading function for subplots 
##
###
def ShadeSubplots(n):
    min = np.minimum(W_P_BalkedOEI, W_P_Cruise)
    min = np.minimum(min, W_P_Takeoff)
    min = np.minimum(min, W_P_Ceiling)
    min = np.minimum(min, W_P_BalkedAEO)
    min = np.minimum(min, W_P_EnRouteClimb)
    min = np.minimum(min, W_P_SSClimb)
    min = np.minimum(min, W_P_TakeoffClimb)
    min = np.minimum(min, W_P_TransitionClimb)
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].fill_between(W_SArray, min, 0, color='grey')
    axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].fill_betweenx(W_PArray, np.minimum(W_SLanding, W_S_Stall), 100*np.ones(graphPoints), color = 'white')
    DesignW_S, DesignW_P = PlotPoint(min)
    return DesignW_S, DesignW_P

# Function to find best point in shaded area
##
###
def PlotPoint(min):
    DesignW_S = np.floor(np.minimum(W_SLanding, W_S_Stall)[1])
    DesignW_P = np.floor(np.interp(DesignW_S, W_SArray, min))
    return DesignW_S, DesignW_P

# Function to initialize subplot
##
### 
def InitializeSubplots(Components):
    n = len(Components)
    global fig
    global axs 
    fig, axs = plt.subplots(int(np.ceil(n/2)), 2, figsize=(15, 20))
    return fig



















# Read excel file 
df1 = pd.read_excel('Aircraft Data.xlsx')
print(df1)

# initialize number of points on graph
graphPoints = 400

# Create txt file 
file1 = open('text.txt', 'w')
file2 = open('text1.txt', 'w')

j = 0
# Payload for aircraft
pilots = df1["Pilots"][j]
passengers = df1["Passengers"][j]
# Attendants will be expressed as a function of crew since we need 1 for every 50 passengers
attendants = np.ceil(passengers/50)
# Establish Flight Profile
R = df1["Range \n(nmi)"][j]*5280 # Range in ft
L_D_max = df1["L/D"][j] # Max lift to drag ratio
L_D = L_D_max*0.94 # Actual lift to drag ratio 
E = df1["Endurance \n(hours)"][j]*3600 # Endurance in seconds
pe = df1["Propellor \nEfficiency"][j] # Propeller efficiency
v = df1["Speed (knots)"][j]*1.668 # Speed in ft/s 
# Extract airplane parameters
N_Prop = df1["No. Propellers"][j] # Number of propellers
D_P = df1["Propeller \nDiameter (ft)"][j] # Diameter of propeller in ft
N_Engine = df1["No. Engines"][j] # Number of required engines
AR = df1["AR"][j] # Aspect ratio
C_L_maxL = df1["Clmax \nLanding"][j] # Max Cl of lift 
C_L_maxCr = df1["Clmax \nCruise"][j] # Max Cl for cruise
C_L_maxTO = df1["Clmax \nTakeoff"][j] # Max Cl for takeoff

W_0, WL_WTO, WCruise_WTakeoff = Calc_weight(j, pilots, attendants, passengers, R, L_D, E, pe, v, D_P, N_Prop, N_Engine)

#
# Drag polar calcs
#
pp = PdfPages('testin.pdf')
pp1 = PdfPages('testin1.pdf')
fig1 = plt.figure()

# Estimate wetted area from historical data
c = -0.0866 # Regression constant from Table 3.5 Roskam vol 1 for twin propeller plane
d = 0.8099 # Regression constant from Table 3.5 Roskam vol 1 for twin propeller plane
S_wet = (10**c)*(W_0)**d # Wetted area
y = 10**(0.31739*np.log10(W_0) + 0.26161) # Equation estimated from Torenbeek, 1990, Fig. 7.4, takeoff weight divided by wing loading
S = W_0/y # Value estimated from Torenbeek, 1990, Fig. 7.4, takeoff weight divided by wing loading
C_f = 0.0026 # Value obtained from Raymer, Table 12.3
C_Do = C_f*S_wet/S

deltaTakeoffFlaps = df1['delta\nTakeoff\nFlaps'][j]
deltaLandingFlaps = df1['delta\nLanding\nFlaps'][j]
deltaLandingGear = df1['delta\nLanding\nGear'][j]
# Clean configuration
C_L_Clean = np.linspace(-1.1, 1.1, 50)
e_Clean = df1["e Clean"][j] 
C_D_Clean = C_Do + C_L_Clean**2/(np.pi*AR*e_Clean) 
C_D_o_Clean = C_Do # Clean config parasite drag
print('Clean config drag polar is:', round(C_D_o_Clean, 4), '+', round(1/(np.pi*AR*e_Clean), 4), 'C_L^2', file = file2)
plt.plot(C_D_Clean, C_L_Clean)

# Takeoff configuration 
C_L_Takeoff = np.linspace(-2.0, 2, 50)
e_Takeoff = df1["e Takeoff"][j] 
C_D_Takeoff = (C_Do+deltaTakeoffFlaps) + C_L_Takeoff**2/(np.pi*AR*e_Takeoff) # Takeoff flaps, gear up, add 0.02 to parasite drag coefficient, Roskam vol 1, Table 3.6
C_D_o_ToGUp = C_Do + deltaTakeoffFlaps # Takeoff config gear up
print('Takeoff config gear up drag polar is:', round(C_D_o_ToGUp, 4), '+', round(1/(np.pi*AR*e_Takeoff), 4), 'C_L^2', file = file2)
plt.plot(C_D_Takeoff, C_L_Takeoff)
C_D_Takeoff = (C_Do+deltaTakeoffFlaps+deltaLandingGear) + C_L_Takeoff**2/(np.pi*AR*e_Takeoff) # Takeoff flaps, gear down, add 0.02+0.025 to parasite drag coefficient, Roskam vol 1, Table 3.6
C_D_o_ToGDo = C_Do + deltaTakeoffFlaps + deltaLandingGear # Takeoff config gear down
print('Takeoff config gear down drag polar is:', round(C_D_o_ToGDo, 4), '+', round(1/(np.pi*AR*e_Takeoff), 4), 'C_L^2', file = file2)
plt.plot(C_D_Takeoff, C_L_Takeoff)

# Landing configuration
C_L_Landing = np.linspace(-2.5, 2.5, 50)
e_Landing = df1["e Landing"][j]
C_D_Landing = (C_Do + deltaLandingFlaps) + C_L_Landing**2/(np.pi*AR*e_Landing) # Landing flaps, gear up, add 0.075 to parasite drag coefficient, Roskam vol 1, Table 3.6
C_D_o_LaGUp = C_Do + deltaLandingFlaps # Landing config gear up
print('Landing config gear up drag polar is:', round(C_D_o_LaGUp, 4), '+', round(1/(np.pi*AR*e_Landing), 4), 'C_L^2', file = file2)
plt.plot(C_D_Landing, C_L_Landing)
C_D_Landing = (C_Do + deltaLandingFlaps + deltaLandingGear) + C_L_Landing**2/(np.pi*AR*e_Landing) # Landing flaps, gear down, add 0.075+0.025 to parasite drag coefficient, Roskam vol 1, Table 3.6
C_D_o_LaGDo = C_Do + deltaLandingFlaps + deltaLandingGear # Landing config gear down
print('Landing config gear down drag polar is:', round(C_D_o_LaGDo, 4), '+', round(1/(np.pi*AR*e_Landing), 4), 'C_L^2', file = file2)
plt.plot(C_D_Landing, C_L_Landing)

plt.legend(['Clean', 'Takeoff flaps, gear up', 'Takeoff flaps, gear down', 'Landing flaps, gear up', 'Landing flaps, gear down'])
plt.xlabel(r'$C_{D}$')
plt.ylabel(r'$C_{L}$')
plt.title('Drag Polar Graphs')
plt.savefig('Drag Polar1.png')
pp.savefig(fig1, bbox_inches="tight")

#
# Plotting W/P to W/S graph section
#
fig2 = plt.figure(figsize=(12, 10.5))
rho = 0.002048 # Density at 5000 ft
rho_SL = 0.00237 # Density at sea level
rho_28000 = (0.0009567 + 0.002048)/2 # Density at FL280 !!!!!!!MIDPOINT NOW
W_SArray = np.linspace(0.001, 100, graphPoints)
W_PArray = np.linspace(0, 35, graphPoints)

# Stall constraint
V_stall = 238/1.3 # Stall speed [ft/s]
W_S_Stall = (0.5*rho*(V_stall**2)*C_L_maxL)*np.ones(graphPoints)

# Takeoff constraint
P_TO = 5000 # Takeoff power [hp]
BFL = 4500 # Balanced field length
S_TOG = BFL*0.6 # Takeoff groundrun
eta_p = 0.8 # Propulsive Efficiency
mu_G = 0.025 # Roskam Table 3.2 ground coefficient
l_p = 5.75 # Roskam pg 101
k_1 = 0.0376 # Roskam pg 101
C_D_o = C_D_o_ToGDo
k_2 = l_p*((rho/rho_SL)*(N_Engine*D_P**2/P_TO))**(1/3)
W_P_Takeoff_Main = ((1/(k_2))*(((k_1*W_SArray/(S_TOG*rho) + 0.72*C_D_o)/C_L_maxTO) + mu_G))**(-1)

# Landing constraint
S_FL = 4500 # Landing field length [ft]
S_a = 1000 # Obstacle clearing distance [ft]
W_SLanding = ((S_FL - S_a)/80)*(C_L_maxL*(rho/rho_SL) / WL_WTO)*np.ones(graphPoints)

# Climb constraints
#
# Takeoff climb
k_s = 1.2
G = 0.012
C_D_o = C_D_o_ToGUp
C_L_maxCL = C_L_maxTO
V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
W_P_TakeoffClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Takeoff*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)) / (550*eta_p/V))**(-1))*np.ones(graphPoints)
# Transition climb
k_s = 1.19
G = 0
C_D_o = C_D_o_ToGDo
C_L_maxCL = C_L_maxTO
V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
W_P_TransitionClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Takeoff*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)) / (550*eta_p/V))**(-1))*np.ones(graphPoints)
# Second segment climb
k_s = 1.2
G = 0.024
C_D_o = C_D_o_ToGUp
C_L_maxCL = C_L_maxTO
V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
W_P_SSClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Takeoff*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)) / (550*eta_p/V))**(-1))*np.ones(graphPoints)
# En-route climb
k_s = 1.25
G = 0.012
C_D_o = C_D_o_Clean
C_L_maxCL = C_L_maxCr
V = np.sqrt(2*W_SArray/(rho_28000*C_L_maxCL))
W_P_EnRouteClimb_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Clean*AR*k_s**2)) + G)*(N_Engine/0.8/(N_Engine-1)/0.94) / (550*eta_p/V))**(-1))*np.ones(graphPoints)
# Balked landing climb (AEO)
k_s = 1.3
G = 0.032
C_D_o = C_D_o_LaGDo
C_L_maxCL = C_L_maxL
V = np.sqrt(2*W_SArray/(rho*C_L_maxCL))
W_P_BalkedAEO_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Landing*AR*k_s**2)) + G)*(WL_WTO/0.8) / (550*eta_p/V))**(-1))*np.ones(graphPoints)
# Balked landing climb (OEI)
k_s = 1.5
G = 0.021
C_D_o = (C_D_o_ToGDo + C_D_o_LaGDo)/2
print('Average parasite drag for landing and takeoff gear down:', C_D_o, file = file2)
C_L_maxCL = C_L_maxL
V = np.sqrt(2*W_SArray/(rho*C_L_maxCL))
W_P_BalkedOEI_Main = ((((k_s**2*C_D_o/C_L_maxCL) + (C_L_maxCL/(np.pi*e_Landing*AR*k_s**2)) + G)*(N_Engine*WL_WTO/0.8/(N_Engine-1)) / (550*eta_p/V))**(-1))*np.ones(graphPoints)

# Cruise constraint
V_cruise = 464 # Cruise speed [ft/s]
C_D_o = C_D_o_Clean # Minimum drag from drag polar for clean configuration
q = 0.5*V_cruise**2*rho_28000
W_P_Cruise_Main = ((q*V_cruise*(C_D_o + (((W_SArray**2)*(WCruise_WTakeoff**2)) / ((q**2)*np.pi*AR*e_Clean)))) / (550*eta_p*W_SArray) * (1.125))**(-1)

# Ceiling constraint
C_D_o = C_D_o_Clean
CL_maxLD = 0.75
V = np.sqrt(2*W_SArray/(rho_28000*CL_maxLD))
W_P_Ceiling_Main = ((2*np.sqrt(C_D_o/(np.pi*e_Clean*AR))*(V/(550*eta_p)))**(-1))*np.ones(graphPoints)

plt.plot(W_S_Stall, W_PArray, label = 'Stall', color = 'C8')
plt.plot(W_SArray, W_P_Takeoff_Main, label = 'Takeoff', color = 'C1')
plt.plot(W_SLanding, W_PArray, label = 'Landing', color = 'C0')
plt.plot(W_SArray, W_P_TakeoffClimb_Main, label = 'Takeoff climb', color = 'red')
plt.plot(W_SArray, W_P_TransitionClimb_Main, label = 'Transition climb', color = 'brown', linestyle = 'dashed')
plt.plot(W_SArray, W_P_SSClimb_Main, label = 'Second segment climb', color = 'purple')
plt.plot(W_SArray, W_P_EnRouteClimb_Main, label = 'En-route climb', color = 'darkgreen', linestyle = 'dotted')
plt.plot(W_SArray, W_P_BalkedAEO_Main, label = 'Balked landing climb (AEO)', color = 'black')
plt.plot(W_SArray, W_P_BalkedOEI_Main, label = 'Balked landing climb (OEI)', color = 'C2')
plt.plot(W_SArray, W_P_Cruise_Main, label = 'Cruise', color = 'cyan')
plt.plot(W_SArray, W_P_Ceiling_Main, label = 'Ceiling', color = 'grey')

DesignW_S, DesignW_P = Shade()

plt.plot(DesignW_S, DesignW_P, 'ro', markersize = '4')

plt.legend(fontsize='8', loc='center left')
plt.xlabel(r'$\frac{W}{S}$' '[lb/ft^2]')
plt.ylabel(r'$\frac{W}{P}$' ' [lb/hp]')
plt.ylim([0,25])
plt.xlim([0, 100])
plt.title('Graph of Power Loading against Wing Loading')
plt.savefig('Power Loading against Wing Loading.png',bbox_inches="tight")
pp1.savefig(fig2, bbox_inches="tight")

# Let A = matrix for battery discharging and B = matrix for battery charging 

# get matrix for serial powertrain
architecture,state='serial','propulsion'
A,Components = getpowertrain(architecture,state)
architecture,state='serial','charge-both-propulsion'
B,Components = getpowertrain(architecture,state)

# Solve for constraints
ConstraintSolve()
fig3 = InitializeSubplots(Components)
# Plot constraints for each component
for i in range(len(Components)):
    PlotPowerLoadingMatrix(i)
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines[0:11], labels[0:11], loc = 'lower center', ncol = 5, bbox_to_anchor=(0.5,0.05), bbox_transform=fig.transFigure, fontsize="large")
plt.savefig('Component Power Loadings1A.png',bbox_inches="tight")
pp1.savefig(fig3, bbox_inches="tight")

# get matrix for parallel powertrain
architecture,state='parallel','propulsion'
A,Components = getpowertrain(architecture,state)
architecture,state='parallel','charge-both-propulsion'
B,Components = getpowertrain(architecture,state)

# Solve for constraints
ConstraintSolve()
fig4 = InitializeSubplots(Components)
# Plot constraints for each component
for i in range(len(Components)):
    PlotPowerLoadingMatrix(i)
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines[0:11], labels[0:11], loc = 'lower center', ncol = 5, bbox_to_anchor=(0.5,0.05), bbox_transform=fig.transFigure, fontsize="large")
plt.savefig('Component Power Loadings2A.png',bbox_inches="tight")
pp1.savefig(fig4, bbox_inches="tight")
pp.close()
pp1.close()

#
# Updated plane weight estimate
#
    
e_f = 12000*3600 # Specific energy of fuel in J/kg
e_bat = 750*3600 # Specific energy of battery in J/kg
    
# Attendants will be expressed as a function of crew since we need 1 for every 50 passengers
attendants = np.ceil(passengers/50)
W_crew = (190+30)*(pilots+attendants)
W_passengers = (200+40)*(passengers)
W_Payload = W_crew + W_passengers

# Set tolerance
sigma = 10**-6
delta = 2*sigma

while delta > sigma:
    #print('payload', W_Payload)
    #W_EverythingElse = 2.1296*W_0**-0.159*W_0
    #W_EverythingElse = 1.1318*W_0**-0.098*W_0
    W_EverythingElse = 2.0834*W_0**-0.162*W_0
    #print('everything else', W_EverythingElse)
    # Estimate wing weight
    Wing_S_ref = W_0/DesignW_S
    W_Wing = Wing_S_ref * 10
    #print('wing', W_Wing)
    
    # Estimate weight of electric motor
    W_EM = W_0/EM_W_P/3.16
    # Estimate weight of gas turbine
    W_Propulsion = 1.35*(1583+0.24*W_0/10)
    W_PowerTrain = W_EM + W_Propulsion*N_Engine + 397

    W_E = W_PowerTrain + W_Wing + W_EverythingElse
    #print('W_E', W_E)

    # Range equation
    #
    # Convert units to metric 
    g = 9.81
    W_Payload = W_Payload/2.205*g
    W_E = W_E/2.205*g
    R1 = R/3.281
    # Specify phi value
    phi = 0.2

    # Specify eta values
    eta_1 = 0.4
    eta_2 = 0.96
    eta_3 = 0.8
    
    A = eta_3*e_f/g*L_D*(eta_1 + eta_2*phi/(1-phi))
    B = phi + e_bat/e_f*(1-phi)
    E0_total = ((np.exp(R1/A)*(W_Payload + W_E)) - (W_Payload + W_E))*e_bat/(g*(B-phi*np.exp(R1/A)))
    E0_fuel = (1-phi)*E0_total
    E0_bat = (phi)*E0_total
    W_f = g*E0_fuel/e_f * 0.225
    W_bat = g*E0_bat/(e_bat*0.8) * 0.225

    #print('fuel', W_f)
    #print('battery', W_bat)

    # Convert weights back to lbs
    W_Payload = W_Payload*2.205/g
    W_E = W_E*2.205/g

    W0_new = W_f + W_bat + W_Payload + W_E
    delta = abs(W0_new-W_0) / abs(W0_new)
    #print('W_0', W0_new)
    #print('')
    W_0 = W0_new




print('')

#R1 = A*np.log((W_E+W_Payload+g/e_bat*E_total*(B))/(W_E+W_Payload+g*phi*E_total/e_bat))

print('revised payload weight is:', W_Payload, file = file1)
print('revised takeoff weight is:', W_0, file = file1)
print('revised empty weight is:', W_E, file = file1)
print('revised battery weight is:', W_bat, file = file1)
print('revised fuel weight is:', W_f, file = file1)





file1.close()
file2.close()
# Convert text files to pdfs
#
#
pdf1 = FPDF()
pdf1.add_page()
pdf1.set_font("Arial", size = 11)
f = open("text1.txt", "r")
for x in f:
    pdf1.cell(200, 10, txt = x, ln = 1, align = 'C')
pdf1.output("pdf1.pdf") 
pdf2 = FPDF()
pdf2.add_page()
pdf2.set_font("Arial", size = 11)
f = open("text.txt", "r")
for x in f:
    pdf2.cell(200, 10, txt = x, ln = 1, align = 'C')
pdf2.output("pdf2.pdf") 
f.close() 

# 
# Merge pdfs
#
pdfs = ['testin.pdf', 'pdf1.pdf', 'testin1.pdf', 'pdf2.pdf']
merger = PdfMerger()
for pdf in pdfs:
    merger.append(pdf)
merger.write("test results ATR.pdf")
merger.close()

#
# Clean up data files
#
os.remove('text.txt')
os.remove('text1.txt')
os.remove('testin.pdf')
os.remove('testin1.pdf')
os.remove('pdf2.pdf')
os.remove('pdf1.pdf')
#plt.show()