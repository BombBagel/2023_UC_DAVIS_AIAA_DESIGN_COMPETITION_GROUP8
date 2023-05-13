import os
mainpath = os.getcwd()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from powertrain_matrix import getpowertrain
from matplotlib.backends.backend_pdf import PdfPages
from fpdf import FPDF
from pypdf import PdfMerger
#from neutral_point import getNeutralPoint
from empennage import getempennage
from NP_CG_SM_rev1 import NP_CG_SM
import Constraint_Functions as CF
import weight_components 
import time
from tkinter import *
from sys import exit
from Direct_op_cost import getOperatingCost

# -----------------Functions-------------------
# Function for estimating weight
## 
###
def Calc_weight(A, B, pilots, attendants, passengers, R, L_D, E, pe, v, D_P, N_Prop, N_Engine, Single_Run = False):
    # Calculate weight of crew and payload
    W_crew = (190+30)*(pilots+attendants)
    W_passengers = (200+40)*(passengers)

    # Establish fuel fractions at different flight segments
    W1_W0 = 0.970 # Takeoff fuel fraction
    W2_W1 = 0.985 # Climb fuel fraction
    W3_W2 = np.exp((-R*((0.5*v) / (550*pe*3600))) / (v*L_D)) # Cruise fuel fraction
    W4_W3 = 0.990 # Descent fuel fraction
    W5_W4 = np.exp((-E*((0.6*v) / (550*pe*3600))) / (0.866*L_D)) # Loiter fuel fraction
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
    
    # Calculate landing fuel fraction for landing constraint
    W1_W0 = 0.970 # Takeoff fuel fraction
    W2_W1 = 0.985 # Climb fuel fraction
    W3_W2 = 1 # No cruise segment
    ELanding = 0.25*3600 # 15 minute loiter time
    W4_W3 = 0.990 # Descent fuel fraction
    W5_W4 = np.exp((-ELanding*((0.6*v) / (550*pe*3600))) / (0.866*L_D)) # Loiter fuel fraction
    WL_WTO = (W1_W0*W2_W1*W3_W2*W4_W3*W5_W4)*1.06 # Landing weight fraction
    if Single_Run == True:
        print('Landing Weight Ratio is:', WL_WTO, file = file1)

    # Recalculate fuel fraction since range has changed
    W3_W2 = np.exp((-R*((0.5*v) / (550*pe*3600))) / (v*L_D)) # Cruise fuel fraction
    W5_W4 = np.exp((-E*((0.6*v) / (550*pe*3600))) / (0.866*L_D)) # Loiter fuel fraction
    W6_W5 = 0.995 # Landing fuel fraction

    # Calculate fuel fraction
    Wf_W0 = (1 - W1_W0*W2_W1*W3_W2*W4_W3*W5_W4*W6_W5)*1.06 # Multiplied by 1.06 to account for trapped and reserved fuel
    WCruise_WTakeoff = W1_W0*W2_W1*W3_W2*1.06 # Cruise weight fraction
    if Single_Run == True:
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
    if Single_Run == True:
        print('Max weight of hybrid electric concept is:', W0, 'lbs', file = file1)
    return W0, WL_WTO, WCruise_WTakeoff, W5_W4

# Function for finding cost of aircraft
## Takes in the row number, i and WE which is the empty weight
### Returns the total cost and prints out the individual production costs
def Calc_cost(v, WE, D_P, N_Prop, N_Engine, E_Bat, Single_Run = False):
    W_Airframe = WE*0.45 # Airframe weight (lbs) (We assume is 45% of empty weight)
    V_H = v/1.668 # Maximum Level Airspeed (KTAS)
    Q = 50 # Number of aircraft to be produced over a 5 year period
    Q_M = 5 # Number of aircraft produced in 1 month
    Q_Proto = 2 # Number of prototype aircraft to be produced 
    CPI = 1.31 # Consumer price index for consideration of inflation since 2012
    # Cost of engineering
    R_Eng = 92 # Rate of engineering labor in US$ per hour (92$/h is recommended) 
    C_Eng = 0.083 * (W_Airframe**0.791) * (V_H**1.521) * (Q**0.183) * 0.67 * 2.00 * 1.03 * 1.03 * 1.66 * R_Eng * CPI
    # Cost of tooling
    R_Tool = 61 # Rate of tooling labor in US$ per hour (61$/h is recommended) 
    C_Tool = 2.1036 * (W_Airframe**0.764) * (V_H**0.899) * (Q**0.178) * (Q_M**0.066) * 2.00 * 0.95 * 1.02 * 1.01 * 1.10 * R_Tool * CPI
    # Cost of manufacturing
    R_Mfg = 53 # Rate of manufacturing labor in US$ per hour (53$/h is recommended) 
    C_Mfg = 20.2588 * (W_Airframe**0.74) * (V_H**0.543) * (Q**0.524) * 0.75 * 1.25 * 1.01 * 1.10 * R_Mfg * CPI
    # Cost of development support
    C_Dev = 0.06458 * (W_Airframe**0.873) * (V_H**1.89) * (Q_Proto**0.346) * CPI * 0.5 * 1.5 * 1.01 * 1.03 * 1.05
    # Cost of flight test operations
    C_Ft = 0.009646 * (W_Airframe**1.16) * (V_H**1.3718) * (Q_Proto**1.281) * CPI * 0.5 * 1.5
    # Cost of quality control
    C_Qc = 0.13 * (C_Mfg) * 0.5 * 1.5 * 1.5
    # Cost of materials
    C_Mat = 24.896 * (W_Airframe**0.689) * (V_H**0.624) * (Q*0.792) * CPI * 0.75 * 1.02 * 1.01 * 1.05 
    # Cost of electric motors
    P_Em = 3500 # Rated power P (hp)
    C_Em = 174 * N_Engine * P_Em * CPI
    # Cost of power management system 
    P_EmTotal = P_Em * N_Engine
    C_Pms = 150 * P_EmTotal * CPI
    # Cost of battery (assuming 200$/kWh)
    E_Bat = E_Bat/1000 # Convert to kWh
    C_Bat = 200 * E_Bat * CPI
    # Cost of Propeller
    P_SHP = 7500 # Shaft horsepower
    C_Prop = 210 * N_Prop * CPI * D_P**2 * ((P_SHP/D_P)**0.12)
    if Single_Run == True:
        print('Cost of engineering for concept is', C_Eng/Q, '$', file = costfile)
        print('Cost of tooling for concept is', C_Tool/Q, '$', file = costfile)
        print('Cost of manufacturing for concept is', C_Mfg/Q, '$', file = costfile)
        print('Cost of development for concept is', C_Dev/Q, '$', file = costfile)
        print('Cost of flight test operations for concept is', C_Ft/Q, '$', file = costfile)
        print('Cost of quality control for concept is', C_Qc/Q, '$', file = costfile)
        print('Cost of materials for concept is', C_Mat/Q, '$', file = costfile)
        print('Cost of electric motors for concept is', C_Em/Q, '$', file = costfile)
        print('Cost of power management system for concept is', C_Pms/Q, '$', file = costfile)
        print('Cost of battery for concept is', C_Bat, '$', file = costfile)
        print('Cost of propeller for concept is', C_Prop/Q, '$', file = costfile)

    TotalCost = (C_Eng/Q + C_Tool/Q + C_Mfg/Q + C_Dev/Q + C_Ft/Q + C_Mat/Q + C_Em/Q + C_Pms/Q + C_Bat + C_Prop/Q)*1.15
    return TotalCost, C_Bat

# Function to generate drag polar
# 
# returns parasite drag for different configurations
def dragPolar(deltaTakeoffFlaps,deltaLandingFlaps,deltaLandingGear,e_Clean,e_Takeoff,e_Landing, Single_Run = False):
    if Single_Run == True:
        fig = plt.figure()
    elif Single_Run == False:
        fig = 0
    # Estimate wetted area from historical data
    c = -0.0866 # Regression constant from Table 3.5 Roskam vol 1 for twin propeller plane
    d = 0.8099 # Regression constant from Table 3.5 Roskam vol 1 for twin propeller plane
    S_wet = (10**c)*(W_0)**d # Wetted area
    y = 10**(0.31739*np.log10(W_0) + 0.26161) # Equation estimated from Torenbeek, 1990, Fig. 7.4, takeoff weight divided by wing loading
    S = W_0/y # Value estimated from Torenbeek, 1990, Fig. 7.4, takeoff weight divided by wing loading
    C_f = 0.0026 # Value obtained from Raymer, Table 12.3
    C_Do = C_f*S_wet/S

    # Clean configuration
    C_L_Clean = np.linspace(-1.7, 1.7, 50)
    C_D_Clean = C_Do + C_L_Clean**2/(np.pi*AR*e_Clean) 
    C_D_o_Clean = C_Do # Clean config parasite drag
    if Single_Run == True:
        print('Clean config drag polar is:', round(C_D_o_Clean, 4), '+', round(1/(np.pi*AR*e_Clean), 4), 'C_L^2', file = file2)
        plt.plot(C_D_Clean, C_L_Clean)

    # Takeoff configuration 
    C_L_Takeoff = np.linspace(-2.0, 2, 50) 
    C_D_Takeoff = (C_Do+deltaTakeoffFlaps) + C_L_Takeoff**2/(np.pi*AR*e_Takeoff) # Takeoff flaps, gear up, add 0.02 to parasite drag coefficient, Roskam vol 1, Table 3.6
    C_D_o_ToGUp = C_Do + deltaTakeoffFlaps # Takeoff config gear up
    if Single_Run == True:
        print('Takeoff config gear up drag polar is:', round(C_D_o_ToGUp, 4), '+', round(1/(np.pi*AR*e_Takeoff), 4), 'C_L^2', file = file2)
        plt.plot(C_D_Takeoff, C_L_Takeoff)
    C_D_Takeoff = (C_Do+deltaTakeoffFlaps+deltaLandingGear) + C_L_Takeoff**2/(np.pi*AR*e_Takeoff) # Takeoff flaps, gear down, add 0.02+0.025 to parasite drag coefficient, Roskam vol 1, Table 3.6
    C_D_o_ToGDo = C_Do + deltaTakeoffFlaps + deltaLandingGear # Takeoff config gear down
    if Single_Run == True:
        print('Takeoff config gear down drag polar is:', round(C_D_o_ToGDo, 4), '+', round(1/(np.pi*AR*e_Takeoff), 4), 'C_L^2', file = file2)
        plt.plot(C_D_Takeoff, C_L_Takeoff)

    # Landing configuration
    C_L_Landing = np.linspace(-2.6, 2.6, 50)
    C_D_Landing = (C_Do + deltaLandingFlaps) + C_L_Landing**2/(np.pi*AR*e_Landing) # Landing flaps, gear up, add 0.075 to parasite drag coefficient, Roskam vol 1, Table 3.6
    C_D_o_LaGUp = C_Do + deltaLandingFlaps # Landing config gear up
    if Single_Run == True:
        print('Landing config gear up drag polar is:', round(C_D_o_LaGUp, 4), '+', round(1/(np.pi*AR*e_Landing), 4), 'C_L^2', file = file2)
        plt.plot(C_D_Landing, C_L_Landing)
    C_D_Landing = (C_Do + deltaLandingFlaps + deltaLandingGear) + C_L_Landing**2/(np.pi*AR*e_Landing) # Landing flaps, gear down, add 0.075+0.025 to parasite drag coefficient, Roskam vol 1, Table 3.6
    C_D_o_LaGDo = C_Do + deltaLandingFlaps + deltaLandingGear # Landing config gear down
    if Single_Run == True:
        print('Landing config gear down drag polar is:', round(C_D_o_LaGDo, 4), '+', round(1/(np.pi*AR*e_Landing), 4), 'C_L^2', file = file2)
    
    # Plotting drag polars
    if Single_Run == True:
        plt.plot(C_D_Landing, C_L_Landing)
        plt.legend(['Clean', 'Takeoff flaps, gear up', 'Takeoff flaps, gear down', 'Landing flaps, gear up', 'Landing flaps, gear down'])
        plt.xlabel(r'$C_{D}$')
        plt.ylabel(r'$C_{L}$')
        plt.title('Drag Polar Graphs')
        plt.savefig('Drag Polar.png')
    return fig, C_D_o_ToGDo, C_D_o_ToGUp, C_D_o_Clean, C_D_o_LaGDo

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
    W_P_Takeoff_Matrix = ComponentPowerLoading(W_P_Takeoff_Main, MatrixA)
    W_P_TakeoffClimb_Matrix = ComponentPowerLoading(W_P_TakeoffClimb_Main, MatrixA)
    W_P_TransitionClimb_Matrix = ComponentPowerLoading(W_P_TransitionClimb_Main, MatrixA)
    W_P_SSClimb_Matrix = ComponentPowerLoading(W_P_SSClimb_Main, MatrixA)
    W_P_EnRouteClimb_Matrix = ComponentPowerLoading(W_P_EnRouteClimb_Main, MatrixA)
    W_P_BalkedOEI_Matrix = ComponentPowerLoading(W_P_BalkedOEI_Main, MatrixA)
    W_P_BalkedAEO_Matrix = ComponentPowerLoading(W_P_BalkedAEO_Main, MatrixA)
    W_P_Cruise_Matrix = ComponentPowerLoading(W_P_Cruise_Main, MatrixB)
    W_P_Ceiling_Matrix = ComponentPowerLoading(W_P_Ceiling_Main, MatrixB)

# Plotting function
##
###
def PlotPowerLoadingMatrix(n, Single_Run = False): 
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
    if Single_Run == True:
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
        DesignW_S, DesignW_P = ShadeSubplots(n, True)
        axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].plot(DesignW_S, DesignW_P, 'ro', markersize = '4')
        axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_title('Graph of Power Loading against Wing Loading for ' + Components[i])
        axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_xlabel(r'$W_{0}/S_{w}$' '[lb/ft^2]')
        axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].set_ylabel(r'$W_{0}/P$' '[hp/ft^2]')
    elif Single_Run == False:
        DesignW_S, DesignW_P = ShadeSubplots(n)
    #axs[int(np.ceil(n/3)-1),int(n-3*(np.floor(n/3))-1)]..ylabel(r"$\frac{{W}}{{{}}}$".format(Components[i]) + "[lb/hp]")
    if Components[i] == 'P_e1':
        global EM_W_P
        EM_W_P = DesignW_P

# Shading function for main graph
## Also returns optimal point
###
def Shade(Single_Run = False):
    min = np.minimum(W_P_BalkedOEI_Main, W_P_Cruise_Main)
    min = np.minimum(min, W_P_Takeoff_Main)
    min = np.minimum(min, W_P_Ceiling_Main)
    min = np.minimum(min, W_P_BalkedAEO_Main)
    min = np.minimum(min, W_P_EnRouteClimb_Main)
    min = np.minimum(min, W_P_SSClimb_Main)
    min = np.minimum(min, W_P_TakeoffClimb_Main)
    min = np.minimum(min, W_P_TransitionClimb_Main)
    if Single_Run == True:
        plt.fill_between(W_SArray, min, 0, color='grey')
        plt.fill_betweenx(W_PArray, np.minimum(W_SLanding, W_S_Stall), 100*np.ones(graphPoints), color = 'white')
    DesignW_S, DesignW_P = PlotPoint(min)
    return DesignW_S, DesignW_P

# Shading function for subplots 
##
###
def ShadeSubplots(n, Single_Run = False):
    min = np.minimum(W_P_BalkedOEI, W_P_Cruise)
    min = np.minimum(min, W_P_Takeoff)
    min = np.minimum(min, W_P_Ceiling)
    min = np.minimum(min, W_P_BalkedAEO)
    min = np.minimum(min, W_P_EnRouteClimb)
    min = np.minimum(min, W_P_SSClimb)
    min = np.minimum(min, W_P_TakeoffClimb)
    min = np.minimum(min, W_P_TransitionClimb)
    if Single_Run == True:
        axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].fill_between(W_SArray, min, 0, color='grey')
        axs[int(np.ceil(n/2)-1),int(n-2*(np.ceil(n/2)-1)-1)].fill_betweenx(W_PArray, np.minimum(W_SLanding, W_S_Stall), 100*np.ones(graphPoints), color = 'white')
    DesignW_S, DesignW_P = PlotPoint(min)
    return DesignW_S, DesignW_P

# Function to find best point in shaded area
##
###
def PlotPoint(min):
    DesignW_S = np.floor(np.minimum(W_SLanding, W_S_Stall)[1]) - 5
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























#
# -----------------Main code-------------------
#

#
# ------------------GUI Set Up-------------------
#
root = Tk()
root.title("Set up")
root.geometry('350x200')
# Create options for first dropdown menu to choose analyzed data
RunOptions = ['Single Run',
           'Trade Study']
Clicked = StringVar()
Clicked.set('Single Run')
# Show additional options if two output options is picked
def showSecondVariable(choice):
    if choice == 'Two Outputs':
        drop2.pack()
        drop3.pack()
    elif choice == 'Single Output':
        drop2.pack()
        drop3.forget()
# Create options for second dropdown menu to choose number of anayzed data
options1 = ['Single Output',
            'Two Outputs']
Clicked1 = StringVar()
Clicked1.set('Single Output')
drop1 = OptionMenu( root, Clicked1, *options1, command=showSecondVariable)
# Create options for third dropdown menu to choose first analyzed data
options2 = ['Fuel Weight',
           'Battery Weight',
           'Takeoff Weight']
Clicked2 = StringVar()
Clicked2.set('Fuel Weight')
drop2 = OptionMenu( root , Clicked2 , *options2 )
# Create options for third dropdown menu to choose first analyzed data
Clicked3 = StringVar()
Clicked3.set('Fuel Weight')
drop3 = OptionMenu( root , Clicked3 , *options2)
# Show additional options if trade study option is picked
def showRunCase(choice):
    if choice == 'Trade Study':
        drop1.pack()
        Clicked1.set('Single Output')
        drop2.pack()
    elif choice == 'Single Run':
        drop1.pack_forget()
        drop2.pack_forget()
        drop3.pack_forget()
# Create first dropdown menu to choose run type
drop = OptionMenu( root , Clicked , *RunOptions, command=showRunCase )
drop.pack()
# Close GUI and lock in selection
def proceed():
    global RunSelection 
    global OutputData1
    global OutputData2
    global OutputNumber
    RunSelection = Clicked.get()
    OutputNumber = Clicked1.get()
    if RunSelection == 'Trade Study':
        if OutputNumber == 'Single Output':
            OutputData1 = Clicked2.get()
        elif OutputNumber == 'Two Outputs':
            OutputData1 = Clicked2.get()
            OutputData2 = Clicked3.get()
    root.destroy()
# Create button to confirm selection
button = Button( root , text = "Proceed" , command = proceed ).pack(side='bottom')
root.mainloop()
print(RunSelection)
if 'OutputData2' in locals():
    if OutputData1 == OutputData2:
        print('Choose two different data sets')
        exit(0)
#print(OutputData1)

#
# --------------------Set up--------------------
#
os.chdir(mainpath)
# Read excel file containing aircraft data
df1 = pd.read_excel('Aircraft Data.xlsx')

# Create folder to store results
if not os.path.exists('Results'):
   os.makedirs('Results')
cwd = os.getcwd()
os.chdir(os.path.join(cwd, 'Results'))

# Extract data from dataframe containing aircraft data
j = 0
# Payload for aircraft
A = df1["Raymer \nCoefficient A"][j]
B = df1["Raymer \nCoefficient B"][j]
pilots = df1["Pilots"][j]
passengers = df1["Passengers"][j]
# Attendants will be expressed as a function of crew since we need 1 for every 50 passengers
attendants = np.ceil(passengers/50)
Bat_SpeEnergy = df1["Battery Specific Energy [wh/kg]"][j] / 2.205
# Establish Flight Profile
R = df1["Range \n(nmi)"][j]*5280 # Range in ft
L_D = df1["L/D"][j] # Max lift to drag ratio
E = df1["Endurance \n(hours)"][j]*3600 # Endurance in seconds
pe = df1["Propellor \nEfficiency"][j] # Propeller efficiency
eta_p = df1["Propellor \nEfficiency"][j] # Propeller efficiency
v = df1["Speed (knots)"][j]*1.668 # Speed in ft/s 
# Extract airplane parameters
N_Prop = df1["No. Propellers"][j] # Number of propellers
D_P = df1["Propeller \nDiameter (ft)"][j] # Diameter of propeller in ft
N_Engine = df1["No. Engines"][j] # Number of required engines
AR = df1["AR"][j] # Aspect ratio
C_L_maxL = df1["Clmax \nLanding"][j] # Max Cl of lift 
C_L_maxCr = df1["Clmax Cruise"][j] # Max Cl for cruise
C_L_maxTO = df1["Clmax \nTakeoff"][j] # Max Cl for takeoff
deltaTakeoffFlaps = df1['delta\nTakeoff\nFlaps'][j]
deltaLandingFlaps = df1['delta\nLanding\nFlaps'][j]
deltaLandingGear = df1['delta\nLanding\nGear'][j]
e_Clean = df1["e Clean"][j]
e_Takeoff = df1["e Takeoff"][j]
e_Landing = df1["e Landing"][j]
c_root = df1["Wing Root Chord [ft]"][j]
c_tip = df1["Wing Tip Chord [ft]"][j]
b_w = df1["Wing Span [ft]"][j]
S_w = df1["Wing Area [ft^2]"][j]
L_ht = df1["Wing Aerodynamic Center to Horizontal Tail Aerodynamic Center [ft]"][j]
L_vt = df1["Wing Aerodynamic Center to Vertical Tail Aerodynamic Center [ft]"][j]
Battery_Percentage = df1["Battery Percentage"][j]

# If case for a single run
if RunSelection == 'Single Run':
    # Create txt file to write data 
    file1 = open('text.txt', 'w') # file to write weight data
    costfile = open('costfile.txt', 'w') # file to write cost data
    file2 = open('text1.txt', 'w') # file to write preliminary design estimate data
    file4 = open('text2.txt', 'w') # file to write empennage suggested sizing data

    # Create pdf files to write data
    pp = PdfPages('Graph.pdf')
    pp1 = PdfPages('Graph1.pdf')

    #
    # -----------------Calculations-------------------
    #

    # Section to carry our preliminary weight calculations
    # 
    W_0, WL_WTO, WCruise_WTakeoff, W5_W4 = Calc_weight(A, B, pilots, attendants, passengers, R, L_D, E, pe, v, D_P, N_Prop, N_Engine, True)

    # Section to carry our drag polar calculations
    #
    fig1,C_D_o_ToGDo, C_D_o_ToGUp, C_D_o_Clean, C_D_o_LaGDo = dragPolar(
        deltaTakeoffFlaps,deltaLandingFlaps,deltaLandingGear,e_Clean,e_Takeoff,e_Landing, True)
    pp.savefig(fig1, bbox_inches="tight")

    # Plotting W/P to W/S graph section
    #
    # initialize number of points on graph
    graphPoints = 800
    rho = 0.002048 # Density at 5000 ft
    rho_SL = 0.00237 # Density at sea level
    rho_28000 = (0.0009567 + 0.002048)/2 # Density at FL280 !!!!!!!MIDPOINT NOW
    W_SArray = np.linspace(0.001, 200, graphPoints)
    W_PArray = np.linspace(0.001, 300, graphPoints)
    # Solve stall constraint
    W_S_Stall = CF.stallConstraint(rho,C_L_maxL,graphPoints)
    # Takeoff constraint
    W_P_Takeoff_Main = CF.takeoffConstraint(C_D_o_ToGDo,rho_SL,N_Engine,D_P,W_SArray,C_L_maxTO)
    # Landing constraint
    W_SLanding = CF.landingConstraint(C_L_maxL,rho,rho_SL,WL_WTO,graphPoints)
    # Solve climb constraints
    W_P_TakeoffClimb_Main,W_P_TransitionClimb_Main,W_P_SSClimb_Main,W_P_EnRouteClimb_Main,W_P_BalkedAEO_Main,W_P_BalkedOEI_Main = CF.climbConstraints(
        C_D_o_ToGUp,C_L_maxTO,W_SArray,rho_28000,e_Takeoff,AR,N_Engine,eta_p,graphPoints,C_D_o_ToGDo,C_D_o_Clean,C_L_maxCr,
                        e_Clean,C_D_o_LaGDo,C_L_maxL,e_Landing,rho,WL_WTO,Single_Run = True,file = file2)
    # Cruise constraint
    W_P_Cruise_Main = CF.cruiseConstraint(v,C_D_o_Clean,rho_28000,W_SArray,WCruise_WTakeoff,AR,e_Clean,eta_p)
    # Ceiling constraint
    W_P_Ceiling_Main = CF.ceilingConstraint(C_D_o_Clean,W_SArray,rho_28000,e_Clean,AR,eta_p,graphPoints)
    # initialize figure
    fig2 = plt.figure(figsize=(12, 10.5))
    # Plot values
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
    # Shade feasible region on graph
    DesignW_S, DesignW_P = Shade(True)
    # Plot out design point
    plt.plot(DesignW_S, DesignW_P, 'ro', markersize = '4')
    # Label details on graph
    plt.legend(fontsize='12', loc='center left')
    plt.xlabel(r'$\frac{W_{0}}{S_{w}}$' '[lb/ft^2]', fontsize='14')
    plt.ylabel(r'$\frac{W_{0}}{P}$' ' [lb/hp]', fontsize='14')
    plt.ylim([0,25])
    plt.xlim([0, 100])
    plt.title('Graph of Power Loading against Wing Loading')
    plt.savefig('Power Loading against Wing Loading.png',bbox_inches="tight")
    pp1.savefig(fig2, bbox_inches="tight")

    # Section to carry out component power loading and wing loading graphs
    #
    # Let MatrixA = matrix for battery discharging and MatrixB = matrix for battery charging 
    # get matrix for serial powertrain
    architecture,state='serial','propulsion'
    MatrixA,Components = getpowertrain(architecture,state)
    architecture,state='serial','charge-both-propulsion'
    MatrixB,Components = getpowertrain(architecture,state)
    # Solve for constraints
    ConstraintSolve()
    fig3 = InitializeSubplots(Components)
    # Plot constraints for each component
    for i in range(len(Components)):
        PlotPowerLoadingMatrix(i, True)
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines[0:11], labels[0:11], loc = 'lower center', ncol = 5, bbox_to_anchor=(0.5,0.04), bbox_transform=fig.transFigure, fontsize="large")
    plt.savefig('Component Power Loadings1.png',bbox_inches="tight")
    pp1.savefig(fig3, bbox_inches="tight")
    # get matrix for parallel powertrain
    architecture,state='parallel','propulsion'
    MatrixA,Components = getpowertrain(architecture,state)
    architecture,state='parallel','charge-both-propulsion'
    MatrixB,Components = getpowertrain(architecture,state)
    # Solve for constraints
    ConstraintSolve()
    fig4 = InitializeSubplots(Components)
    # Plot constraints for each component
    for i in range(len(Components)):
        PlotPowerLoadingMatrix(i, True)
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines[0:11], labels[0:11], loc = 'lower center', ncol = 5, bbox_to_anchor=(0.5,0.04), bbox_transform=fig.transFigure, fontsize="large")
    plt.savefig('Component Power Loadings2.png',bbox_inches="tight")
    pp1.savefig(fig4, bbox_inches="tight")

    # Close pdf files since all data has been written
    pp.close()
    pp1.close()

    # Section for preliminary plane weight estimate
    #
    e_f = 12000*3600 # Specific energy of fuel in J/kg
    e_bat = Bat_SpeEnergy*3600*2.205 # Specific energy of battery in J/kg   
    W_crew = (190+30)*(pilots+attendants) # Weight of crew
    W_passengers = (200+40)*(passengers) # Weight of passengers
    W_Payload = W_crew + W_passengers # Weight of payload
    # Set tolerance
    sigma = 10**-6
    delta = 2*sigma
    while delta > sigma:
        W1_W0 = 0.99811 # Takeoff fuel fraction
        W2_W1 = 0.992 # Climb fuel fraction
        W4_W3 = 0.990 # Descent fuel fraction
        W6_W5 = 0.995 # Landing fuel fraction

        # Find weight of fuel needed for takeoff
        W_fuelTakeoff = W_0 - (W1_W0 * W_0)
        #print(W_fuelTakeoff)
        W_batTakeoff = W_fuelTakeoff*Battery_Percentage*5443.1084/Bat_SpeEnergy
        W_fuelTakeoff = W_fuelTakeoff * (1-Battery_Percentage) # Reduce to account for battery providing power
        # Find weight of fuel needed for climb
        W_fuelClimb = (W1_W0 * W_0) - (W2_W1 * W1_W0 * W_0)
        W_batClimb = W_fuelClimb*Battery_Percentage*5443.1084/Bat_SpeEnergy
        W_fuelClimb = W_fuelClimb * (1-Battery_Percentage) # Reduce to account for battery providing power

        # Find everything else weight (From historical data)
        W_EverythingElse = 2.0834*W_0**-0.162*W_0

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

        # Range equation
        #
        # Convert units to metric 
        g = 9.81
        W_Payload = W_Payload/2.205*g
        W_E = W_E/2.205*g
        R1 = R/3.281*0.55 # Multiply by 0.55 to account for cruise distance 
        # Specify phi value
        phi = 0 # this is 0 since we are not using any batery in cruise
        # Specify eta values
        eta_1 = 0.4
        eta_2 = 0.96
        eta_3 = 0.8
        A1 = eta_3*e_f/g*L_D*(eta_1 + eta_2*phi/(1-phi))
        B1 = phi + e_bat/e_f*(1-phi)
        E0_total = ((np.exp(R1/A1)*(W_Payload + W_E)) - (W_Payload + W_E))*e_bat/(g*(B1-phi*np.exp(R1/A1)))
        E0_fuel = (1-phi)*E0_total
        E0_bat = (phi)*E0_total
        W_fuelCruise = g*E0_fuel/e_f * 0.225
        #print(W_fuelCruise)
        W_batCruise = g*E0_bat/(e_bat*0.8) * 0.225 
        
        # Convert weights back to lbs
        W_Payload = W_Payload*2.205/g
        W_E = W_E*2.205/g

        # Find cruise weight fraction
        W3_W2 = ((W2_W1 * W1_W0 * W_0) - W_fuelCruise)/(W2_W1 * W1_W0 * W_0)
        
        # Find weight of fuel needed for descent
        W_fuelDescent = (W3_W2 * W2_W1 * W1_W0 * W_0) - (W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0)
        #print(W_fuelDescent)
        W_batDescent = W_fuelDescent*Battery_Percentage*5443.1084/Bat_SpeEnergy*0.95 # 0.95 to account for battery being charged
        W_fuelDescent = W_fuelDescent * (1-Battery_Percentage) # Reduce to account for battery providing power, 0.95 for battery cruise battery charging
        # Find weight of fuel needed for Loiter
        W_fuelLoiter = (W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0) - (W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0)
        W_batLoiter = W_fuelLoiter*Battery_Percentage*5443.1084/Bat_SpeEnergy*0.95 # 0.95 to account for battery being charged
        W_fuelLoiter = W_fuelLoiter * (1-Battery_Percentage) # Reduce to account for battery providing power, 0.95 for battery cruise battery charging
        # Find weight of fuel needed for landing
        W_fuelLanding = (W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0) - (W6_W5 * W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0)
        W_batLanding = W_fuelLanding*Battery_Percentage*5443.1084/Bat_SpeEnergy*0.95 # 0.95 to account for battery being charged
        W_fuelLanding = W_fuelLanding * (1-Battery_Percentage) # Reduce to account for battery providing power, 0.95 for battery cruise battery charging

        W_bat = W_batTakeoff + W_batClimb + W_batDescent + W_batLoiter + W_batLanding
        W_fuel = W_fuelTakeoff + W_fuelClimb + W_fuelCruise + W_fuelDescent + W_fuelLoiter + W_fuelLanding 
        
        W0_new = W_fuel + W_bat + W_Payload + W_E 

        delta = abs(W0_new-W_0) / abs(W0_new)
        W_0 = W0_new
    print('')
    print('revised payload weight is:', W_Payload, file = file1)
    print('revised takeoff weight is:', W_0, file = file1)
    print('revised empty weight is:', W_E + W_bat, file = file1)
    print('revised battery weight is:', W_bat, file = file1)
    print('revised powertrain weight is:', W_PowerTrain, file = file1)
    print('revised total fuel weight is:', W_fuel, file = file1)
    print('revised takeoff fuel weight is:', W_fuelTakeoff, file = file1)
    print('revised climb fuel weight is:', W_fuelClimb, file = file1)
    print('revised cruise fuel weight is:', W_fuelCruise, file = file1)
    print('revised descent fuel weight is:', W_fuelDescent, file = file1)
    print('revised loiter fuel weight is:', W_fuelLoiter, file = file1)
    print('revised landing fuel weight is:', W_fuelLanding, file = file1)
    print('wing loading is:', DesignW_S, file = file1)
    print('power loading is:', DesignW_P, file = file1)
    print('wing area is:', Wing_S_ref, file = file1)

    # Section to recalculate cost weight obtained from prelminary weight estimate
    #
    Bat_Energy = W_bat * Bat_SpeEnergy
    print('Batter energy is', Bat_Energy, file = file1)
    Unit_Cost, C_Bat = Calc_cost(v, W_E, D_P, N_Prop, N_Engine, Bat_Energy, True)
    print('Total cost of hybrid electric concept is:', Unit_Cost, '$', file = costfile)

    # Section to create Power vs Area graph 
    # 
    plt.figure(figsize=(12, 10.5))
    plt.plot(W_0/W_S_Stall, W_0/W_PArray, label = 'Stall', color = 'C8')
    plt.plot(W_0/W_SArray, W_0/W_P_Takeoff_Main, label = 'Takeoff', color = 'C1')
    plt.plot(W_0/W_SLanding, W_0/W_PArray, label = 'Landing', color = 'C0')
    plt.plot(W_0/W_SArray, W_0/W_P_TakeoffClimb_Main, label = 'Takeoff climb', color = 'red')
    plt.plot(W_0/W_SArray, W_0/W_P_TransitionClimb_Main, label = 'Transition climb', color = 'brown', linestyle = 'dashed')
    plt.plot(W_0/W_SArray, W_0/W_P_SSClimb_Main, label = 'Second segment climb', color = 'purple')
    plt.plot(W_0/W_SArray, W_0/W_P_EnRouteClimb_Main, label = 'En-route climb', color = 'darkgreen', linestyle = 'dotted')
    plt.plot(W_0/W_SArray, W_0/W_P_BalkedAEO_Main, label = 'Balked landing climb (AEO)', color = 'black')
    plt.plot(W_0/W_SArray, W_0/W_P_BalkedOEI_Main, label = 'Balked landing climb (OEI)', color = 'C2')
    plt.plot(W_0/W_SArray, W_0/W_P_Cruise_Main, label = 'Cruise', color = 'cyan')
    plt.plot(W_0/W_SArray, W_0/W_P_Ceiling_Main, label = 'Ceiling', color = 'grey')
    # Plot out design point on power vs area graph
    plt.plot(W_0/DesignW_S, W_0/DesignW_P, 'ro', markersize = '4')
    # Shade feasible area for power vs area graph
    max = np.maximum(W_0/W_P_BalkedOEI_Main, W_0/W_P_Cruise_Main)
    max = np.maximum(max, W_0/W_P_Takeoff_Main)
    max = np.maximum(max, W_0/W_P_Ceiling_Main)
    max = np.maximum(max, W_0/W_P_BalkedAEO_Main)
    max = np.maximum(max, W_0/W_P_EnRouteClimb_Main)
    max = np.maximum(max, W_0/W_P_SSClimb_Main)
    max = np.maximum(max, W_0/W_P_TakeoffClimb_Main)
    max = np.maximum(max, W_0/W_P_TransitionClimb_Main)
    plt.fill_between(W_0/W_SArray, max, np.ones(len(max))*100000, color='grey')
    plt.fill_betweenx(W_0/W_PArray, 100*np.ones(graphPoints), np.maximum(W_0/W_SLanding, W_0/W_S_Stall), color = 'white')
    # Fill in details for power vs area graph
    plt.legend(fontsize='12', loc='upper right')
    plt.xlabel(r'${S}$' '[ft^2]')
    plt.ylabel(r'${P}$' ' [hp]')
    plt.ylim([0, W_0/5])
    plt.xlim([W_0/(DesignW_S+30), W_0/40])
    plt.title('Graph of Power vs Wing Area')
    plt.savefig('Power against Wing Area.png',bbox_inches="tight")

    # Section to calculate required empennage sizes
    # 
    print("Empennage Sizing", file = file4)
    S_ht,S_vt, cbar_w = getempennage(c_root,c_tip,b_w,S_w,L_ht,L_vt)
    print("Recommended horizontal tail area is: ", S_ht, file = file4)
    print("Recommended vertical tail area is: ", S_vt, file = file4)
    print("", file = file4)
    print("Mean wing chord is: ", cbar_w, file = file4)


    # Section to calculate weights for components
    #
    (W_w, W_ht, W_vt, W_fus, W_mlg, W_nlg, W_ng, W_eCon, W_start, W_fs, W_fc, W_APUin, W_inst, W_hyd,
                W_ele, W_avio, W_fur, W_ac, W_dIce, W_prop) = weight_components.calcCompWeight(W_0, WL_WTO*W_0, W_fuel, W_Propulsion, W_bat)

    # Section to calculate cg, np and static margin
    #
    NP_CG_SM(v, W_w, W_ht, W_vt, W_fus, W_mlg, W_nlg, W_ng, W_eCon, W_start, W_fs, W_fc, W_APUin, 
            W_inst, W_hyd,W_ele, W_avio, W_fur, W_ac, W_dIce, W_prop, W_fuel, W_bat, W_Propulsion)
    
    getOperatingCost(W_fuel,Bat_Energy,W_0,Unit_Cost,C_Bat)

    #
    # -------------------File clean up--------------------------
    #
    # Close text files
    costfile.close()
    file1.close()
    file2.close()
    file4.close()
    # Convert text files to pdfs
    #
    # converting text1 to pdf1
    pdf1 = FPDF()
    pdf1.add_page()
    pdf1.set_font("Arial", size = 11)
    f = open("text1.txt", "r")
    for x in f:
        pdf1.cell(200, 10, txt = x, ln = 1, align = 'C')
    pdf1.output("pdf1.pdf") 
    # converting text to pdf2
    pdf2 = FPDF()
    pdf2.add_page()
    pdf2.set_font("Arial", size = 11)
    f = open("text.txt", "r")
    for x in f:
        pdf2.cell(200, 10, txt = x, ln = 1, align = 'C')
    pdf2.output("pdf2.pdf") 
    f.close() 
    # converting text2 to pdf3
    pdf3 = FPDF()
    pdf3.add_page()
    pdf3.set_font("Arial", size = 11)
    f = open("text2.txt", "r")
    for x in f:
        pdf3.cell(200, 10, txt = x, ln = 1, align = 'C')
    pdf3.output("pdf3.pdf") 
    f.close() 
    # converting costfile to costpdf
    costpdf = FPDF()
    costpdf.add_page()
    costpdf.set_font("Arial", size = 11)
    f = open("costfile.txt", "r")
    for x in f:
        costpdf.cell(200, 7.5, txt = x, ln = 1, align = 'C')
    costpdf.output("costpdf.pdf") 
    f.close() 
    # Merge pdfs
    pdfs = ['Graph.pdf','pdf1.pdf','Graph1.pdf','pdf2.pdf','costpdf.pdf','pdf3.pdf','pdf4.pdf','NP_CG_SM.pdf','Graph2.pdf', 'Direct_op_cost.pdf']
    merger = PdfMerger()
    for pdf in pdfs:
        merger.append(pdf)
    merger.write("test results.pdf")
    merger.close()
    # Clean up data files
    for i in pdfs:
        os.remove(i)
    textfiles = ['text.txt','text1.txt','text2.txt','text3.txt','costfile.txt', 'NP_CG_SM.txt']
    for i in textfiles:
        os.remove(i)



# Section for trade study case
#
#
elif RunSelection == 'Trade Study':
    # Identify column for variable that will be changed
    VariableCounter = 0 # Variable to keep count of number of changed variables
    for column in df1:
        if pd.isna(df1[column][1]) == False:
            #print(df1[column])
            if VariableCounter == 0:
                Variables = np.linspace((df1[column][1]), (df1[column][2]), int(df1[column][3]))
                VariableName = column
                VariableCounter += 1
            elif VariableCounter == 1:
                Variables2 = np.linspace((df1[column][1]), (df1[column][2]), int(df1[column][3]))
                VariableName2 = column
                VariableCounter += 1
    if VariableCounter == 1:
        rowNum = len(Variables)
        colNum = 1
        Results = np.ones(len(Variables))
    elif VariableCounter == 2:
        xx, yy = np.meshgrid(Variables, Variables2)
        rowNum = xx.shape[0]
        colNum = yy.shape[1]
        Results = np.ones([rowNum, colNum])
        Results2 = np.ones([rowNum, colNum])
        print(xx)
        print(yy)
    print(VariableName)
    print(Variables)
    if 'VariableName2' in locals():
        print(VariableName2)
        print(Variables2)
    if OutputNumber == 'Two Outputs':
        if VariableCounter == 1:
            print('Not enough variables for a contour plot')
            exit(0)
    for k in range(rowNum):
        #print(Variables[k])
        # Identify what variable is being changed
        if VariableName == 'Battery Percentage':
            Battery_Percentage = Variables[k]
        elif VariableName == 'Speed (knots)':
            v = Variables[k]*1.668 # Speed in ft/s 
        elif VariableName == 'Battery Specific Energy [wh/kg]':
            Bat_SpeEnergy = Variables[k] / 2.205
        elif VariableName == 'AR':
            AR = Variables[k] 
        elif VariableName == 'L/D':
            L_D = Variables[k]
        elif VariableName == 'Clmax Cruise':
            C_L_maxCr = Variables[k]
        for l in range(colNum):
            #print(Variables2[l])
            if VariableCounter == 2:
                if VariableName2 == 'Battery Percentage':
                    Battery_Percentage = Variables2[l]
                elif VariableName2 == 'Speed (knots)':
                    v = Variables2[l]*1.668 # Speed in ft/s 
                elif VariableName2 == 'Battery Specific Energy [wh/kg]':
                    Bat_SpeEnergy = Variables2[l] / 2.205
                elif VariableName2 == 'AR':
                    AR = Variables2[l]
                elif VariableName2 == 'Clmax Cruise':
                    C_L_maxCr = Variables[l]

            #
            # -----------------Calculations-------------------
            #

            # Section to carry our preliminary weight calculations
            # 
            W_0, WL_WTO, WCruise_WTakeoff, W5_W4 = Calc_weight(A, B, pilots, attendants, passengers, R, L_D, E, pe, v, D_P, N_Prop, N_Engine)

            # Section to carry our drag polar calculations
            #
            fig1,C_D_o_ToGDo, C_D_o_ToGUp, C_D_o_Clean, C_D_o_LaGDo = dragPolar(
                deltaTakeoffFlaps,deltaLandingFlaps,deltaLandingGear,e_Clean,e_Takeoff,e_Landing)

            # Plotting W/P to W/S graph section
            #
            # initialize number of points on graph
            graphPoints = 800
            rho = 0.002048 # Density at 5000 ft
            rho_SL = 0.00237 # Density at sea level
            rho_28000 = (0.0009567 + 0.002048)/2 # Density at FL280 !!!!!!!MIDPOINT NOW
            W_SArray = np.linspace(0.001, 200, graphPoints)
            W_PArray = np.linspace(0.001, 300, graphPoints)
            # Solve stall constraint
            W_S_Stall = CF.stallConstraint(rho,C_L_maxL,graphPoints)
            # Takeoff constraint
            W_P_Takeoff_Main = CF.takeoffConstraint(C_D_o_ToGDo,rho_SL,N_Engine,D_P,W_SArray,C_L_maxTO)
            # Landing constraint
            W_SLanding = CF.landingConstraint(C_L_maxL,rho,rho_SL,WL_WTO,graphPoints)
            # Solve climb constraints
            W_P_TakeoffClimb_Main,W_P_TransitionClimb_Main,W_P_SSClimb_Main,W_P_EnRouteClimb_Main,W_P_BalkedAEO_Main,W_P_BalkedOEI_Main = CF.climbConstraints(
                C_D_o_ToGUp,C_L_maxTO,W_SArray,rho_28000,e_Takeoff,AR,N_Engine,eta_p,graphPoints,C_D_o_ToGDo,C_D_o_Clean,C_L_maxCr,
                                e_Clean,C_D_o_LaGDo,C_L_maxL,e_Landing,rho,WL_WTO)
            # Cruise constraint
            W_P_Cruise_Main = CF.cruiseConstraint(v,C_D_o_Clean,rho_28000,W_SArray,WCruise_WTakeoff,AR,e_Clean,eta_p)
            # Ceiling constraint
            W_P_Ceiling_Main = CF.ceilingConstraint(C_D_o_Clean,W_SArray,rho_28000,e_Clean,AR,eta_p,graphPoints)
            # initialize figure
            #fig2 = plt.figure(figsize=(12, 10.5))
            # Plot values
            #plt.plot(W_S_Stall, W_PArray, label = 'Stall', color = 'C8')
            #plt.plot(W_SArray, W_P_Takeoff_Main, label = 'Takeoff', color = 'C1')
            # plt.plot(W_SLanding, W_PArray, label = 'Landing', color = 'C0')
            # plt.plot(W_SArray, W_P_TakeoffClimb_Main, label = 'Takeoff climb', color = 'red')
            # plt.plot(W_SArray, W_P_TransitionClimb_Main, label = 'Transition climb', color = 'brown', linestyle = 'dashed')
            # plt.plot(W_SArray, W_P_SSClimb_Main, label = 'Second segment climb', color = 'purple')
            # plt.plot(W_SArray, W_P_EnRouteClimb_Main, label = 'En-route climb', color = 'darkgreen', linestyle = 'dotted')
            # plt.plot(W_SArray, W_P_BalkedAEO_Main, label = 'Balked landing climb (AEO)', color = 'black')
            # plt.plot(W_SArray, W_P_BalkedOEI_Main, label = 'Balked landing climb (OEI)', color = 'C2')
            # plt.plot(W_SArray, W_P_Cruise_Main, label = 'Cruise', color = 'cyan')
            # plt.plot(W_SArray, W_P_Ceiling_Main, label = 'Ceiling', color = 'grey')
            # Shade feasible region on graph
            DesignW_S, DesignW_P = Shade()
            # Plot out design point
            #plt.plot(DesignW_S, DesignW_P, 'ro', markersize = '4')
            # Label details on graph
            # plt.legend(fontsize='12', loc='center left')
            # plt.xlabel(r'$\frac{W_{0}}{S_{w}}$' '[lb/ft^2]', fontsize='14')
            # plt.ylabel(r'$\frac{W_{0}}{P}$' ' [lb/hp]', fontsize='14')
            # plt.ylim([0,25])
            # plt.xlim([0, 100])
            # plt.title('Graph of Power Loading against Wing Loading')

            # Section to carry out component power loading and wing loading graphs
            #
            # Let MatrixA = matrix for battery discharging and MatrixB = matrix for battery charging 
            # get matrix for serial powertrain
            architecture,state='serial','propulsion'
            MatrixA,Components = getpowertrain(architecture,state)
            architecture,state='serial','charge-both-propulsion'
            MatrixB,Components = getpowertrain(architecture,state)
            # Solve for constraints
            ConstraintSolve()
            #fig3 = InitializeSubplots(Components)
            # Plot constraints for each component
            for i in range(len(Components)):
                PlotPowerLoadingMatrix(i)
            #lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            #lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
            #fig.legend(lines[0:11], labels[0:11], loc = 'lower center', ncol = 5, bbox_to_anchor=(0.5,0.04), bbox_transform=fig.transFigure, fontsize="large")
            # get matrix for parallel powertrain
            architecture,state='parallel','propulsion'
            MatrixA,Components = getpowertrain(architecture,state)
            architecture,state='parallel','charge-both-propulsion'
            MatrixB,Components = getpowertrain(architecture,state)
            # Solve for constraints
            ConstraintSolve()
            #fig4 = InitializeSubplots(Components)
            # Plot constraints for each component
            for i in range(len(Components)):
                PlotPowerLoadingMatrix(i)
            #lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            #lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
            #fig.legend(lines[0:11], labels[0:11], loc = 'lower center', ncol = 5, bbox_to_anchor=(0.5,0.04), bbox_transform=fig.transFigure, fontsize="large")

            # Section for preliminary plane weight estimate
            #
            e_f = 12000*3600 # Specific energy of fuel in J/kg
            e_bat = Bat_SpeEnergy*3600*2.205 # Specific energy of battery in J/kg   
            W_crew = (190+30)*(pilots+attendants) # Weight of crew
            W_passengers = (200+40)*(passengers) # Weight of passengers
            W_Payload = W_crew + W_passengers # Weight of payload
            # Set tolerance
            sigma = 10**-6
            delta = 2*sigma
            while delta > sigma:
                W1_W0 = 0.99811 # Takeoff fuel fraction
                W2_W1 = 0.992 # Climb fuel fraction
                W4_W3 = 0.990 # Descent fuel fraction
                W6_W5 = 0.995 # Landing fuel fraction

                # Find weight of fuel needed for takeoff
                W_fuelTakeoff = W_0 - (W1_W0 * W_0)
                #print(W_fuelTakeoff)
                W_batTakeoff = W_fuelTakeoff*Battery_Percentage*5443.1084/Bat_SpeEnergy
                W_fuelTakeoff = W_fuelTakeoff * (1-Battery_Percentage) # Reduce to account for battery providing power
                # Find weight of fuel needed for climb
                W_fuelClimb = (W1_W0 * W_0) - (W2_W1 * W1_W0 * W_0)
                W_batClimb = W_fuelClimb*Battery_Percentage*5443.1084/Bat_SpeEnergy
                W_fuelClimb = W_fuelClimb * (1-Battery_Percentage) # Reduce to account for battery providing power

                # Find everything else weight (From historical data)
                W_EverythingElse = 2.0834*W_0**-0.162*W_0

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

                # Range equation
                #
                # Convert units to metric 
                g = 9.81
                W_Payload = W_Payload/2.205*g
                W_E = W_E/2.205*g
                R1 = R/3.281*0.55 # Multiply by 0.55 to account for cruise distance 
                # Specify phi value
                phi = 0 # this is 0 since we are not using any batery in cruise
                # Specify eta values
                eta_1 = 0.4
                eta_2 = 0.96
                eta_3 = 0.8
                A1 = eta_3*e_f/g*L_D*(eta_1 + eta_2*phi/(1-phi))
                B1 = phi + e_bat/e_f*(1-phi)
                E0_total = ((np.exp(R1/A1)*(W_Payload + W_E)) - (W_Payload + W_E))*e_bat/(g*(B1-phi*np.exp(R1/A1)))
                E0_fuel = (1-phi)*E0_total
                E0_bat = (phi)*E0_total
                W_fuelCruise = g*E0_fuel/e_f * 0.225
                #print(W_fuelCruise)
                W_batCruise = g*E0_bat/(e_bat*0.8) * 0.225 
                
                # Convert weights back to lbs
                W_Payload = W_Payload*2.205/g
                W_E = W_E*2.205/g

                # Find cruise weight fraction
                W3_W2 = ((W2_W1 * W1_W0 * W_0) - W_fuelCruise)/(W2_W1 * W1_W0 * W_0)
                
                # Find weight of fuel needed for descent
                W_fuelDescent = (W3_W2 * W2_W1 * W1_W0 * W_0) - (W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0)
                #print(W_fuelDescent)
                W_batDescent = W_fuelDescent*Battery_Percentage*5443.1084/Bat_SpeEnergy*0.95 # 0.95 to account for battery being charged
                W_fuelDescent = W_fuelDescent * (1-Battery_Percentage) # Reduce to account for battery providing power, 0.95 for battery cruise battery charging
                # Find weight of fuel needed for Loiter
                W_fuelLoiter = (W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0) - (W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0)
                W_batLoiter = W_fuelLoiter*Battery_Percentage*5443.1084/Bat_SpeEnergy*0.95 # 0.95 to account for battery being charged
                W_fuelLoiter = W_fuelLoiter * (1-Battery_Percentage) # Reduce to account for battery providing power, 0.95 for battery cruise battery charging
                # Find weight of fuel needed for landing
                W_fuelLanding = (W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0) - (W6_W5 * W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0 * W_0)
                W_batLanding = W_fuelLanding*Battery_Percentage*5443.1084/Bat_SpeEnergy*0.95 # 0.95 to account for battery being charged
                W_fuelLanding = W_fuelLanding * (1-Battery_Percentage) # Reduce to account for battery providing power, 0.95 for battery cruise battery charging

                W_bat = W_batTakeoff + W_batClimb + W_batDescent + W_batLoiter + W_batLanding
                W_fuel = W_fuelTakeoff + W_fuelClimb + W_fuelCruise + W_fuelDescent + W_fuelLoiter + W_fuelLanding 
                
                W0_new = W_fuel + W_bat + W_Payload + W_E 

                delta = abs(W0_new-W_0) / abs(W0_new)
                W_0 = W0_new
            if VariableCounter == 1:
                if OutputData1 == 'Fuel Weight':
                    Results[k] = W_fuel
                elif OutputData1 == 'Battery Weight':
                    Results[k] = W_bat
                elif OutputData1 == 'Takeoff Weight':
                    Results[k] = W_0
            if VariableCounter == 2:
                if OutputData1 == 'Fuel Weight':
                    Results[l,k] = W_fuel
                elif OutputData1 == 'Battery Weight':
                    Results[l,k] = W_bat
                elif OutputData1 == 'Takeoff Weight':
                    Results[l,k] = W_0
                if OutputNumber == 'Two Outputs':    
                    if OutputData2 == 'Fuel Weight':
                        Results2[l,k] = W_fuel
                    elif OutputData2 == 'Battery Weight':
                        Results2[l,k] = W_bat
                    elif OutputData2 == 'Takeoff Weight':
                        Results2[l,k] = W_0
    # Remove units from variable name if it exists
    for line in (VariableName.split()):
        if '[' in line:
            VariableNamenoUnits = VariableName.replace(' ' + line, '')
            break
        else:
            VariableNamenoUnits = VariableName
    # Plot for one independent variable
    if VariableCounter == 1:
        plt.plot(Variables, Results)
        plt.xlabel(VariableName)
        # Determine what unit to assign to y axis label
        for line in (OutputData1.split()):
            if line == 'Weight':
                OutputData1withUnits = OutputData1 + ' [lbf]'
                break
        plt.ylabel(OutputData1withUnits)
        plt.title('Graph of ' + OutputData1 + ' against ' + VariableNamenoUnits)
    # Plot for two independent variables
    if VariableCounter == 2:
        # Remove units from variable name if it exists
        for line in (VariableName2.split()):
            if '[' in line:
                VariableName2noUnits = VariableName2.replace(' ' + line, '')
                break
            else:
                VariableName2noUnits = VariableName2 
        # Determine what unit to assign to y axis label
        for line in (OutputData1.split()):
            if line == 'Weight':
                OutputData1withUnits = OutputData1 + ' [lbf]'
                break
        plt.figure(figsize=(12, 10.5))
        contour1 = plt.contour(xx, yy, Results, levels = 25, colors = 'peru')
        plt.clabel(contour1, inline=True, fontsize=8)
        if OutputNumber == 'Two Outputs':
            # Determine what unit to assign to y axis label
            for line in (OutputData2.split()):
                if line == 'Weight':
                    OutputData2withUnits = OutputData2 + ' [lbf]'
                    break
            contour2 = plt.contour(xx, yy, Results2, levels = 10, colors = 'blue') 
            plt.clabel(contour2, inline=10, fontsize=8)
            h1,l1 = contour1.legend_elements()
            h2,l1 = contour2.legend_elements()
            plt.legend([h1[0], h2[0]], [OutputData1withUnits, OutputData2withUnits])
            plt.title('Effect of ' + VariableNamenoUnits + ' and ' + VariableName2noUnits + ' on ' + OutputData1 + ' and ' + OutputData2 )
        elif OutputNumber == 'Single Output':
            h1,l1 = contour1.legend_elements()
            plt.legend([h1[0]], [OutputData1withUnits])
            plt.title('Effect of ' + VariableNamenoUnits + ' and ' + VariableName2noUnits + ' on ' + OutputData1)
        plt.xlabel(VariableName)
        plt.ylabel(VariableName2)      
        print(Results) 
plt.show()