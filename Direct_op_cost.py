import numpy as np
import os
import pandas as pd
from fpdf import FPDF

# -----------------Set up-------------------
# Save directory of where codes are stored
mainpath = os.getcwd()
# Used to change to main directory after startup
os.chdir(mainpath)
# Create folder to store results
if not os.path.exists('Results'):
   os.makedirs('Results')
# Read excel file containing aircraft data
df1 = pd.read_excel('Aircraft Data.xlsx')
j = 0
range = df1["Range \n(nmi)"][j]

# Direct Operating Cost (DOC) Calculations
# NOTE: RFP wants DOC per airplane flight hour
    
x = 1   # NOTE: These are variables that need to be determined
        # Or if variable already exist, then needs to be replaced!

# NOTE: Depending on the needed inputs, unit conversions may need to be added/adjusted

def getOperatingCost(W_fuel,Bat_Energy,W_0,Unit_Cost,C_Bat):
    # initialize text file to write operating cost data
    Direct_op_costfile = open('Direct_op_cost.txt', 'w')

    # Crew Cost [Finger et al.]
    R_crew = 60 # Flight crew cost rate ($/hr)
    N_crew = 2 
    t = 2.25   # Flight time (hr)
    C_crew = R_crew*N_crew*1.5*t
    print('Crew Cost is: ', C_crew/t, file=Direct_op_costfile)

    # Attendents cost [Metabook, Liebeck et al.]
    def getCEF(by,ty): #base year, then year
        # CEF must be determined for each different base year!
        t_CEF = 5.17053+0.104981*(by-2006)  
        b_CEF = 5.17053+0.104981*(by-2006)
        CEF = t_CEF/b_CEF   # Cost escalation factor  
        return CEF
    CEF = getCEF(1999,2035) 
    C_attd = 60*1*CEF*t # only 1 attendant 
    print('Attendents Cost is: ', C_attd/t, file=Direct_op_costfile)

    # Fuel Cost [Metabook, Kroo]
    P_f = 6.75   # Price per gallon of jet fuel ($/gal)
    rho_f = 6.75 # Density of jet fuel (lb/gal)
    C_fuel = 1.02*W_fuel*P_f/rho_f 
    print('Fuel Cost is: ', C_fuel/t, file=Direct_op_costfile)

    # Electricity Cost [Finger et al.]
    T_E = 0.40 # Price per kWh  ($/kWh)
    eta_charge = 0.9    # Charging efficiency
    C_E = T_E*Bat_Energy/eta_charge/1000
    print('Electricity Cost is: ', C_E/t, file=Direct_op_costfile)

    # Oil Cost [Metabook, not sure who, I assume Kroo?]
    W_oil = 0.0125*W_fuel*t/100    # Weight of oil (lb)
    P_o = 5 # Price per gallon of oil ($/gal)
    rho_o = 7.5   # Density of oil (lb/gal)
    C_oil = 1.02*W_oil*P_o/rho_o
    print('Oil Cost is: ', C_oil/t, file=Direct_op_costfile)

    # Landing Fees [Metabook, Liebeck et al.]
    C_airport = CEF*1.5*W_0/1000  
    print('Landing Fees Cost is: ', C_airport/t, file=Direct_op_costfile) 

    # Navigation Fees [Metabook, Liebeck et al.]
    CEF = getCEF(1989,2035)
    C_nav = 0.5*CEF*(1.852*range/t)*np.sqrt(0.00045359237*W_0/50)
    print('Navigation Fees Cost is: ', C_nav/t, file=Direct_op_costfile)

    # Maintenance Cost [Finger et al.]
    R_ap = 60   # Rate of certified mechanic ($/hr)
    F_mf = 0.2  # ratio between flight hours and maintenance hours
    C_maint = R_ap*F_mf*t
    print('Maintenance Cost is: ', C_maint/t, file=Direct_op_costfile)

    # Aircraft Depreciation Cost [Metabook] 
    U = (1.510**3)*(3.4546*t + 2.994 - (((12.289*(t**2))-(5.6626*t)+8.964)**0.5))    
    # Annual utilization (hr) [Roskam, part VIII] 
    K_dep = 0.1 
    # Depreciation Factor 
    years = 20 # Number of years in service 
    C_ACD = Unit_Cost*(1-K_dep)*t/(years*U)
    # Aircraft Depreciation Cost [Finger et al.]
    #n_flights = 1000   # Number of flights over the depreciation period
    #C_ACD = Unit_Cost/n_flights
    print('Aircraft Depreciation Cost is (Fixed Cost): ', C_ACD/(years*U), file=Direct_op_costfile)

    # Battery Depreciation Cost [Finger et al.]
    n_cyc = 1000    # Assumed life cycle of 1000
    C_BatD = C_Bat*Bat_Energy/(n_cyc*Bat_Energy*(1-((1-(0.2))/2))) # 0.2 represents percentage of (E_bat end of life/ E_bat)
    print('Battery Depreciation Cost is (Fixed Cost): ', C_BatD/(years*U), file=Direct_op_costfile)

    # Aircraft Insurance Cost [Finger et al.]
    C_ins = 500+(0.015*Unit_Cost)
    print('Aircraft Insurance Cost is (Fixed Cost): ', C_ins/(years*U), file=Direct_op_costfile)

    TotalDOC = (C_crew + C_attd + C_fuel + C_E + C_oil + C_airport + C_nav + C_maint + C_ACD + C_BatD + C_ins)
    print('Total DOC is: ', TotalDOC, file=Direct_op_costfile)

    # Close text files
    Direct_op_costfile.close()

    # Convert NP_CG_SMfile to pdf
    #
    #
    Direct_op_costfilepdf = FPDF()
    Direct_op_costfilepdf.add_page()
    Direct_op_costfilepdf.set_font("Arial", size = 11)
    f = open("Direct_op_cost.txt", "r")
    for x in f:
        Direct_op_costfilepdf.cell(200, 10, txt = x, ln = 1, align = 'C')
    Direct_op_costfilepdf.output("Direct_op_cost.pdf") 
    f.close()

