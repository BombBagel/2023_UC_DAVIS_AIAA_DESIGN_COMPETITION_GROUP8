import numpy as np
import pandas as pd
from fpdf import FPDF
import os

# -----------------Set up-------------------
# Save directory of where codes are stored
mainpath = os.getcwd()
# Used to change to main directory after startup
os.chdir(mainpath)
# Read excel file containing aircraft data
df1 = pd.read_excel('Aircraft Data.xlsx')
# Extract data from dataframe
j = 0
AR = df1["AR"][j] # Aspect ratio
pilots = df1["Pilots"][j]
passengers = df1["Passengers"][j]
S_w = df1["Wing Area [ft^2]"][j]
t_c_root = df1["t_c nondimensional thickness of wing at root"][j]
taper_w = df1["Wing Taper Ratio"][j]
sweep_w = df1["Wing Sweep Angle [deg]"][j]*np.pi/180
S_csw = df1["Wing Control Surface Area (including flaps) [ft^2]"][j]
N_z = df1["Limit Load Factor"][j]*1.5
S_ht = df1["Horizontal Tail Area [ft^2]"][j]
L_vt = df1["Wing Quarter MAC to Tail Quarter MAC [ft]"][j]
sweep_ht = df1["Horizontal Tail Sweep Angle [deg]"][j]*np.pi/180
AR_h = df1["Horizontal Tail AR"][j]
S_e = df1["Elevator Area [ft^2]"][j]
S_vt = df1["Vertical Tail Area [ft^2]"][j]
sweep_vt = df1["Vertical Tail Sweep Angle [deg]"][j]*np.pi/180
AR_v = df1["Vertical Tail AR"][j]
S_f = df1["Fuselage Wetted Area [ft]"][j]
l_f= df1["Fuselage Length [ft]"][j]
L = df1["Fuslage Structural Length [ft]"][j] # Assume it is 97% of fuselage length
B_w = df1["Wing Span [ft]"][j]
t = df1["Fuselage Structural Depth [ft]"][j]
N_l = df1["Ultimate Landing Load Factor"][j]
l_m = df1["Length of Extended Main Landing Gear [in]"][j]
N_mw = df1["Number of Main Wheels"][j]
N_mss = df1["Number of Main Gear Shock Struts"][j]
l_n = df1["Length of Extended Nose Gear [in]"][j]
N_nw = df1["Number of Nose Wheels"][j]
N_Lt = df1["Nacelle Length [ft]"][j]
N_w = df1["Nacelle Width [ft]"][j]
S_n = df1["Nacelle Wetted Area [ft^2]"][j]
N = df1["No. Engines"][j]
L_ec = df1["Routing Distance from Engine Front to Cockpit [ft]"][j]
W_en = df1["Engine Weight"][j]
V_t = df1["Total Volume of Fuel [gal]"][j]
V_i = df1["Integral Tanks Volume [gal]"][j]
V_p = df1["Self Sealing Tanks Volume [gal]"][j]
N_t = df1["Number of Tanks"][j]
N_f = df1["Number of Separate Functions Performed by Surface Controls"][j]
W_APU = df1["Uninstalled APU Weight"][j]
R_kva = df1["System Electrical Rating [kVa]"][j]
L_a = df1["Electrical Routing Distance from Generators to Avionics [ft]"][j]
N_gen = df1["Number of Generators"][j]
# Attendants will be expressed as a function of crew since we need 1 for every 50 passengers
attendants = np.ceil(passengers/50)

def calcCompWeight(W_dg, W_L, W_fuel, W_PowerTrain, W_bat):
      # initialize text file to write component weight data
      file3 = open('text3.txt', 'w')

      #
      # Wing weight
      #
      W_w = 0.0051*(W_dg*N_z)**0.557*(S_w)**0.649*AR**0.5*(t_c_root)**(-0.4)*(1+taper_w)**0.1*(np.cos(sweep_w))**(-1.0)*S_csw**0.1
      print('Wing weight is:', W_w, file = file3)

      K_y = 0.3*L_vt          # Aircraft pitching radius of gyration (formula suggested by Raymer)
      #
      # Horizontal tail weight
      #
      W_ht = 0.0379*1*(1)**(-0.25)*W_dg**0.639*N_z**0.10*S_ht**0.75*L_vt**-1.0*K_y**0.704*(np.cos(sweep_ht))**-1.0*AR_h**0.166*(1+S_e/S_ht)**0.1
      print('Horizontal tail weight is:', W_ht, file = file3)

      K_z = L_vt                # Aircraft yawing radius of gyration 
      #
      # Vertical tail weight 
      #
      W_vt = 0.0026*(1+1)**0.225*W_dg**0.556*N_z**0.536*L_vt**-0.5*S_vt**0.5*K_z**0.875*(np.cos(sweep_vt))**-1*AR_v**0.35*t_c_root**-0.5
      print('Vertical tail weight is:', W_vt, file = file3)

      K_ws = 0.75*(1+2*taper_w)/(1+taper_w)*(B_w/L)*np.tan(sweep_w) # Wing sweep factor
      #
      # Fuselage weight
      #
      W_fus = 0.3280*1.06*1.12*(W_dg*N_z)**0.5*L**0.25*S_f**0.302*(1+K_ws)**0.04*(L/t)**0.1
      print('Fuselage weight is:', W_fus, file = file3)

      V_stall = 238/1.3 # Stall speed 
      #
      # Main landing gear weight
      #
      W_mlg = 0.0106*1*W_L**0.888*N_l**0.25*l_m**0.4*N_mw**0.321*N_mss**-0.5*V_stall**0.1
      print('Main landing gear weight is:', W_mlg, file = file3)

      #
      # Nose landing gear weight
      #
      W_nlg = 0.032*1*W_L**0.646*N_l**0.2*l_n**0.5*N_nw**0.45
      print('Nose landing gear weight is:', W_nlg, file = file3)

      K_p = 1.4
      K_tr = 1
      W_ec = 2.331*W_en**0.901*K_p*K_tr
      #
      # Nacelle group weight
      #
      W_ng = 0.6724*1*N_Lt**0.1*N_w**0.294*N_z**0.119*W_ec**0.611*N**0.984*S_n**0.224
      print('Nacelle group weight is:', W_ng, file = file3)

      #
      # Engine controls weight
      #
      W_eCon = 5.0*N + 0.8*L_ec
      print('Engine controls weight is:', W_eCon, file = file3)

      #
      # Starter weight equation
      #
      W_start = 49.19*(N*W_en/1000)**0.541
      print('Starter weight is:', W_start, file = file3)

      #
      # Fuel system weight equation 
      #
      W_fs = 2.405*V_t**0.606*(1+V_i/V_t)**(-1.0)*(1+V_p/V_t)*N_t**0.5
      print('Fuel system weight is:', W_fs, file = file3)

      #N_m = 0                       # Number of surface controls driven by mechanical actuation instead of hydraulics
      #S_cs = 8.38*3.779 + S_e + S_csw # Total control surface area
      #I_yaw = 123                     # Yawing moment of inertia
      #
      # Flight control weight equation
      #
      #W_fc = 145.9*N_f**0.554*(1+N_m/N_f)**(-1.0)*S_cs**0.2*(I_yaw*10**-6)**0.07
      W_fc = 0.64*(W_dg)**(2/3)*1.2
      print('Flight control weight is:', W_fc, file = file3)

      #
      # Installed APU weight equation
      #
      W_APUin = 2.2*W_APU
      print('Installed APU weight is:', W_APUin, file = file3)

      N_c = pilots + attendants # Number of crew
      #
      # Instruments weight equation
      #
      W_inst = 4.509*K_p*0.793*N_c**0.541*N*(l_f+B_w)**0.5
      print('Instrument weight is:', W_inst, file = file3)

      #
      # Hydraulics weight equation
      #
      W_hyd = 0.2673*N_f*(l_f+B_w)**0.937 
      print('Hydraulics weight is:', W_hyd, file = file3)

      #
      # Electrical weight equation
      #
      W_ele = 7.291*R_kva**0.782*L_a**0.346*N_gen**0.1
      print('Electrical weight is:', W_ele, file = file3)

      W_uav = 1000 # Uninstalled avionics weight (rough estimate)
      #
      # Avionics weight equation
      #
      W_avio = 1.73*W_uav**0.983
      print('Avionics weight is:', W_avio, file = file3)


      W_c = 30*3 + 40*50 # Max cargo weight
      #
      # Furnishing weight equation
      #
      #W_fur = 0.0577*N_c**0.1*W_c**0.393*S_f**0.75
      W_fur = 0.211*(W_dg-3057)**0.91
      print('Furnishing weight is:', W_fur, file = file3)

      #
      # Seat weight 
      #
      #W_seat = 39*25 + 14
      #print('Seat weight is:', W_seat)

      N_pers = pilots + passengers + attendants # Number of crew and passengers
      V_pr = 2370 # Volume of pressurized section
      #
      # Air conditioning weight
      #
      W_ac = 62.36*N_pers**0.25*(V_pr/1000)**0.604*W_uav**0.10
      print('Air conditioning is:', W_ac, file = file3)

      #
      # Anti ice weight equation
      #
      W_dIce = 0.002*W_dg
      print('Anti ice weight is:', W_dIce, file = file3)

      #
      # Propeller Weight
      #
      #W_prop = 24*(6)**0.391*(13)*(5000/1000)**0.782
      W_prop = 0.108*(13*5000*6**0.5)**0.782*2
      print('Propeller weight is:', W_prop, file = file3)

      print('Total weight from above is:', W_w + W_ht + W_vt + W_fus + W_mlg + W_nlg + W_ng + W_eCon + W_start + W_fs + W_fc + W_APUin + W_inst + W_hyd\
            + W_ele + W_avio + W_fur + W_ac + W_dIce, file = file3)

      print('Payload weight: 12660', file = file3)
      print('Battery weight:', W_bat, file = file3)
      print('Powertrain weight:', W_PowerTrain, file = file3)
      print('Fuel weight:', W_fuel, file = file3)
      print('Total weight is:', W_w + W_ht + W_vt + W_fus + W_mlg + W_nlg + W_ng + W_eCon + W_start + W_fs + W_fc + W_APUin + W_inst + W_hyd\
            + W_ele + W_avio + W_fur + W_ac + W_dIce + 12660 + W_bat + W_PowerTrain + W_fuel + W_prop, file = file3)

      # close text files
      file3.close()

      # Convert text files to pdfs
      #
      #
      pdf4 = FPDF()
      pdf4.add_page()
      pdf4.set_font("Arial", size = 11)
      f = open("text3.txt", "r")
      for x in f:
            pdf4.cell(200, 10, txt = x, ln = 1, align = 'C')
      pdf4.output("pdf4.pdf") 
      f.close()

      return (W_w, W_ht, W_vt, W_fus, W_mlg, W_nlg, W_ng, W_eCon, W_start, W_fs, W_fc, W_APUin, W_inst, W_hyd,
            W_ele, W_avio, W_fur, W_ac, W_dIce, W_prop)

