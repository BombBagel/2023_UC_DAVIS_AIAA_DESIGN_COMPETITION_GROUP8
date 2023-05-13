import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
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
TipChord = df1['Wing Tip Chord [ft]'][j]
RootChord = df1['Wing Root Chord [ft]'][j]
wingTaperRatio = TipChord/RootChord
Alpha = df1['Wing Sweep Angle [deg]'][j] * np.pi/180 # Wing sweep (Average over sections 1, 2, 3, neglect 4: 22+15+20/3)
S_w = df1['Wing Area [ft^2]'][j] # Wing reference area [ft^2]
AR = df1['AR'][j] # Wing AR
b = df1['Wing Span [ft]'][j]  # Wingspan [ft] (used projected span due to trapezoidal assumption)
c = 2/3*RootChord*((1+wingTaperRatio+wingTaperRatio**2)/(1+wingTaperRatio)) # Wing Mean aerodynamic chord
Y = (b/6)*((1+(2*wingTaperRatio))/(1+wingTaperRatio))   # Spanwise location of Wing MAC
# Horizontal Tail
TipChord_h = df1['Horizontal Tail Tip Chord [ft]'][j] 
RootChord_h = df1['Horizontal Tail Root Chord [ft]'][j]
H_TaperRatio = TipChord_h/RootChord_h
S_h = df1['Horizontal Tail Area [ft^2]'][j] 
Alpha_h = df1['Horizontal Tail Sweep Angle [deg]'][j] * np.pi/180
AR_h = df1['Horizontal Tail AR'][j] # Horizontal tail AR 
b_h = df1['Horizontal Tail Span [ft]'][j] # Htail span [ft]
c_h = 2/3*RootChord_h*((1+H_TaperRatio+H_TaperRatio**2)/(1+H_TaperRatio)) # Htail Mean aerodynamic chord
Y_h = (b_h/6)*((1+(2*H_TaperRatio))/(1+H_TaperRatio))   # Spanwise location of Htail Wing MAC
# Vertical Tail
TipChord_v = df1['Vertical Tail Tip Chord [ft]'][j]  # Vertical tail [ft]
RootChord_v = df1['Vertical Tail Root Chord [ft]'][j]
V_TaperRatio = TipChord_v/RootChord_v  
Alpha_v = df1['Vertical Tail Sweep Angle [deg]'][j] * np.pi/180 # Vertical tail sweep
AR_v = df1['Vertical Tail AR'][j] # Vertical tail AR 
b_v = df1['Vertical Tail Span [ft]'][j] # Vtail span [ft]
c_v = 2/3*RootChord_v*((1+V_TaperRatio+V_TaperRatio**2)/(1+V_TaperRatio)) # Vtail Mean aerodynamic chord
Y_v= (b_v/6)*((1+(2*V_TaperRatio))/(1+V_TaperRatio))   # Spanwise location of Vtail MAC
x_wLE = df1['Distance from nose to wing leading edge [ft]'][j]
x_hLE = df1['Distance from nose to horizontal tail leading edge [ft]'][j]
x_vLE = df1['Distance from nose to vertical tail leading edge [ft]'][j] # Distance from nose to vertical tail leading edge [ft] (obtained from openvsp)

def NP_CG_SM(v,W_w,W_h,W_vt,W_f,W_mlg,W_nlg,W_n,W_eCon,W_start,W_fs,W_fc,W_APUin,W_inst,
             W_hyd,W_ele,W_avio, W_fur, W_ac, W_dIce, W_prop, W_fuel, W_bat, W_pt):
    def getNeutralPoint(v):
        Ta = 418.8 # Atmospheric pressure at cruise altitude

        l_h = x_hLE + (0.25*c_h) + (Y_h*np.tan(Alpha_h)) - x_cg_0  # Distance of horizontal tail AC from CG (appromximated geometrically assuming wing is trapezoidal and uniform, Raymer)
        x_ac_w = x_wLE + 0.25*c + Y*np.tan(Alpha)          # Distance of Wing AC from nose 
        print('The wing aerodynamic center from the nose is:', round(x_ac_w,2), file = NP_CG_SMfile)
        a_a = np.sqrt(1.4*53.35*32.17*Ta) 
        M = v/a_a
        WingLiftDerivative = (2*np.pi*AR) / (2+np.sqrt((AR/eta)**2*(1+(np.tan(Alpha))**2-M**2)+4))
        TailCleanLiftDerivative = (2*np.pi*AR_h) / (2+np.sqrt((AR_h/eta)**2*(1+(np.tan(Alpha_h))**2-M**2)+4))
        TailLiftDerivative = TailCleanLiftDerivative*(1 - (2/np.pi/AR)*WingLiftDerivative)
        FuselagePitchingMoment = Kf*d_f**2*l_f/S_w/c*(WingLiftDerivative)**(-1)
        Xnp = ((l_h*S_h/c/S_w)*TailLiftDerivative*(WingLiftDerivative)**(-1)-FuselagePitchingMoment)*c #Neutral point, distance from the Wing AC
        Xnp_0 = x_ac_w+Xnp    # Distance of NP from the nose
        x_cg = x_cg_0-x_ac_w  # reference Location of CG from the Wing AC
        x_ac_h = x_cg_0 + l_h # reference of htail AC from the nose
        print('The neutral point from the nose is: ', round(Xnp_0,2), file = NP_CG_SMfile)
        # print('x_cg from wing AC = ', x_cg, file = NP_CG_SMfile)
        # print('x_ac of wing = ', x_ac, file = NP_CG_SMfile)
        # print('x_ac_h = ', x_ac_h, file = NP_CG_SMfile)
        return Xnp_0
    
    # initialize text file to write NP_CG_SM data
    NP_CG_SMfile = open('NP_CG_SM.txt', 'w')
     
    # Insert variables here
    eta = 0.97
    Kf = 0.487 # Empirical factor depending on wing quarter chord position (Change base on x_ac/l_f)
    d_f = 9.1 # Fuselage diameter [ft]
    l_f = 83 # Fuselage length [ft]

    print('Wing MAC =',c, file = NP_CG_SMfile)

    ##########[ X_cg Calculations ]###################################################################################
    x_cg_w =  x_wLE +  Y*np.tan(Alpha) + (0.4*c)  #(ft) cg of wing from nose
    print('Wing center of gravity = ', x_cg_w, file = NP_CG_SMfile)

    x_cg_h =  x_hLE +  Y_h*np.tan(Alpha_h) + (0.4*c_h)  #(ft) cg of htail from nose
    print('Horizontal tail center of gravity = ', x_cg_h, file = NP_CG_SMfile)

    x_cg_v =  x_vLE +  Y_v*np.tan(Alpha_v) + (0.4*c_v)  #(ft) cg of vtail from nose
    print('Vertical tail center of gravity = ', x_cg_v, file = NP_CG_SMfile)

    x_cg_f = 0.45*l_f #fuselage cg is 40-50% of length
    print('Fuselage center of gravity = ', x_cg_f, file = NP_CG_SMfile)

    x_cg_mlg = 45   #Located at centroid of landing gear
    print('Main landing gear center of gravity = ', x_cg_mlg, file = NP_CG_SMfile)

    x_cg_nlg = 5.5435   #Located at centroid of landing g
    print('Nose landing gear center of gravity = ', x_cg_nlg, file = NP_CG_SMfile)

    x_cg_n = 25.4955+2.5 + (0.4*14.8333)   #cg of nacelle from nose
    print('Nacelle center of gravity = ', x_cg_n, file = NP_CG_SMfile)

    # x_cg of above weights will use "All-else empty" distance
    x_cg_AEE = 0.45*l_f #CG location of "All-else empty"
    print('"All else empty" center of gravity = ', x_cg_AEE, file = NP_CG_SMfile)

    #x_cg_fc is same as wing
    print('Flight controls center of gravity = ', x_cg_w, file = NP_CG_SMfile)

    x_cg_avi = 4
    print('Avionics center of gravity = ', x_cg_avi, file = NP_CG_SMfile)

    x_cg_apu = 80
    print('APU center of gravity = ', x_cg_apu, file = NP_CG_SMfile)

    x_cg_prop = 25.487+2.5  #cg of propeller
    print('Propeller center of gravity = ', x_cg_prop, file = NP_CG_SMfile)

    p = 0  #Percentage of passenger baggage that is loaded into the cargo containers. 
    W_payload = 12220 -(p*(50*40)) #Payload weight (Total weight of payload (200+40)*50 passengers + (190+30)*3 - cockpit/pilot crew - baggage stored in the rear)
    x_cg_payload = 37.4   #Location of seats CG
    print('Payload center of gravity = ', x_cg_payload, file = NP_CG_SMfile)

    W_baggage = p*(50*40)    #Weight of passenger baggage that is loaded into the cargo containers.
    x_cg_baggage = (58.5+65.2)/2

    W_pilot = 380+60    #Pilots weight and baggage
    x_cg_pilot = 8      #Pilot seat CG
    print('Cockpit crew center of gravity = ', x_cg_pilot, file = NP_CG_SMfile)

    x_cg_bat = ((17.5+(12/2))+(48+(10/2)))/2 #Average cg of battery 1 and 2
    print('Battery center of gravity = ', x_cg_bat, file = NP_CG_SMfile)

    x_cg_pt = 28 + ((3.3/2)+(6.992/2))/2 #cg of powertrain
    print('Powertrain center of gravity = ', x_cg_pt, file = NP_CG_SMfile)

    #x_cg_fuel is same as wing
    print('Fuel center of gravity = ', x_cg_w, file = NP_CG_SMfile)


    # x_MAC = x_wLE + Y*np.tan(Alpha) #ignore
    # print(c)
    # print(x_MAC)
    # print(Y)

    def getTW():
        # Calculates the total weight
        W_total = (W_w+W_h+W_vt+W_f+W_mlg+W_nlg+W_n+W_bat+W_pt+W_prop+ W_fuel+W_eCon+W_start+
                   W_fs+W_fc+W_APUin+W_inst+W_hyd+W_ele+W_avio+W_fur+W_ac+W_dIce+W_payload+W_baggage+W_pilot)
        print(W_total,file = NP_CG_SMfile)
        return W_total

    ######[ x_cg calculation (Reference Distance of CG from the nose [ft] Note: This is reference from nose!!!!) ]#########################################################
    def getCG():
        # Must get total weight first!
        # Calculates the cg location of the aircraft from the nose
        x_cg_0 = ((W_w*x_cg_w)+(W_h*x_cg_h)+(W_vt*x_cg_v)+(W_f*x_cg_f)+(W_mlg*x_cg_mlg)+(W_nlg*x_cg_nlg)+(W_n*x_cg_n)+(W_bat*x_cg_bat)
        +(W_pt*x_cg_pt)+(W_prop*x_cg_prop)+(W_fuel*x_cg_w)+(W_payload*x_cg_payload)+ (W_baggage*x_cg_baggage) + (W_pilot*x_cg_pilot) 
        +(W_APUin*x_cg_apu)+(W_avio*x_cg_avi)+(W_fc*x_cg_w)+((W_eCon+W_start+W_fs+W_inst+W_hyd+W_ele+W_fur+W_ac+W_dIce)*x_cg_AEE)) /(W_total)
        print('Airplane X_cg from the nose is: ', round(x_cg_0,2), file = NP_CG_SMfile)
        return x_cg_0 

    def getSM():
        # Calculates the static margin
        SM = (Xnp_0-x_cg_0)/c
        print('Static Margin = ', round(SM*100,2),'%', file = NP_CG_SMfile)   #CG must be in front of the neutral point for stability (positive SM = stability, ideally around 5% to 10%)
        return SM

    # Fully loaded case (Empty+Crew+Fuel+Payload)
    print('\n[Fully loaded]', file = NP_CG_SMfile)
    W_total = getTW()
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()

    W_total_a = W_total
    x_cg_0_a = x_cg_0
    Xnp_0_a = Xnp_0
    SM_a = SM

    # Fully loaded case, but 50% of passenger baggage in the rear (Empty+Crew+Fuel+Payload)
    p = 0.5 
    W_payload = 12220 -(p*(50*40))
    W_baggage = p*(50*40)

    print('\n[Fully loaded, 50 percent baggage in rear]', file = NP_CG_SMfile)
    W_total = getTW()
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()

    W_total_a2 = W_total
    x_cg_0_a2 = x_cg_0
    Xnp_0_a2 = Xnp_0
    SM_a2 = SM

    # Fully loaded case, but 100% of passenger baggage in the rear (Empty+Crew+Fuel+Payload)
    p = 1 
    W_payload = 12220 -(p*(50*40))
    W_baggage = p*(50*40)

    print('\n[Fully loaded, 100 percent baggage in rear]', file = NP_CG_SMfile)
    W_total = getTW()
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()

    W_total_a3 = W_total
    x_cg_0_a3 = x_cg_0
    Xnp_0_a3 = Xnp_0
    SM_a3 = SM

    # No payload (Empty+Crew+Fuel)
    W_payload = 0
    W_baggage = 0
    W_pilot = 380+60

    print('\n[Empty+Crew+Fuel]', file = NP_CG_SMfile)
    W_total = getTW()
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()

    W_total_b = W_total
    x_cg_0_b = x_cg_0
    Xnp_0_b = Xnp_0
    SM_b = SM

    # No payload, No Fuel (Empty+Crew)
    W_payload = 0
    W_baggage = 0
    W_pilot = 380+90
    W_fuel = 0

    print('\n[Empty+Crew]', file = NP_CG_SMfile)
    W_total = getTW() #updates the total weight
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()
    W_total_c = W_total
    x_cg_0_c = x_cg_0
    Xnp_0_c = Xnp_0
    SM_c = SM

    #Empty (no crew, no payload, no fuel)
    W_pilot = 0
    W_payload = 0
    W_baggage = 0
    W_fuel = 0

    print('\n[Empty]', file = NP_CG_SMfile)
    W_total = getTW() #updates the total weight
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()
    W_total_d = W_total
    x_cg_0_d = x_cg_0
    Xnp_0_d = Xnp_0
    SM_d = SM

    # Out of fuel (Empty+Crew+Payload)
    p = 0
    W_payload = 12220 -(p*(50*40))
    W_baggage = p*(50*40) 
    W_pilot = 380+60
    W_fuel = 0

    print('\n[Loaded but No fuel]', file = NP_CG_SMfile)
    W_total = getTW() #updates the total weight
    x_cg_0 = getCG()
    Xnp_0 = getNeutralPoint(v)
    SM = getSM()
    W_total_e = W_total
    x_cg_0_e = x_cg_0
    Xnp_0_e = Xnp_0
    SM_e = SM

    X = np.array([x_cg_0_a, x_cg_0_a2, x_cg_0_a3, x_cg_0_e, x_cg_0_b, x_cg_0_c, x_cg_0_d])
    W = np.array([W_total_a, W_total_a2, W_total_a3, W_total_e, W_total_b, W_total_c, W_total_d])

    X_lim_right = (Xnp_0 - 0.05*c)    #CG right hand limit based on 5% SM limit
    #X_lim_left = (X_lim_right - 0.08*c) #CG left hand limit based on 8% of chord away from right hand limit
    X_lim_left = (Xnp_0 - 0.1*c) #CG left hand limit based on 10% SM 
    print(X_lim_right, X_lim_left, file = NP_CG_SMfile)

    # Create pdf to save graph
    pp = PdfPages('Graph2.pdf')
    
    plt.figure(figsize=(10,6))
    plt.plot(X/c,W, marker ='o', markerfacecolor = 'orange')
    plt.axvline(x = X_lim_right/c, linestyle = '--', color = 'gray', label = 'Aft c.g. limit')
    plt.axvline(x = X_lim_left/c, linestyle = '--', color = 'gray', label = 'Forward c.g. limit')
    plt.title('C.G. Location for Weight Scenarios')
    plt.xlabel(r'$x_{cg}/ \bar{c}$')
    plt.ylabel('Gross Weight (lb)')
    labels = ['Fully loaded', 'Fully loaded 2','Fully loaded 3','Out of fuel', 'No payload', 'Crew only', 'Empty']
    for i in range(7):
        plt.annotate(labels[i], xy = (X[i]/c, W[i]))
    plt.text((X_lim_right/c)+0.0003, 56000, "5% Static Margin", rotation=90, color = 'gray')
    plt.text((X_lim_left/c)-0.0022, 56000, "10% Static Margin", rotation=90, color = 'gray')
    plt.savefig('CG Travel Graph.png',bbox_inches="tight")
    pp.savefig(bbox_inches="tight")
    
    # Close pdf files after graphs have been saved
    pp.close()

    # Close text files
    NP_CG_SMfile.close()

    # Convert NP_CG_SMfile to pdf
    #
    #
    NP_CG_SMfilepdf = FPDF()
    NP_CG_SMfilepdf.add_page()
    NP_CG_SMfilepdf.set_font("Arial", size = 11)
    f = open("NP_CG_SM.txt", "r")
    for x in f:
        NP_CG_SMfilepdf.cell(200, 10, txt = x, ln = 1, align = 'C')
    NP_CG_SMfilepdf.output("NP_CG_SM.pdf") 
    f.close()