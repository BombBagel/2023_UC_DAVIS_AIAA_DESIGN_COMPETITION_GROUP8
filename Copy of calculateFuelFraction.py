import numpy as np

# -----------------Parameters-------------------
eta = 0.8 # Propulsive efficiency
P = 9670 # Takeoff power in horsepower
W0 = 77363 # Takeoff weight in lbf

# -----------------Main Code--------------------
# Calculate warm up and takeoff fuel fraction
t = 60 
W1_W0 = (1-t*0.5/3600/eta*P/W0)*(1-t*15*0.3/3600/eta*P*0.05/W0)
print(W1_W0)

# Calculate climb fuel fraction
W = W1_W0*W0 
S = 923
e = 0.85
AR = 13
C_d_0 = 0.0156*2
rho = 0.0009567
rho_SL = 0.00237
# Find best rate of climb
V = np.sqrt(2*W/rho_SL/S*np.sqrt(1/(np.pi*e*AR)/3/C_d_0))
# Specify number of segments
i = 50
delta_h = 28000/i
D = rho*V**2/2*S*C_d_0
Wi1_Wi = np.exp(-0.5/3600*delta_h/(V*(1-D/(P*550/eta/V))))
for j in range(i):
    W = Wi1_Wi*W0 
    # Find best rate of climb
    V = np.sqrt(2*W/rho_SL/S*np.sqrt(1/(np.pi*e*AR)/3/C_d_0))
    # Specify number of segments
    i = 10000
    delta_h = 28000/i
    D = rho*V**2/2*S*C_d_0
    Wi1_Wi = np.exp(-0.5/3600*delta_h/(V*(1-D/(P*550/eta/V))))*Wi1_Wi
W2_W1 = Wi1_Wi
print(W2_W1)
W_after_cruise = W2_W1*W1_W0*W0 
W3_W2 = (W_after_cruise-1050)/W_after_cruise
print(W3_W2)

W4_W3 = 0.990

# Loiter fuel fraction
E = 0.25*3600
V = 250*1.688
L_D = 17
W5_W4 = np.exp((-E*((0.5*V) / (550*eta*3600))) / (0.866*L_D))

W6_W5 = 0.995

fuelWeight_taxi_takeoff_climb = (1-W2_W1*W1_W0)*W0*0.8
print(fuelWeight_taxi_takeoff_climb)
fuelWeight_Descent_Loiter_Landing = -(W6_W5*W5_W4*W4_W3*W3_W2*W2_W1*W1_W0 - W4_W3*W3_W2*W2_W1*W1_W0)*W0*0.8
print(fuelWeight_Descent_Loiter_Landing)
totalFuel = fuelWeight_Descent_Loiter_Landing + fuelWeight_taxi_takeoff_climb + 1050
print(totalFuel)