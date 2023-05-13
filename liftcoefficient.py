import numpy as np
kts_to_mps = 0.514444
MTOW = 64482*4.44822 #N

rho_SL = 1.23 #kg/m3
rho_28000 = 0.496 #kg/m3
S = 730*0.092903 #m^2

v_cruise = 275*kts_to_mps #kts

#cl_TO = MTOW/(0.5*rho_SL*S*77.1667**2)  # kg m s-2 / kg m-3 m2 m2 s-2
#cl_Land =  MTOW/(0.5*rho_SL*S*v_Land**2)

clmax_TO = 2
clmax_cruise = 1.55
clmax_appoach = 2.6
cl_TO = clmax_TO*(1/1.13)**2
cl_cruise =  MTOW/(0.5*rho_28000*S*v_cruise**2) 
cl_approach = clmax_appoach*(1/1.23)**2

print("Cl take off: " +str(cl_TO))
print("Cl cruise: " +str(cl_cruise))
print("Cl approach: " +str(cl_approach))

v_to = np.sqrt(MTOW/(cl_TO*.5*rho_SL*S))*3.28084


print(v_to)



