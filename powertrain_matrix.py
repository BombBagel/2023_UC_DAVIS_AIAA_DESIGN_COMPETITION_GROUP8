import numpy as np
import matplotlib.pyplot as plt
import sys
 
# available states: propulsion, charge-both-propulsion, charge-no-P1-propulsion, charge-no-P2-propulsion
# available architectures: serial-parallel, serial, parallel, conventional
# supplied power ratio is the hybridzation ratio
def getpowertrain(architecture,state,supplied_power_ratio):
    component = ['GT','GB','P1','EM1','PM','EM2','P2','SuppliedPR','ShaftPR','P'] # components of the architecture
    variables = ['P_f','P_gt','P_gb','P_s1','P_e1','P_bat','P_e2','P_s2','P_p1','P_p2'] # calculated power values
    
    # efficiencies and constants
    n_GT = 0.40
    n_GB = 0.96
    n_P1 = 0.9
    n_EM1 = 0.96
    n_PM = 0.99
    n_EM2 = 0.96
    n_P2 = 0.85
    shaft_power_ratio = 0.9

    # unchanged component equations
    GT = [-n_GT,1,0,0,0,0,0,0,0,0]
    P1 = [0,0,0,-n_P1,0,0,0,0,1,0]
    EM2 = [0,0,0,0,0,0,-n_EM2,1,0,0]
    P2 = [0,0,0,0,0,0,0,-n_P2,0,1]
    Shaft_Ratio = [0,0,0,-shaft_power_ratio,0,0,0,(1-shaft_power_ratio),0,0]
    Supplied_Ratio = [-supplied_power_ratio,0,0,0,0,(1-supplied_power_ratio),0,0,0,0]
    P = [0,0,0,0,0,0,0,0,1,1]
    
    # changes component equations depending on state/architecture
    if state == 'propulsion':
        # propulsion
        if (architecture == 'serial-parallel') or (architecture == 'serial') or (architecture == 'conventional'):
            # power comes from battery to P2 and power comes from gearbox to P2
            GB = [0,-n_GB,1,1,0,0,0,0,0,0]
            EM1 = [0,0,-n_EM1,0,1,0,0,0,0,0]
            PM = [0,0,0,0,-n_PM,-n_PM,1,0,0,0]
        elif architecture == 'parallel':
            # power comes from battery and fuel to P1
            GB = [0,-n_GB,-n_GB,1,0,0,0,0,0,0]
            EM1 = [0,0,1,0,-n_EM1,0,0,0,0,0]
            PM = [0,0,0,0,1,-n_PM,1,0,0,0]
        else:
            sys.exit("unknown architecture")
    elif state == 'charge-both-propulsion':
        # battery recharging with propulsion of P1 and P2
        GB = [0,-n_GB,1,1,0,0,0,0,0,0]
        EM1 = [0,0,-n_EM1,0,1,0,0,0,0,0]
        PM = [0,0,0,0,-n_PM,1,1,0,0,0]
    elif state == 'charge-no-P2-propulsion':
        # battery recharging but no electrical propulsion (P2)
        GB = [0,-n_GB,1,1,0,0,0,0,0,0]
        EM1 = [0,0,-n_EM1,0,1,0,0,0,0,0]
        PM = [0,0,0,0,-n_PM,1,1,0,0,0]
        removeindex = ['EM2','P2','ShaftPR']
        removevars = ['P_e2','P_s2','P_p2']
    elif state == 'charge-no-P1-propulsion':
        # battery recharging but no fuel propulsion (P1)
        GB = [0,-n_GB,1,1,0,0,0,0,0,0]
        EM1 = [0,0,-n_EM1,0,1,0,0,0,0,0]
        PM = [0,0,0,0,-n_PM,1,1,0,0,0]
        removeindex = ['P1','ShaftPR']
        removevars = ['P_s1','P_p1']
    else:
        sys.exit("unknown state")
    A = np.array([GT,GB,P1,EM1,PM,EM2,P2,Supplied_Ratio,Shaft_Ratio,P])
    
    # decided which variables are in each architecture
    if architecture == 'serial-parallel':
        vars = ['P_f','P_gt','P_gb','P_s1','P_e1','P_bat','P_e2','P_s2','P_p1','P_p2']
        index = ['GT','GB','P1','EM1','PM','EM2','P2','ShaftPR','SuppliedPR','P']
    elif architecture == 'serial':
        vars = ['P_f','P_gt','P_gb','P_e1','P_bat','P_e2','P_s2','P_p2']
        index = ['GT','GB','EM1','PM','EM2','P2','SuppliedPR','P']
    elif architecture == 'parallel':
        vars = ['P_f','P_gt','P_gb','P_s1','P_e1','P_bat','P_p1']
        index = ['GT','GB','P1','EM1','PM','SuppliedPR','P']
    elif architecture == 'conventional':
        vars = ['P_f','P_gt','P_s1','P_p1']
        index = ['GT','GB','P1','P']
    else:
        sys.exit("unknown architecture")

    # remove variables and commponents if needed depending on state
    if 'removeindex' in locals(): # checks if 'removeindex' is a variable
        for i,v in zip(removeindex,removevars):
            try:
                index.remove(i)
            except ValueError:
                pass
            try:
                vars.remove(v)
            except ValueError:
                pass

    # extracted the wanted variables
    rowkeep = index.copy()
    columnkeep = vars.copy()
    i = -1
    for ind, var in zip(index,vars):
        i += 1
        rowkeep[i] = component.index(ind)
        columnkeep[i] = variables.index(var)
    A = A[rowkeep,:] # rows
    A = A[:,columnkeep] # columns
    return A,vars


'''## example inputs and outputs
# serial-parallel powertrain
architecture,state,supplied_power_ratio = 'serial-parallel','charge-both-propulsion',0.2
A,vars = getpowertrain(architecture,state,supplied_power_ratio)
print('\n----------------------------------------------------')
print('{}, {}'.format(architecture,state))
print('----------------------------------------------------')
print("A = \n{}".format(A))
print("\nvariables = \n{}".format(vars))


# serial powertrain
architecture,state,supplied_power_ratio = 'serial','charge-no-P1-propulsion',0.2
A,vars = getpowertrain(architecture,state,supplied_power_ratio)
print('\n----------------------------------------------------')
print('{}, {}'.format(architecture,state))
print('----------------------------------------------------')
print("A = \n{}".format(A))
print("\nvariables = \n{}".format(vars))


# parallel powertrain
architecture,state,supplied_power_ratio = 'parallel','charge-no-P2-propulsion',0.2
A,vars = getpowertrain(architecture,state,supplied_power_ratio)
print('\n----------------------------------------------------')
print('{}, {}'.format(architecture,state))
print('----------------------------------------------------')
print("A = \n{}".format(A))
print("\nvariables = \n{}".format(vars))


# conventional powertrain
architecture,state,supplied_power_ratio = 'conventional','propulsion',0.2
A,vars = getpowertrain(architecture,state,supplied_power_ratio)
print('\n----------------------------------------------------')
print('{}, {}'.format(architecture,state))
print('----------------------------------------------------')
print("A = \n{}".format(A))
print("\nvariables = \n{}".format(vars))'''

