import os
import subprocess as sp
import numpy as np
import re
import matplotlib.pyplot as plt

cwd = os.getcwd()
cwd = os.path.join(cwd, 'avl')
os.chdir(cwd)
avlfile = os.path.join(cwd, 'avl.exe')


def getAVLFIle(cl,config):

    cl_array = [ ]
    filenames = []
    for i in range(len(cl)):
        cl_array = cl_array + ["a", "c", str(cl[i]),"x", "FT", "Cl" + str(cl[i]) + ".txt"]
        filenames = filenames + ["Cl" + str(cl[i]) + ".txt"]
   


    fileArray = []
    avlProcess = sp.Popen(avlfile,
                        shell=False,
                        stdin=sp.PIPE,
                        stdout=open("AVLsession.log", 'w'),         # Put the output of this terminal into the open log file
                        stderr=sp.PIPE)
    allCommands = ["load",
                config,
                "mass",
                "EP.mass",
                "mset",
                "0",
                "oper",
                "d3 pm", 
                "0",
                ] + cl_array 
    avlProcess.stdin.write(str.encode("\n".join(allCommands)))
    avlProcess.stdin.flush()

    # Close the input
    avlProcess.communicate()
    
    return(filenames)
    
def readData(file):
    datafile = open(file)
    for line in datafile:
        if 'CDind' in line:
            line = re.sub("=", "= ", line)
            line = re.sub('\s+', ' ', line)
            data = np.array(line.split())
            data_loc = np.where(data == "CDind")[0]
            Cdind = float(data[data_loc+2])
    return Cdind

# Specify number of data points
N = 20

# Get data for cruise config
CL_cruise = np.linspace(-0.9,1.5,N) # Input cruise config CL
config = "EPCruiseFinal.avl"
files = getAVLFIle(CL_cruise, config)
CDindArray = np.ones(len(CL_cruise))
for i in range(len(CDindArray)):
    CDindArray[i] = readData(files[i])
    os.remove(files[i])
Cdo_cruise = 0.020352323157863724
plt.plot(CDindArray+Cdo_cruise,CL_cruise)

# Get data for takeoff config
CL_takeoff = np.linspace(-1.2,1.95,N) # Input takeoff config CL
config = "EPTOFinal.avl"
files = getAVLFIle(CL_takeoff, config)
CDindArray = np.ones(len(CL_takeoff))
for i in range(len(CDindArray)):
    CDindArray[i] = readData(files[i])
    os.remove(files[i])
Cdo_TakeoffGearUp = 0.07277269152189407
Cdo_TakeoffGearDown = 0.07930746264865464
plt.plot(CDindArray + Cdo_TakeoffGearUp,CL_takeoff)
plt.plot(CDindArray + Cdo_TakeoffGearDown,CL_takeoff)

# Get data for landing config
CL_landing = np.linspace(-1.4,2.3,N) # Input landing config CL
config = "EPLandFinal.avl"
files = getAVLFIle(CL_landing, config)
CDindArray = np.ones(len(CL_landing))
for i in range(len(CDindArray)):
    CDindArray[i] = readData(files[i])
    os.remove(files[i])
Cdo_LandingGearUp = 0.10705693750239244
Cdo_LandingGearDown = 0.11359170862915301
plt.plot(CDindArray + Cdo_LandingGearUp,CL_landing)
plt.plot(CDindArray + Cdo_LandingGearDown,CL_landing)
plt.title("Drag Polars")
plt.xlabel(r'$C_{D}$')
plt.ylabel(r'$C_{L}$')
plt.legend(["Clean", "Takeoff Gear Up", "Takeoff Gear Down", "Landing Gear Up", "Landing Gear Down"])
plt.savefig('AVL Drag Polars.png',bbox_inches="tight")

plt.show()

