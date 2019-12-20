""""
Runner for a grid of astrophysical parameters and several values of Delta
Pablo Villanueva Domingo
Last update: 6/11/19
"""
import numpy as np
import fileinput
import time, sys, os

time_ini = time.time()

# Function which replaces strings
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)

# Routine for solving the thermal history and compute the 21 cm signal
def Routine_Solver(xi_ion,xi_heat,Tvir,Deltavec):

    print("{:.2e} {:.2e} {:.2e}".format(xi_ion, xi_heat, Tvir))
    replaceAll(paramfile,"xi_ion = ","xi_ion = {:2e} \n".format(xi_ion))
    replaceAll(paramfile,"xi_heat = ","xi_heat = {:2e} \n".format(xi_heat))
    replaceAll(paramfile,"Tvir = ","Tvir = {:2e} \n".format(Tvir))

    for i, D in enumerate(Deltavec):
        print("Delta = {:2e}".format(D)+", "+str(i+1)+" of "+str(len(Deltavec)))
        replaceAll(solverfile,"Delta = ","Delta = {:2e} \n".format(D))
        os.system("python "+solverfile)

paramfile = "Source/astro_params.py"
solverfile = "21cmEvolver.py"

# Astrophysical parameters grid
xi_ionvals = [0.]
xi_heatvals = [0.]#[0.1,1.,10.]#np.logspace(-3,2,num=50)#[0.1,1.,10.]
#xi_heatvals = [0.1,1.,10.]
Tvirvals = [1.e4]

xi_ionvals = [0.]
xi_heatvals = [0.1,1.,10.]
Tvirvals = [1.e3,1.e4]

# numDelta = 30 is enough between 1.e-2 and 30
minDelta, maxDelta, numDelta = 1.e-2, 10., 30

# For Delta constant:
Deltavec = [1.]#np.logspace(np.log10(minDelta),np.log10(maxDelta),num=numDelta)
# For Delta growing linearly:
#Deltavec = np.linspace(0.1,1.9,num=50)

#--- MAIN ---#

# Run thermal history for different astro params and different Deltas
print("Solving grid \n")
print("xi_ion   xi_heat     Tvir")
for i, xi_ion in enumerate(xi_ionvals):
    for n, xi_heat in enumerate(xi_heatvals):
        for j, Tvir in enumerate(Tvirvals):

            Routine_Solver(xi_ion,xi_heat,Tvir,Deltavec)

print("Minutes elapsed:",(time.time()-time_ini)/60.)
