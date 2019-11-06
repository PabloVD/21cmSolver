""""
Runner for a grid of astrophysical parameters and several values of Delta
Pablo Villanueva Domingo
Last update: 6/11/19
"""
import numpy as np
import fileinput
import time, sys, os

time_ini = time.time()

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)

def routine(xi_ion,xi_heat,Tvir,Deltavec):

    print "{:.2e} {:.2e} {:.2e}".format(xi_ion, xi_heat, Tvir)
    replaceAll(paramfile,"xi_ion = ","xi_ion = {:2e} \n".format(xi_ion))
    replaceAll(paramfile,"xi_heat = ","xi_heat = {:2e} \n".format(xi_heat))
    replaceAll(paramfile,"Tvir = ","Tvir = {:2e} \n".format(Tvir))

    for i, D in enumerate(Deltavec):
        print "Delta = {:2e}".format(D)+", "+str(i+1)+" of "+str(len(Deltavec))
        replaceAll(solverfile,"Delta = ","Delta = {:2e} \n".format(D))
        os.system("python "+solverfile)

def Deltamin(zmax,zmin):
    return (zmax-zmin)/(1.+zmax)

paramfile = "Source/astro_params.py"
solverfile = "21cmEvolver.py"

# Astrophysical parameters grid
xi_ionvals = [1.]
xi_heatvals = [0.1,1.,10.]
Tvirvals = [1.e3,1.e4]

# For Delta constant:
Deltavec = np.logspace(np.log10(5.e-2),np.log10(7.),num=100)
# For Delta growing linearly:
#Deltavec = np.linspace(0.1,1.9,num=50)

#--- MAIN ---#

print "xi_ion   xi_heat     Tvir"
for i, xi_ion in enumerate(xi_ionvals):
    for n, xi_heat in enumerate(xi_heatvals):
        for j, Tvir in enumerate(Tvirvals):

            routine(xi_ion,xi_heat,Tvir,Deltavec)
            routine(xi_ion,xi_heat,Tvir,[1.])

print "Minutes elapsed:",(time.time()-time_ini)/60.
