""""
21cm evolver
Pablo Villanueva Domingo
Last update: 6/11/19
"""

# Notes: warning error with Tvir=1.e4 !
# error with high heating in caseA
# Test ionization, T, lyalpha, dTb with Furlanetto
# Check eq. x_e, correction LLS

from scipy.integrate import solve_ivp
from Source.functions import *
import os, time

#--- INITIAL ARRANGEMENTS ---#

T_init =  Tk_ad(z_init)
xe_init = 1.e-4
zvec = np.linspace(z_init,z_end,num=50)
lnavec = np.log(1./(1.+zvec))
Delta = 1.000000e+00 
Ttilde_init = TtoTtilde(T_init,xe_init,z_init,1.)

#--- MAIN ---#

# Solving evolution equations
sol = solve_ivp(lambda lna, y: EvolutionEquations(lna, y, Delta), [lnavec[0],lnavec[-1]], [0., xe_init, Ttilde_init], t_eval = lnavec)
Q, xe, Tt = sol.y[0], sol.y[1], sol.y[2]
Tk = TtildetoT(Tt,xe,zvec,Delta)
xHIIbar = Q + (1.-Q)*xe

# Effective cooling at T = 1.e4 K
for tt, T in enumerate(Tk):
    if T>1.e4:  Tk[tt] = 1.e4
# Correct possible overionization
for tt, x in enumerate(xHIIbar):
    if x>1.:  xHIIbar[tt] = 1.

Ja = Jalpha(Tvir,zvec)
Tspin = Ts(zvec,y_coll(zvec,1.-xHIIbar,Tk,Delta),y_lya(zvec,Ja,Tk),Tk)

if not os.path.exists("Outputs"):   os.system("mkdir Outputs")
nameoutput = "Outputs/Evolution_xi_ion_{:.2e}_xi_heat_{:.2e}_tvir_{:.2e}_Delta_{:.2e}.dat".format(xi_ion,xi_heat,Tvir,Delta)
np.savetxt(nameoutput, np.transpose([zvec,xHIIbar,xe,Tk,Tspin,T21(zvec,1.-xHIIbar,Tspin,Delta),Ja]), fmt="%1.3e", header="\t\t z \t\t\t Q \t\t xe \t\t Tk [K] \t\t Ts [K] \t\t dTb [mK] \t J_alpha [cm^-2 s^-1 Hz^-1] ")
