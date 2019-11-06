""""
Initializing code
Pablo Villanueva Domingo
Last update: 18/12/18
"""

import numpy as np
from Source.astro_params import *

# Init
index, C2, a2, C3, a3 = np.loadtxt("External_tables/stellar_spectra.dat",unpack=True)
nu_n = 4./3.*(1.-index**(-2))
nmax = len(index)
B2 = np.zeros(nmax)
Ntot2 = np.sum(Pop2*C2)
zeta2 = Ntot2*Msun/(mu*m_p)
#print Ntot2, Pop2*C2[0] # numbers from BL05
for n in range(nmax-1):
    B2[n] = Pop2*C2[n]*(a2[n]+1)/(nu_n[n+1]**(a2[n]+1)-nu_n[n]**(a2[n]+1))  # index as 21cmFAST
    #B2[n] = Pop2*C2[n]*(a2[n])/(nu_n[n+1]**(a2[n])-nu_n[n]**(a2[n]))    # index-1 as B&L05
#B2[-1]=C2[-1]   # last element is defined different
B2[-1]=Pop2*C2[-1]*(a2[-1]+1)/((4./3.*(1.-24**(-2)))**(a2[-1]+1)-nu_n[-1]**(a2[-1]+1))  # index as 21cmFAST
#B2[-1]=Pop2*C2[-1]*(a2[-1])/((4./3.*(1.-24**(-2)))**(a2[-1])-nu_n[-1]**(a2[-1]))    # index-1 as B&L05

ind_n, frec_n = np.loadtxt("External_tables/frecycle.dat",unpack=True)
