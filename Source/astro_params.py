"""
Astrophysical parameters
Pablo Villanueva Domingo
Last update: 6/11/19
"""

import numpy as np
from Source.constants import *

#--- ASTROPHYSICAL PARAMETERS ---#

# Fraction of baryons into stars
fstar = 0.01

# Ionization efficiency, Nion*fstar*fesc
# Nion: number of ionizing photons per baryons
# fesc: fraction of photons escaping to the IGM
xi_ion = 1.000000e+00 

# Minimum energy for X-rays, eV
hnu_thresh = 0.5e3

# Spectral index X-rays
a_spec = 1.2

# Heating efficiency, deviation from L/SFR in Furlanetto '06
xi_heat = 1.000000e+01 

# X ray luminosity over SFR, erg/s/(Msun/yr)
LSFR = 3.4e40*xi_heat

# Number of X ray photons per baryon
NX = LSFR/(hnu_thresh/kelvintoeV/(1./kelvintoerg*year_sec*mu*m_p/Msun))/(a_spec/(a_spec-1))

# Minimum virial temperature
Tvir = 1.000000e+04 

# Number of ionizing photons per baryon of Pop II stars
Pop2 = 4361

# Cutoff mass for integrals
cutoff_mass = 1e20

# Initial redshift
z_init = 35.

# Final redshift
z_end = 10.

#--- APPROXIMATIONS ---#

# Use analytical approximation for Jalpha_anal
Ja_anal = 1

# Use on-the-spot approximation for x-ray heating
use_ots = 1

# Use strong Ly alpha coupling approximation
strongcoupling = 0

# Use linear growth for the overdensities
lingrowth = 0
