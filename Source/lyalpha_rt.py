""""
Ly alpha radiative transfer module
Functions below from Furlanetto & Pritchard 2006
Pablo Villanueva Domingo
Last update: 17/9/19
"""

import numpy as np
from scipy import integrate
from Source.constants import *
from scipy import special

def Hz(z):  # put in another module
    return H0*np.sqrt(Omega_m*(1.+z)**3. + Omega_l)

def DeltanuDalpha(Tk):  # Hz
    return nualpha*np.sqrt(2.*Tk*kelvintoeV/m_pc2)

def a_voigt(Tk):
    return Aalpha/(4.*np.pi*DeltanuDalpha(Tk))

def taugp(z):       # from 21cmFAST
    return 1.342881e-7*n0_h*(1.+z)**3./Hz(z)

def gammaGP(z):
    return 1./taugp(z)

def eta_lya(Tk):
    return hplanck*nualpha/m_pc2/eVtoerg*(nualpha/DeltanuDalpha(Tk))

def voigt(x,Tk):
    # taken form here https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/
    #return a_voigt(Tk)/np.pi**(3./2.)*...
    #sigma = DeltanuDalpha(Tk)/np.sqrt(2.)
    #return np.real(special.wofz((x + 1j*a_voigt(Tk))/sigma/np.sqrt(2.))) / sigma/np.sqrt(2.*np.pi)

    # here, x=(nu-nu0)/DeltanuDalpha(Tk), not like in the link
    z = x + 1j*a_voigt(Tk)
    return np.real(special.wofz(z)) /np.sqrt(np.pi)

def sigmavoigt(x,Tk):
    y = np.linspace(0,x)
    return integrate.simps(1./voigt(y,Tk), y)

xlim = 50.

def J1(x,z,Tk): # over Jinf; injected with x>0
    return Jat0(z,Tk)*np.exp( -2.*eta_lya(Tk)*x -2.*gammaGP(z)*sigmavoigt(x,Tk) )

def J2(x,z,Tk): # over Jinf; continuous back. and injected with x<0
    logyy = np.linspace(-5.,np.log(xlim))
    yy = np.exp(logyy)
    deltaJ = 2.*eta_lya(Tk)*integrate.simps( np.exp( -2.*eta_lya(Tk)*yy -2.*gammaGP(z)*( sigmavoigt(x,Tk) - sigmavoigt(x-yy,Tk) ) )*yy , logyy )
    return 1.-deltaJ

def Jinj(x,z,Tk):
    if x>0: return J1(x,z,Tk)
    else:   return J2(x,z,Tk)

def Jat0(z,Tk): # over Jinf
    logyy = np.linspace(-5.,np.log(xlim))
    yy = np.exp(logyy)
    deltaJat0 = 2.*eta_lya(Tk)*integrate.simps( np.exp( -2.*eta_lya(Tk)*yy -2.*gammaGP(z)*( sigmavoigt(yy,Tk) ) )*yy , logyy )
    return 1.-deltaJat0

def Ic_lyalpha(z,Tk):
    beta = eta_lya(Tk)*(4.*a_voigt(Tk)/np.pi/gammaGP(z))**(1./3.)
    return (4./np.pi)**(-1./6.)*np.pi**(3./2.)*(a_voigt(Tk)/gammaGP(z))**(1./3.)*beta*( special.airy(-beta)[0]**2. + special.airy(-beta)[2]**2. )

def Ii_lyalpha(z,Tk):
    beta = eta_lya(Tk)*(4.*a_voigt(Tk)/np.pi/gammaGP(z))**(1./3.)
    A0, A1, A2 = -0.6979, 2.5424, -2.5645
    #y = np.linspace(1.e-5,1.e3)
    #integ = integrate.simps(1./np.sqrt(y)*np.exp( -np.pi*gammaGP/6./a_voigt(Tk)*y**3. - 2.*eta_lya(Tk)*y )*special.erfc(np.sqrt(np.pi*gammaGP/2./a_voigt(Tk)*y**3.)), y )
    #return eta_lya(Tk)*np.sqrt(a_voigt(Tk)/2./gammaGP)*integ -
    return (a_voigt(Tk)/gammaGP(z))**(1./3.)*( A0 + A1*beta + A2*beta**2. )

def S_alpha_exact(z,Tk):    # using continuous
    x = np.linspace(-xlim,xlim,num=100)
    return integrate.simps( voigt(x,Tk)*J2(x,z,Tk), x)

def S_alpha_wing(z,Tk):
    alf = eta_lya(Tk)*(3.*a_voigt(Tk)/2./np.pi/gammaGP(z))**(1./3.)
    OneMinusSalpha = 4.*alf/9.*( 3.**(2./3.)*np.pi*special.airy(-2.*alf/3.**(1./3.))[2] + (3.*alf**2.)*special.hyp1f2(1.,4./3.,5./3.,-8.*alf**3./27.)[0] )
    return 1.-OneMinusSalpha

def S_alpha_wing_1st(z,Tk):
    alf = eta_lya(Tk)*(3.*a_voigt(Tk)/2./np.pi/gammaGP(z))**(1./3.)
    return 1.-4.*np.pi/(3.*np.sqrt(3.)*special.gamma(2./3.))*alf

# Hirata fits

def S_alpha_Hirata(z,Tk,TS):    # Tk, TS in K
    xi = (1.e-7*taugp(z))**(1./3.)*Tk**(-2./3.)
    return (1. -0.0631789*Tk**(-1.) +0.115995*Tk**(-2.) -0.401403*Tk**(-1.)*TS**(-1.) +0.336463*Tk**(-2.)*TS**(-1.) )*(1.+2.98394*xi + 1.53583*xi**2. + 3.85289*xi**3.)**(-1.)

def T_col_Hirata(Tk,TS):        # Tk, TS in K
    Tcinv = Tk**(-1.) + 0.405535*Tk**(-1.)*( TS**(-1.) - Tk**(-1.) )
    return Tcinv**(-1.)



sigmavoigt = np.vectorize(sigmavoigt)
J2 = np.vectorize(J2)
Jinj = np.vectorize(Jinj)
S_alpha_wing = np.vectorize(S_alpha_wing)
S_alpha_exact = np.vectorize(S_alpha_exact)
