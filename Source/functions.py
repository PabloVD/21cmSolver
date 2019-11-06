""""
Functions for the 21 cm code
Pablo Villanueva Domingo
Last update: 6/11/19
"""

import numpy as np
import math
from scipy import integrate
from scipy.integrate import solve_ivp
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from Source.constants import *
from Source.astro_params import *
from Source.initialize import *
from Source.lyalpha_rt import *
from scipy import special

cosmo = cosmology.setCosmology('planck18')
rhoc = cosmo.rho_c(0)*1e9   # in M_sun Mpc^-3 h^2

#--- Miscelanea ---#

def scinot(x):      # write in Latex scientific notation
    exp = int(math.floor(math.log10(abs(x))))
    prefactor = x / 10**exp
    if prefactor == 1.:
        return "10^{"+str(exp)+"}"
    else:
        return "{:.0f}".format(prefactor)+" \\times 10^{"+str(exp)+"}"

def Hz(z):
    return H0*np.sqrt(Omega_m*(1.+z)**3. + Omega_l)

def dtdz(z):
    return -1./Hz(z)/(1.+z)

def n_H(z):     # cm^-3
    return n0_h*(1.+z)**3.

def Tk_ad(z):
    return Tcmb0*(1.+zdec)*((1.+z)/(1.+zdec))**2.

def Mmin(Tvir,z):       # Msun/h, tvir in K
    oz = Omega_m*(1.+z)**3./(Omega_m*(1.+z)**3. + 1.-Omega_m)
    d = oz -1.
    Deltac = 18.*np.pi**2. +82.*d - 39.*d**2.
    return 1.e8*(Omega_m*Deltac/(oz*18.*np.pi**2.))**(-1./2.)*(Tvir/1.98e4)**(3./2.)*((1.+z)/10.)**(-3./2.)

def dndlnM(M,z):         # dndlnM units of (h/Mpc)^3
    return mass_function.massFunction(M,z, q_out = "dndlnM", model = "sheth99")

def Fcoll(Tvir,z):
    logM = np.linspace(np.log(Mmin(Tvir,z)),np.log(cutoff_mass))
    integral = integrate.simps(np.exp(logM)*dndlnM(np.exp(logM),z), logM )
    return integral/Omega_m/rhoc

def dFcolldz(Tvir,z):
    dz = 0.01
    return (Fcoll(Tvir,z+dz)-Fcoll(Tvir,z-dz))/2./dz

def sigma(nu,a,nuT,beta,s):
    if nu<nuT:
        return 0.
    else:
        return a*( beta*(nu/nuT)**(-s) + (1.-beta)*(nu/nuT)**(-s-1.))

def Intsigma(nuth,s,numin,nu0,alpha):
    return (nuth/nu0)**(s)*(nu0/numin)**(alpha+s-1.)*alpha/(alpha + s -1.)*(1. - (nuth/numin)*(alpha + s -1.)/(alpha + s) )

def sigmaeff(a,b,nuth,s,numin,nu0,alpha):
    return a*( b*Intsigma(nuth,s,numin,nu0,alpha) + (1.-b)*Intsigma(nuth,s+1.,numin,nu0,alpha) )

def sigmaefftot(xe,numin,nu0,alpha):
    sigHI = f_H*(1.-xe)*sigmaeff(alist[0],blist[0],nuTlist[0],slist[0],numin,nu0,alpha)
    sigHeI = f_He*(1.-xe)*sigmaeff(alist[1],blist[1],nuTlist[1],slist[1],numin,nu0,alpha)
    sigHeII = f_He*xe*sigmaeff(alist[2],blist[2],nuTlist[2],slist[2],numin,nu0,alpha)
    return sigHI + sigHeI + sigHeII

def hnu_tauone(z,zprime):    # nu when tau=1, eV, see Mirocha 2014, eq 17
    hnu_mu = 366.5           # eV
    return hnu_mu*(1.+z)**(1./2.)*( 1. - ((1.+z)/(1.+zprime))**(3./2.) )**(1./3.)

def hnumin(z,zprime):
    return  max( hnu_thresh, hnu_tauone(z,zprime) )

hnumin = np.vectorize(hnumin)

def frac_heat(x):   # From Furlanetto 2006, Shull & Steenberg 1985, check x_e or x_HII
    C, a, b = 0.9971, 0.2663, 1.3163
    if x>=1.:
        return C
    elif x > 1.e-4:
        return C*(1.-(1.-x**(a))**(b))
    else:
        return 0.15

#--- Ionization evolution ---#

def alpha_recA(T):   # Case A, from Abel '97, cm^3/s
    logT = np.log(T*kelvintoeV)
    if T<1.: logT = np.log(1.*kelvintoeV)   # fit works for T from 1 to 1e8 K
    caseA = np.exp( -28.6130338 - 0.72411256*logT - 2.02604473e-2*logT**2. - 2.38086188e-3*logT**3. - 3.21260521e-4*logT**4.- 1.42150291e-5*logT**5. + 4.98910892e-6*logT**6. + 5.75561414e-7*logT**7. - 1.85676704e-8*logT**8. - 3.07113524e-9 * logT**9.)
    return caseA

def alpha_recB(T):   # Case B, from Spitzer '78, 21cmFAST, cm^3/s
    #caseB = 2.6e-13*(T/1.e4)**(-0.85)   # Case B, from Bernal et al.
    caseB = 2.59e-13*(T/1.e4)**(-0.75)
    return caseB

def opacity_LLS(z): # Madau 2017, Mpc^-1
    if z<5.5:
        return (37.*((1.+z)/5.)**(-5.4))**(1.)
    else:
        return 0.

def opacity_IGM(z,Q): # Madau 2017, Mpc^-1
    return n_H(z)*sigmaH*(1.-Q)*MpcToCm

def LLS_fac(z,Q): # Madau 2017
    return 1. #1./(1.+opacity_LLS(z)/opacity_IGM(z,Q))

def t_rec(z):   # Recombination time, s
    Clump = 1.
    return (alpha_recA(1.e4)*(1.+Y_He/(4*(1.-Y_He)))*Clump*n_H(z))**(-1.)

#--- Heating rates per baryon density over H(z), \epsilon/n_b/H [K] ---#

def xray_heat_full( z, xe ): # UPDATE
    zprime = np.linspace( z, z_init+1. , num=100)
    int = integrate.simps( sigmaefftot(xe,hnumin(z,zprime),hnu_thresh,a_spec)*(-1.)*dFcolldz(Tvir,zprime)*((1.+zprime)/(1.+z))**(-a_spec), zprime )
    return hnu_thresh/kelvintoeV*xi_heat/(a_spec/(a_spec-1.))*frac_heat(xe)*c*(1.+z)**3.*n0_b*int/Hz(z)

def xray_heat_ots( z, xe ):    # On the spot aprox. xi_heat is f_X of Furlanetto '06
    # old xi_heat = 0.012 and hnu0=0.5 keV correspond to fiducial values of eq. 11 of Furlanetto '06, f_X=1
    # return f_X*frac_heat(xe)*fstar*mu*m_p/Msun*year_sec/kelvintoerg*3.4e40*dFcolldz(Tvir,z)/dtdz(z)/Hz(z)
    #return hnu_thresh/kelvintoeV*xi_heat*frac_heat(xe)*dFcolldz(Tvir,z)/dtdz(z)/Hz(z)
    #LSFR = hnu_thresh/kelvintoeV*xi_heat/(1./kelvintoerg*year_sec*mu*m_p/Msun*fstar)
    return LSFR/kelvintoerg*year_sec*mu*m_p/Msun*fstar*frac_heat(xe)*dFcolldz(Tvir,z)/dtdz(z)/Hz(z)

def xray_heat_ots_complete( z, xe ):    # On the spot aprox. with 2nd term with \nu_th
    return xray_heat_ots(z, xe) - (nuTlist[0]+nuTlist[1]+nuTlist[2])/kelvintoeV*xi_heat/(a_spec/(a_spec-1.))*frac_heat(xe)*dFcolldz(Tvir,z)/dtdz(z)/Hz(z)

def xray_heat( z, xe ):
    if use_ots:
        return xray_heat_ots( z, xe )
    else:
        return xray_heat_full( z, xe )

def heat_21( z, xH, Tk, Delta ):
    #if Tk<1.e-4:    Tk=1.e-4
    if strongcoupling:    Tspin=Tk
    else:               Tspin = Ts(z,y_coll(z,xH,Tk,Delta),y_lya(z,Jalpha(Tvir,z),Tk),Tk)
    #Tspin=Tk
    #if Tspin<Tk:    Tspin=Tk
    heat21 = 3./4.*Tstar*A10*xH*n0_h/n0_b/Hz(z)*(Tcmb0*(1.+z)/Tspin - 1.)
    #print z, Tspin/(Tcmb0*(1.+z)), Tk/(Tcmb0*(1.+z)), heat21
    return heat21

def compton_cool( z, xe, Tk ):
    return 4.*sigmaT*c*rho_cmb0*(1.+z)**4./m_ec2*xe*(Tcmb0*(1.+z)-Tk)/Hz(z)

def heat_lyalpha( z, Tk, Delta ):      # Furlanetto, Pritchard 2006
    #return 4.*np.pi/c*hplanck*nualpha/kelvintoerg*Ic_lyalpha(z,Tk)*DeltanuDalpha(Tk)*Jalpha(Tvir,z)/(n0_b*(1.+z)**3.)
    Jcont = Jalpha_cont( Tvir, z )
    return 4.*np.pi/c*hplanck*nualpha/kelvintoerg*DeltanuDalpha(Tk)/(n0_b*(1.+z)**3.*Delta)*(Ic_lyalpha(z,Tk)*Jcont + Ii_lyalpha(z,Tk)*Jalpha_inject(Tvir,z,Jcont))

def Heat_rate( z, xe, Tk, Delta ):       # Total heating rate
    return  xray_heat(z,xe) + compton_cool(z,xe,Tk) + heat_21(z,1.-xe,Tk,Delta) + heat_lyalpha(z, Tk, Delta)

def TtoTtilde( T, xe, z, Delta ):   # Define Ttilde = constant if no heating
    return T*(1.+xe)/(1.+z)**2./Delta**(2./3.)

def TtildetoT( Ttilde, xe, z, Delta ):  # Transform Ttilde to T
    return Ttilde*Delta**(2./3.)*(1.+z)**2./(1.+xe)

#--- Evolution routine ---#

def EvolutionEquations(lna, y, Delta): # = dy/dlna (=dy/dt/H), y = [Q, x_e, Ttilde], Ttilde=Tk*(1+x_e)*a^2/Delta^(2/3)

    z = np.exp(-lna)-1.

    # use linear growth for delta, otherwise, remain constant
    if lingrowth:   Deltaz = 1. + cosmo.growthFactor(z)/cosmo.growthFactor(z_init)*(Delta-1.)
    else:           Deltaz = Delta
    if Deltaz<0.:   # to check errors
        print z, Deltaz

    Q, xe, Tt = y[0], y[1], y[2]
    if Tt<0.:    Tt = 1.e-4
    T = TtildetoT(Tt,xe,z,Deltaz)

    if T>1.e4:  T = 1.e4        # effective cooling
    if Q>1.: Q=1.               # maximum ionization

    if xi_ion==0.: ion_eq=0.    # If no UV-rays, switch off reionization to avoid numerical problems
    else:
        ion_eq = xi_ion*LLS_fac(z,Q)*(-(1.+z))*dFcolldz(Tvir,z) - Q/t_rec(z)/Hz(z)   # dQ/dlna

    if xi_heat==0.: xe_eq=0.    # If no X-rays, switch off ionizations in neutral IGM to avoid numerical problems
    else:
        xe_eq = NX*fstar*(-(1.+z))*dFcolldz(Tvir,z)*(1.-xe) - alpha_recB(T)*n_H(z)/Hz(z)*xe**2.   # dx_e/dlna, ionized fraction in neutral regions, ots aprox assumed

    #temp_eq = -2*T - T*xe_eq/(1.+xe) + 2./3.*Heat_rate(xi_heat, Tvir, z, xe, T*Deltaz**(2./3.), Deltaz)/(1.+xe)/Deltaz**(2./3.)        # dT/Deltaz^(2/3)/dlna
    temp_eq = 2./3.*Heat_rate( z, xe, T, Deltaz)/Deltaz**(2./3.)/(1.+z)**2.        # dTtilde/dlna

    return [ion_eq, xe_eq, temp_eq]

#--- Lyalfa flux ---#

# Emissivity Lyalfa (Hz^-1) from BL05 (same in 21cmFAST)
def EmiLy(nu_norm): # nu_norm=nu/nualpha
    for n in range(nmax-1):     # Using index as 21cmFAST
        if (nu_norm >= nu_n[n]) and (nu_norm < nu_n[n+1]):
            return B2[n]*nu_norm**a2[n]/nualpha
    return B2[-1]*nu_norm**a2[-1]/nualpha
    """for n in range(nmax-1): # Using index-1 as Barkana and Loeb '05
        if (nu_norm >= nu_n[n]) and (nu_norm < nu_n[n+1]):
            return B2[n]*nu_norm**(a2[n]-1.)/nualpha
    return B2[-1]*nu_norm**(a2[-1]-1.)/nualpha"""

# Ly alpha background, cm^-2 s^-1 Hz^-1
def Jalpha_full(Tvir,z):     # Full computation using BL05 emissivity
    sum = 0.
    for inn, n in enumerate(ind_n):
        zmaxn = (1.+z)*(1.-(n+1)**(-2.))/(1.-(n)**(-2.)) - 1.
        zz = np.linspace(z,zmaxn,num=15)
        integral_n = integrate.simps( -(1.+zz)*dFcolldz(Tvir,zz)*EmiLy(nu_n[inn]*(1.+zz)/(1.+z)) , zz )
        sum+= frec_n[inn]*integral_n
    Ja = c*(1.+z)**2./(4.*np.pi)*n0_b*fstar*sum
    return Ja

def Jalpha_anal(Tvir,z):     # Analytical aprox. using power law nu^-1
    sum = 0.
    for inn, n in enumerate(index):
        zmaxn = (1.+z)*(1.-(n+1.)**(-2.))/(1.-(n)**(-2.)) - 1.
        sum+= frec_n[inn]/nu_n[inn]*( Fcoll(Tvir,z) - Fcoll(Tvir,zmaxn) )
    Ja = c*(1.+z)**3./(4.*np.pi)*n0_b*fstar*Ntot2/np.log(4./3.)/nualpha*sum
    return Ja

def Jalpha_aprox(Tvir,z):     #  Constant emissivity constant, aprox
    sum = 0.
    for inn, n in enumerate(index):
        zmaxn = (1.+z)*(1.-(n+1)**(-2.))/(1.-(n)**(-2.)) - 1.
        sum += frec_n[inn]*( Fcoll(Tvir,z) - (1.+zmaxn)/(1.+z)*Fcoll(Tvir,zmaxn) )
    Ja = c*(1.+z)**3./(4.*np.pi)*n0_b*fstar*Ntot2/(nualpha/3.)*sum
    return Ja

def Jalpha_F06(Tvir,z):     # Furlanetto '06 aprox
    frec_ave = 0.72    # for PopII, 0.63 for PopIII
    Ja = Ntot2*fstar*frec_ave*c*(1.+z)**2./(4.*np.pi*nualpha/3.)*n0_b*Fcoll(Tvir,z)  # Seems to work better with (1+z)^3
    return Ja

def Jalpha_M15(Tvir,z):     # Mirocha '15 aprox
    Ja = Ntot2*fstar*c*(1.+z)**2./(4.*np.pi*nualpha/3.*Hz(z))*n0_b*dFcolldz(Tvir,z)/dtdz(z)
    return Ja

def Jalpha_cont(Tvir,z):    # J_alpha only from continuous bakcground
    inn, n = 0, 2
    zmaxn = (1.+z)*(1.-(n+1.)**(-2.))/(1.-(n)**(-2.)) - 1.
    if Ja_anal:
        sum = frec_n[inn]/nu_n[inn]*( Fcoll(Tvir,z) - Fcoll(Tvir,zmaxn) )
        Ja = c*(1.+z)**3./(4.*np.pi)*n0_b*fstar*Ntot2/np.log(4./3.)/nualpha*sum
    else:
        zz = np.linspace(z,zmaxn,num=15)
        integral_n = integrate.simps( -(1.+zz)*dFcolldz(Tvir,zz)*EmiLy(nu_n[inn]*(1.+zz)/(1.+z)) , zz )
        sum = frec_n[inn]*integral_n
        Ja = c*(1.+z)**2./(4.*np.pi)*n0_b*fstar*sum
    return Ja

def Jalpha_inject(Tvir,z,Jcont):  # J_alpha only from injected source
    return Jalpha(Tvir,z)-Jcont

def Jalpha(Tvir,z):
    if Ja_anal: return Jalpha_anal(Tvir,z)
    else:       return Jalpha_full(Tvir,z)

def J0(z):      # One Lyalpha photon per H atom, cm^-2 s^-1 Hz^-1   NOT USED
    return c*n_H(z)/(4.*np.pi*nualpha)

def x_alpha(Jalph,z):      # Jalph in units of cm^-2 s^-1 Hz^-1 NOT USED
    return Jalph/J0(z)*((1.+z)/21.)**2./0.069

#--- 21cm Signal ---#

def y_lya(z,Jalph,T):   # Ly-alpha coupling coefficient, y_alpha = P_10/A_10*Tstar/T_k
    if strongcoupling:
        return 1.e3
    else:
        return  (16.*np.pi**2.*Tstar*e2*f12*c)/(27.*A10*T*m_ec2*eVtoerg)*S_alpha_wing(z,T)*Jalph

def kappa_eH(T):
    if T<1.:
        ke = 0.
    elif T>=1. and T<= 1e4:
        ke = 10.**( -9.607 + 0.5*np.log10(T)*np.exp(-(np.log10(T))**4.5/1800.) )
    else:
        ke = 10.**( -9.607 + 0.5*np.log10(1.e4)*np.exp(-(np.log10(1.e4))**4.5/1800.) )
    return ke

def kappa_HH(T):
    if T>10.:
        kh = 3.1e-11*T**(0.357)*np.exp(-32./T)
    elif T<=10 and T>1:
        kh = 3.6e-16*T**(3.640)*np.exp(6.035/T)
    else:
        kh = 3.6e-16*1.**(3.640)*np.exp(6.035/1.)
    return kh

def kappa_pH(T):    # My own parameterization, based on Furnaletto & Furnaletto '07
    if T>130.:
        kp = 2.*kappa_HH(T)
    else:
        kp = 2.*kappa_HH(130.)
    return kp

def y_coll(z,xH,T,Delta):   # Collisional coupling coefficient
    return Tstar/T/A10*n_H(z)*Delta*( kappa_HH(T)*xH + kappa_eH(T)*(1.-xH) + kappa_pH(T)*(1.-xH) )

def Ts(z,yc,yly,Tk):    # Spin temperature [K]
    return (Tcmb0*(1.+z) + (yc + yly)*Tk)/(1. + yc + yly)

def Ts_v2(z,yc,yly,Tk,Tc):  # Includes color temperature, not used
    xc, xly = yc*Tk/(Tcmb0*(1.+z)), yly*Tk/(Tcmb0*(1.+z))
    Tsinv = ( (Tcmb0*(1.+z))**(-1.) + xc*Tk**(-1.) + xly*Tc**(-1.) )/( 1. + xc +xly )
    return Tsinv**(-1.)

def T21(z,xH,Ts,Delta): # Global 21 cm signal [mK]
    return 27.*xH*Delta*(1.-(Tcmb0*(1.+z)/Ts))*np.sqrt((1.+z)/10.)

def Ts_iter():  # Iterative TS, to be finished
    Sa = S_alpha_Hirata(z,Tk,TS)
    Tc = T_col_Hirata(Tk,TS)

# Vectorize some functions in order to accept numpy arrays as input
Fcoll = np.vectorize(Fcoll)
dFcolldz = np.vectorize(dFcolldz)
Ts = np.vectorize(Ts)
kappa_eH = np.vectorize(kappa_eH)
kappa_HH = np.vectorize(kappa_HH)
kappa_pH = np.vectorize(kappa_pH)
EmiLy = np.vectorize(EmiLy)
Jalpha_full = np.vectorize(Jalpha_full)
Jalpha_cont = np.vectorize(Jalpha_cont)
