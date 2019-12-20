""""
Plotter code
Pablo Villanueva Domingo
Last update: 6/11/19
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Source.functions import *
import os

colors = ["red","blue","green","purple"]
lines = ["-","--",":","-."]

fig_x, (ax_x) = plt.subplots(1,1)
fig_T, (ax_T, ax_21) = plt.subplots(2, 1, sharex=True)
fig_Ja, (axJa) = plt.subplots(1, 1)
fig_T.subplots_adjust(hspace=0)

# Astrophysical parameters grid
xi_ionvals = [0.]
xi_heatvals = [0.1,1.,10.]
Tvirvals = [1.e3,1.e4]

Delta = 1.000000e+00

#--- MAIN ---#

for i, xi_ion in enumerate(xi_ionvals):
    for n, xi_heat in enumerate(xi_heatvals):
        for j, Tvir in enumerate(Tvirvals):

            nameoutput = "Outputs/Evolution_xi_ion_{:.2e}_xi_heat_{:.2e}_tvir_{:.2e}_Delta_{:.2e}.dat".format(xi_ion,xi_heat,Tvir,Delta)
            tab = np.loadtxt(nameoutput, unpack=True)
            zvec, xHIIbar, xe, Tk, Tspin, dTb, Ja = tab

            ax_x.semilogy(zvec, xHIIbar , color=colors[n], linestyle=lines[j])
            ax_x.semilogy(zvec, xe , color=colors[n], linestyle=lines[j],alpha=0.7, linewidth=0.5)
            ax_T.semilogy(zvec, Tk, color=colors[n], linestyle=lines[j])
            ax_T.semilogy(zvec, Tspin, color=colors[n], linestyle=lines[j], alpha=0.7, linewidth=0.5)
            ax_21.plot(zvec, dTb, color=colors[n], linestyle=lines[j])
            axJa.semilogy(zvec, Ja, color="k", linestyle=lines[j])

#--- PLOTS ---#

customlegend, Jalegend = [], []
for n, xi_heat in enumerate(xi_heatvals):
    customlegend.append( Line2D([0], [0], color=colors[n], lw=4, label="$\\xi_{heat}=$"+"{:.3f}".format(xi_heat)))
for j, TT in enumerate(Tvirvals):
    customlegend.append( Line2D([0], [0], color="black", linestyle=lines[j], label="$T_{vir}^{min}="+scinot(TT)+"$ K"))
    Jalegend.append( Line2D([0], [0], color="black", linestyle=lines[j], label="$T_{vir}^{min}="+scinot(TT)+"$ K"))

#ax_x.semilogy([6.,6.,],[1.e-3,1.],"k--") # End EoR z=6
ax_x.set_xlim(z_end,20.)
#ax_x.set_ylim(1.e-3,1.1)
ax_x.legend(handles=customlegend, loc='upper right',fontsize = 8,framealpha=1.)
ax_x.set_xlabel(r"$z$")
ax_x.set_ylabel(r"$Q_{HII}$")

ax_T.semilogy(zvec, Tk_ad(zvec), color="black", linestyle="-.", alpha=0.7, linewidth=0.5, label="Adiabatic")
ax_T.semilogy(zvec, Tcmb0*(1.+zvec), color="grey", linestyle="--", alpha=0.7, linewidth=0.5, label="CMB")
#ax_T.set_ylim(1.e0,1.e3)
ax_T.set_ylabel(r"$T$ [K]")

ax_21.plot([zvec[0],zvec[-1]],[0.,0.],":",alpha=0.2)
ax_21.set_xlim(z_end,z_init)
#ax_21.set_ylim(-300,10.)
ax_21.set_xlabel(r"$z$")
ax_21.set_ylabel(r"$\delta T_b \; [mK]$")
ax_21.legend(handles=customlegend, loc='lower right',fontsize = 8,framealpha=1.)

#axJa.set_ylim(0.,12.)
axJa.legend(handles=Jalegend, loc='upper right',fontsize = 8,framealpha=1.)
axJa.set_xlabel(r"$z$")
axJa.set_ylabel(r"$J_\alpha \; [cm^{-2} \cdot s^{-1} \cdot Hz^{-1}]$")

if not os.path.exists("Plots"):   os.system("mkdir Plots")
fig_x.savefig("Plots/Ionization_Evolution.pdf", bbox_inches='tight')
fig_T.savefig("Plots/Temperature_Evolution.pdf", bbox_inches='tight')
fig_Ja.savefig("Plots/J_alpha_Evolution.pdf", bbox_inches='tight')
#plt.show()
plt.gcf().clear()
