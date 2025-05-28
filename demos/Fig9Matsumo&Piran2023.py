# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 16:58:13 2025

@author: rohdo
"""
import pandas as pd
import numpy as np
from equipartition import Equipartition
import matplotlib.pyplot as plt
import matplotlib

fileName = "Cendes_et_al_2021.xlsx" # Table 2 data AT 2019dsg

data = pd.read_excel(fileName).to_numpy()

z = data[:, 7]
tdays = data[:, 0]
FpmJy = data[:, 2]
logNup = data[:, 5]

z = z[2:] # drop the lowest two times
tdays = tdays[2:]
FpmJy = FpmJy[2:]
logNup = logNup[2:]

nup10 = 10**logNup/1e10

theta = np.logspace(-0.6, 0.2, 9)

equip_Newtonian = Equipartition(FpmJy, nup10, tdays, z, 0, newtonian = True)
equip_onAxis = [Equipartition(FpmJy, nup10, tdays, z, thet, newtonian = False, onAxis = True) for thet in theta]
equip_offAxis = [Equipartition(FpmJy, nup10, tdays, z, thet, newtonian = False, onAxis = False) for thet in theta]

cmap = plt.cm.jet  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#cmaplist[0] = (.5, .5, .5, 1.0)
# create the new map
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
bounds = np.arange(-0.6, 0.2 + 1e-14, 0.1)
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(4, 1, layout='constrained', figsize=(4,10))
ax1 = ax[0]
ax2 = ax[1]
ax3 = ax[2]
ax4 = ax[3]

xleft = 55
xright = 600

#plt.subplot(4, 1, 1)
#plt.title("AT 2019dsg")
ax1.set_xscale("log")
ax1.set_xlim(xleft, xright)
ax1.set_yscale("log")
ax1.set_ylim(1e-2, 1)
ax1.set_ylabel(r"$\beta_{eq, N}$")
ax1.axhline(0.23, color = "gray")
ax1.plot(tdays, equip_Newtonian.betaeqN(), color = "k", marker = "o")

#plt.subplot(4, 1, 2)
ax2.set_xscale("log")
ax2.set_xlim(xleft, xright)
ax2.set_yscale("log")
ax2.set_ylim(1e-2, 1e2)
ax2.set_ylabel(r"$\Gamma\beta$")
ax2.plot(tdays, equip_Newtonian.gammaBeta(), color = "k", marker = "o")
for i in range(0, len(theta) - 1):
    #ax2.plot(tdays, equip_onAxis[i].gammaBeta(), color = cmap(norm(np.log10(theta[i]))), alpha = 0.5)#; print(equip_onAxis[i].gammaBulk())
    #ax2.plot(tdays, equip_offAxis[i].gammaBeta(), color = cmap(norm(np.log10(theta[i]))), alpha = 0.5)#; print(equip_offAxis[i].gammaBeta())
    ax2.fill_between(tdays, equip_onAxis[i].gammaBeta(), equip_onAxis[i + 1].gammaBeta(), color = cmap(norm(np.log10(theta[i])/2 + np.log10(theta[i + 1])/2)))
    ax2.fill_between(tdays, equip_offAxis[i].gammaBeta(), equip_offAxis[i + 1].gammaBeta(), color = cmap(norm(np.log10(theta[i])/2 + np.log10(theta[i + 1])/2)))

#plt.subplot(4, 1, 3)
ax3.set_xscale("log")
ax3.set_xlim(xleft, xright)
ax3.set_yscale("log")
ax3.set_ylim(3e15, 3e18)
ax3.set_ylabel(r"$R_{eq}$ [cm]")
ax3.plot(tdays, equip_Newtonian.Req(), color = "k", marker = "o")
for i in range(0, len(theta) - 1):
    ax3.fill_between(tdays, equip_onAxis[i].Req(), equip_onAxis[i + 1].Req(), color = cmap(norm(np.log10(theta[i])/2 + np.log10(theta[i + 1])/2)))
    ax3.fill_between(tdays, equip_offAxis[i].Req(), equip_offAxis[i + 1].Req(), color = cmap(norm(np.log10(theta[i])/2 + np.log10(theta[i + 1])/2)))

#plt.subplot(4, 1, 4)
ax4.set_xscale("log")
ax4.set_xlim(xleft, xright)
ax4.set_yscale("log")
ax4.set_xlabel(r"Time since outflow launch : $t$ [day]")
ax4.set_ylim(1e46, 3e51)
ax4.set_ylabel(r"$E_{eq}$ [erg]")
ax4.plot(tdays, equip_Newtonian.energyeq(), color = "k", marker = "o")#, lw = 0.5)
for i in range(0, len(theta) - 1):
    ax4.fill_between(tdays, equip_onAxis[i].energyeq(), equip_onAxis[i + 1].energyeq(), color = cmap(norm(np.log10(theta[i])/2 + np.log10(theta[i + 1])/2)))
    ax4.fill_between(tdays, equip_offAxis[i].energyeq(), equip_offAxis[i + 1].energyeq(), color = cmap(norm(np.log10(theta[i])/2 + np.log10(theta[i + 1])/2)))

#plt.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
plt.tight_layout()
cax = fig.add_axes([0.4, 0.305, 0.5, 0.005])
cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap), cax = cax, orientation = "horizontal", label = r"$\log_{10}(\theta)$", location = "top")
cbar.ax.tick_params(labelsize=7)