# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 04:15:23 2026

@author: rohdo
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 12:34:55 2025

@author: rohdo
"""

from equipartition import Equipartition
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, brentq

pb = 2.5

t = 0.2

tctilde = t * np.array([0, 1, 2])
colors = ["#553CA5", "#96B4DF", "#E69F00"]
linestyles = [":", "-.", "-"]
labels = [r"$\tilde{\theta}_c=0$", r"$\tilde{\theta}_c=\theta=0.2$", r"$\tilde{\theta}_c=2\theta=0.4$"]

tctildedivt = tctilde/t

uonlb = (1 + tctilde)**(-3 * (pb + 6)/(2 * pb + 13))
uoffub = (tctilde/(1 - np.cos(t)))**(3 * (pb + 6)/(2 * (pb + 10)))

u = np.linspace(0, 25, num = 3000)

#plt.title(r"Constraint $f$ for various $\tilde{\theta}_c/\theta$ $(\theta=" + str(t) + ",\overline{p}=" + str(pb) +")$")
plt.xlabel(r"Four velocity $u$")
plt.xlim(0, 25)
plt.ylabel(r"Constraint $f$")
plt.ylim(-0.4, 0.7)
plt.axhline(0, color = "k", linestyle = "--", alpha = 0.2)

f = lambda u : Equipartition.constraint_gammaBeta_corr(u, t, tctilde[2], pb)

ucmin = minimize(f, 1.0).x
uon = brentq(f, uonlb[2], ucmin)
uoff = brentq(f, uoffub[2], ucmin)

for i in range(len(tctildedivt)):
    if i == 2:
        plt.vlines(uonlb[i], ymin = -0.4, ymax = Equipartition.constraint_gammaBeta_corr(uonlb[i], t, tctilde[i], pb), color = colors[i], linestyle = "--")
        plt.annotate(r"$u_{on, lb}$", (uonlb[i] - 0.55, -0.35), rotation = 90)
        plt.vlines(uoffub[i], ymin = -0.4, ymax = Equipartition.constraint_gammaBeta_corr(uoffub[i], t, tctilde[i], pb), color = colors[i], linestyle = "--")
        plt.annotate(r"$u_{off, ub}$", (uoffub[i] - 0.55, -0.35), rotation = 90)
        plt.vlines(ucmin, ymin = -0.4, ymax = Equipartition.constraint_gammaBeta_corr(ucmin, t, tctilde[i], pb), color = colors[i], linestyle = "--")
        plt.annotate(r"$u_{min}$", (ucmin - 0.55, -0.35), rotation = 90)
        plt.plot(uon, 0, color = colors[i], marker = "o")
        plt.annotate(r"$u_{on}$", (uon, 0), xytext = (1, 3), textcoords = "offset points")
        plt.plot(uoff, 0, color = colors[i], marker = "s")
        plt.annotate(r"$u_{off}$", (uoff - 1.6, 0), xytext = (1, 3), textcoords = "offset points")
        
    y = Equipartition.constraint_gammaBeta_corr(u, t, tctilde[i], pb)
    plt.plot(u, y, color = colors[i], linestyle = linestyles[i], label = labels[i])
    
plt.legend(loc = "upper center")
plt.savefig("plotf.svg")
plt.savefig("plotf.png")