# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 10:29:42 2025

@author: rohdo
"""
import numpy as np
from equipartition import Equipartition

t, nupGHz, _, nupGHz_err, FpmJy, _, FpmJy_err = np.loadtxt("apjad675bt3_ascii.txt", dtype = str, skiprows = 6).T

t = np.float64(t); print("day:", t)
nupGHz = np.float64(nupGHz); print("peak freq (GHz):", nupGHz)
nupGHz_err = np.float64(nupGHz_err); #print(nupGHz_err)
FpmJy = np.float64(FpmJy); print("peak flux (mJy):", FpmJy)
FpmJy_err = np.float64(FpmJy_err); #print(FpmJy_err)

tdays = t
nup10 = nupGHz/10
nup10_err = nupGHz_err/10

z = 0.0262
p = 2.8 # ?

equip_newtonian    = Equipartition(FpmJy, nup10, tdays, z,    0, p = p, epse = 0.1, epsB =   0.2, fA = 1, fV = 0.36, newtonian = True, isoNewtonianNe = True)
equip_offAxis_thin = Equipartition(FpmJy, nup10, tdays, z, 1.05, p = p, epse = 0.1, epsB = 0.003, fA = 1, fV =    1, onAxis = False)
equip_offAxis_wide = Equipartition(FpmJy, nup10, tdays, z, 1.57, p = p, epse = 0.1, epsB = 0.005, fA = 1, fV =    1, onAxis = False)

print("cosmology:", equip_newtonian.cosmo)

s =     "newtonian\nt      log(R) log(E) log(B) log(Ne) log(next) beta\n"

for i in range(len(t)):
    s +="{}   {:.2f}   {:.2f}   {:.2f}   {:.2f}   {:.2f}      {:.3f}\n".format(\
        equip_newtonian.tdays[i],\
        np.log10(equip_newtonian.Req())[i],\
        np.log10(equip_newtonian.energyeq())[i],\
        np.log10(equip_newtonian.magField())[i],\
        np.log10(equip_newtonian.Ne())[i],\
        np.log10(equip_newtonian.CNMnumDens())[i],\
        equip_newtonian.betaeqN()[i])
        
s +=    "off axis theta = 1.05\nt      log(R) log(E) log(B) log(Ne) log(next) Gamma\n"

for i in range(len(t)):
    s +="{}   {:.2f}  {:.2f}   {:.2f}    {:.2f}   {:.2f}      {:.1f}\n".format(\
        equip_offAxis_thin.tdays[i],\
        np.log10(equip_offAxis_thin.Req())[i],\
        np.log10(equip_offAxis_thin.energyeq())[i],\
        np.log10(equip_offAxis_thin.magField())[i],\
        np.log10(equip_offAxis_thin.Ne())[i],\
        np.log10(equip_offAxis_thin.CNMnumDens())[i],\
        equip_offAxis_thin.gammaBulk()[i])
        
s +=    "off axis theta = 1.57\nt      log(R) log(E) log(B) log(Ne) log(next) Gamma\n"

for i in range(len(t)):
    s +="{}   {:.2f}  {:.2f}   {:.2f}    {:.2f}   {:.2f}      {:.1f}\n".format(\
        equip_offAxis_wide.tdays[i],\
        np.log10(equip_offAxis_wide.Req())[i],\
        np.log10(equip_offAxis_wide.energyeq())[i],\
        np.log10(equip_offAxis_wide.magField())[i],\
        np.log10(equip_offAxis_wide.Ne())[i],\
        np.log10(equip_offAxis_wide.CNMnumDens())[i],\
        equip_offAxis_wide.gammaBulk()[i])

print(s)

with open("ChristyTable4.txt", "w") as text_file:
    text_file.write(s)