# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 10:29:42 2025

@author: rohdo
"""
import numpy as np
from equipartition import Equipartition

def printlatextable(headers, data):
    """"""
    #headers = ["x 1","x 2","x 3","x 4","x 5","bbar"]
    #data["x 4"] = [1,2,3,1,0,9]
    #data["x 5"] = [3,2,2,0,1,15]
    #data["z"] = [1,9,3,0,0,0]
    
    #for k, v in kwargs.items():
    #    data[k] = v
    
    textabular = f"l|{'c'*len(headers)}"
    texheader = " & " + " & ".join(headers) + "\\\\"
    texdata = "\\hline\n"
    for label in data:
       texdata += f"{label} & {' & '.join(map(str,data[label]))} \\\\\n"
    
    print("\\begin{tabular}{"+textabular+"}")
    print(texheader)
    print(texdata,end="")
    print("\\end{tabular}")

def printchristytable(equip_newtonian, equip_offAxis_thin, equip_offAxis_wide):
    """"""
    headers = ["$t$", "$\\log(R)$ (cm)", "$\\log(E)$ (ergs)", "$\\log(B)$ (G)", "$\\log(N_e)$", "$\\log(n_{\\text{ext}})$ ($\\text{cm}^{-3}$)", "$\\beta$"]
    
    N_T = len(equip_newtonian)
    
    data = dict()
    
    keys = ["Spherical", "$f_A = {:.2f}$".format(equip_newtonian[0].fA), "$f_V = {:.2f}$".format(equip_newtonian[0].fV), "$\epsilon_e = {:.3f}$".format(equip_newtonian[0].epse), "$\epsilon_B = {:.3f}$".format(equip_newtonian[0].epsB), "---", "------"]
    
    for i in range(N_T):
        k = keys[i]
        v = ["${:.0f}".format(equip_newtonian[i].tdays) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_newtonian[i].Req())))        + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_newtonian[i].Req()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_newtonian[i].energyeq())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_newtonian[i].energyeq()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_newtonian[i].magField())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_newtonian[i].magField()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_newtonian[i].Ne())))         + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_newtonian[i].Ne()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_newtonian[i].CNMnumDens()))) + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_newtonian[i].CNMnumDens()), ddof = 1)) + "$",\
             "${:.3f}".format(np.mean(equip_newtonian[i].betaeqN()))            + "\\pm" + "{:.3f}".format(np.std(equip_newtonian[i].betaeqN(), ddof = 1)) + "$"]
        data[k] = v
    
    printlatextable(headers, data)
    
    headers[-1] = "$\\Gamma$"
    
    data = dict()
    
    keys = ["Jet", "$\\theta_{\\text{obs}}=$" + "${:.2f}$".format(equip_offAxis_thin[0].theta), "$f_A = {:.2f}$".format(equip_offAxis_thin[0].fA), "$f_V = {:.2f}$".format(equip_offAxis_thin[0].fV), "$\epsilon_e = {:.3f}$".format(equip_offAxis_thin[0].epse), "$\epsilon_B = {:.3f}$".format(equip_offAxis_thin[0].epsB), "---"]
    
    for i in range(N_T):
        k = keys[i]
        v = ["${:.0f}".format(equip_offAxis_thin[i].tdays) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].Req())))        + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].Req()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].energyeq())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].energyeq()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].magField())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].magField()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].Ne())))         + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].Ne()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].CNMnumDens()))) + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].CNMnumDens()), ddof = 1)) + "$",\
             "${:.3f}".format(np.mean(equip_offAxis_thin[i].gammaBulk()))            + "\\pm" + "{:.3f}".format(np.std(equip_offAxis_thin[i].gammaBulk(), ddof = 1)) + "$"]
        data[k] = v
    
    printlatextable(headers, data)
    
    data = dict()
    
    keys = ["Jet", "$\\theta_{\\text{obs}}=$" + "${:.2f}$".format(equip_offAxis_wide[0].theta), "$f_A = {:.2f}$".format(equip_offAxis_wide[0].fA), "$f_V = {:.2f}$".format(equip_offAxis_wide[0].fV), "$\epsilon_e = {:.3f}$".format(equip_offAxis_wide[0].epse), "$\epsilon_B = {:.3f}$".format(equip_offAxis_wide[0].epsB), "---"]
    
    for i in range(N_T):
        k = keys[i]
        v = ["${:.0f}".format(equip_offAxis_wide[i].tdays) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].Req())))        + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].Req()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].energyeq())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].energyeq()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].magField())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].magField()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].Ne())))         + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].Ne()), ddof = 1)) + "$",\
             "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].CNMnumDens()))) + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].CNMnumDens()), ddof = 1)) + "$",\
             "${:.3f}".format(np.mean(equip_offAxis_wide[i].gammaBulk()))            + "\\pm" + "{:.3f}".format(np.std(equip_offAxis_wide[i].gammaBulk(), ddof = 1)) + "$"]
        data[k] = v
    
    printlatextable(headers, data)
    
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

N = 1000

equip_newtonian    = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z,    0, p = p, epse = 0.1, epsB =   0.2, fA = 1, fV = 0.36, newtonian = True, isoNewtonianNe = True) for i in range(len(t))]
equip_offAxis_thin = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z, 1.05, p = p, epse = 0.1, epsB = 0.003, fA = 1, fV =    1, onAxis = False) for i in range(len(t))]
equip_offAxis_wide = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z, 1.57, p = p, epse = 0.1, epsB = 0.005, fA = 1, fV =    1, onAxis = False) for i in range(len(t))]

print("cosmology:", equip_newtonian[0].cosmo)

s =     "newtonian\nt      log(R)       log(E)       log(B)      log(Ne)      log(next)      beta\n"

for i in range(len(t)):
    s +="{}   {:.2f}±{:.2f}   {:.2f}±{:.2f}   {:.2f}±{:.2f}   {:.2f}±{:.2f}   {:.2f}±{:.2f}      {:.3f}±{:.3f}\n".format(\
        equip_newtonian[i].tdays,\
        np.mean(np.log10(equip_newtonian[i].Req())), np.std(np.log10(equip_newtonian[i].Req()), ddof = 1),\
        np.mean(np.log10(equip_newtonian[i].energyeq())), np.std(np.log10(equip_newtonian[i].energyeq()), ddof = 1),\
        np.mean(np.log10(equip_newtonian[i].magField())), np.std(np.log10(equip_newtonian[i].magField()), ddof = 1),\
        np.mean(np.log10(equip_newtonian[i].Ne())), np.std(np.log10(equip_newtonian[i].Ne()), ddof = 1),\
        np.mean(np.log10(equip_newtonian[i].CNMnumDens())), np.std(np.log10(equip_newtonian[i].CNMnumDens()), ddof = 1),\
        np.mean(equip_newtonian[i].betaeqN()), np.std(equip_newtonian[i].betaeqN(), ddof = 1))
        

s +=    "off axis theta = 1.05\nt      log(R)      log(E)       log(B)       log(Ne)      log(next)      Gamma\n"
 
for i in range(len(t)):
    s +="{}   {:.2f}±{:.2f}  {:.2f}±{:.2f}   {:.2f}±{:.2f}    {:.2f}±{:.2f}   {:.2f}±{:.2f}      {:.1f}±{:.1f}\n".format(\
        equip_offAxis_thin[i].tdays,\
        np.mean(np.log10(equip_offAxis_thin[i].Req())), np.std(np.log10(equip_offAxis_thin[i].Req()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_thin[i].energyeq())), np.std(np.log10(equip_offAxis_thin[i].energyeq()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_thin[i].magField())), np.std(np.log10(equip_offAxis_thin[i].magField()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_thin[i].Ne())), np.std(np.log10(equip_offAxis_thin[i].Ne()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_thin[i].CNMnumDens())), np.std(np.log10(equip_offAxis_thin[i].CNMnumDens()), ddof = 1),\
        np.mean(equip_offAxis_thin[i].gammaBulk()), np.std(equip_offAxis_thin[i].gammaBulk(), ddof = 1))
        
s +=    "off axis theta = 1.57\nt      log(R)       log(E)       log(B)       log(Ne)      log(next)      Gamma\n"

for i in range(len(t)):
    s +="{}   {:.2f}±{:.2f}  {:.2f}±{:.2f}   {:.2f}±{:.2f}    {:.2f}±{:.2f}   {:.2f}±{:.2f}      {:.1f}±{:.1f}\n".format(\
        equip_offAxis_wide[i].tdays,\
        np.mean(np.log10(equip_offAxis_wide[i].Req())), np.std(np.log10(equip_offAxis_wide[i].Req()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_wide[i].energyeq())), np.std(np.log10(equip_offAxis_wide[i].energyeq()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_wide[i].magField())), np.std(np.log10(equip_offAxis_wide[i].magField()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_wide[i].Ne())), np.std(np.log10(equip_offAxis_wide[i].Ne()), ddof = 1),\
        np.mean(np.log10(equip_offAxis_wide[i].CNMnumDens())), np.std(np.log10(equip_offAxis_wide[i].CNMnumDens()), ddof = 1),\
        np.mean(equip_offAxis_wide[i].gammaBulk()), np.std(equip_offAxis_wide[i].gammaBulk()), ddof = 1)

print(s)

with open("ChristyTable4.txt", "w") as text_file:
    text_file.write(s)
    
print()
printchristytable(equip_newtonian, equip_offAxis_thin, equip_offAxis_wide)