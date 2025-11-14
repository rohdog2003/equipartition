# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 10:29:42 2025

@author: rohdo
"""
import numpy as np
from equipartition import Equipartition
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

def percentdifference(tseries1, tseries1err, tseries2, tseries2err):
    """"""
    pdiff = tseries1/tseries2
    pdifferr = np.sqrt((tseries1err/tseries2)**2 + (tseries2err * tseries1/tseries2**2)**2)
    
    return pdiff, pdifferr

colors = ["#006EB2", "#E69F00", "#750000"]

def percentdiffsubplot(subplotnum, title, xlabel, ylabel, t,\
                       tseries1a, tseries1aerr, tseries1b, tseries1berr, tseries1c, tseries1cerr,\
                       tseries2a, tseries2aerr, tseries2b, tseries2berr, tseries2c, tseries2cerr,\
                       plota = True, plotb = True, plotc = True, top = None, bottom = None, flipy = False, legend = False, labels = [], size = 6, loc = "best", ylog = False):
    """""" # TODO replace with colorblind friendly colors
    shiftfactor = 1.1
    alpha = 1
    
    ax = plt.subplot(subplotnum)
    if flipy:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    
    if ylog:
        plt.yscale("log")
    
    plt.xscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.axhline(1, linestyle = "-.", color = "k")
    pdiffa, pdiffaerr = percentdifference(tseries1a, tseries1aerr, tseries2a, tseries2aerr)
    if plota:
        plt.errorbar(t, pdiffa, yerr = pdiffaerr, marker = "s", alpha = alpha, linestyle = ":", mfc = "None", color = colors[0])
    pdiffb, pdiffberr = percentdifference(tseries1b, tseries1berr, tseries2b, tseries2berr)
    if plotb:
        plt.errorbar(t * shiftfactor, pdiffb, yerr = pdiffberr, marker = "d", alpha = alpha, linestyle = "--", mfc = "None", color = colors[1])
    pdiffc, pdiffcerr = percentdifference(tseries1c, tseries1cerr, tseries2c, tseries2cerr)
    if plotc:
        plt.errorbar(t * shiftfactor**2, pdiffc, yerr = pdiffcerr, marker = "o", alpha = alpha, linestyle = "-", mfc = "None", color = colors[2])
    
    if legend:
        plt.axhline(1, linestyle = "-.", color = "k")
        pdiffa, pdiffaerr = percentdifference(tseries1a, tseries1aerr, tseries2a, tseries2aerr)
        if plota:
            plt.errorbar(t, pdiffa, yerr = pdiffaerr, marker = "s", alpha = alpha, linestyle = ":", mfc = "None", color = colors[0], label = labels[0])
        pdiffb, pdiffberr = percentdifference(tseries1b, tseries1berr, tseries2b, tseries2berr)
        if plotb:
            plt.errorbar(t * shiftfactor, pdiffb, yerr = pdiffberr, marker = "d", alpha = alpha, linestyle = "--", mfc = "None", color = colors[1], label = labels[1])
        pdiffc, pdiffcerr = percentdifference(tseries1c, tseries1cerr, tseries2c, tseries2cerr)
        if plotc:
            plt.errorbar(t * shiftfactor**2, pdiffc, yerr = pdiffcerr, marker = "o", alpha = alpha, linestyle = "-", mfc = "None", color = colors[2], label = labels[2])
            
        plt.legend(prop={'size': size}, loc = loc)
        
    if top is not None:
        plt.ylim(top = top)
    if bottom is not None:
        plt.ylim(bottom = bottom)

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

N = 1000 # TODO

equip_newtonian    = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z,    0, p = p, epse = 0.1, epsB =   0.2, fA = 1, fV = 0.36, newtonian = True, isoNewtonianNe = True, factorsFour = True) for i in range(len(t))]
equip_offAxis_thin = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z, 1.05, p = p, epse = 0.1, epsB = 0.003, fA = 1, fV =    1, onAxis = False) for i in range(len(t))]
equip_offAxis_wide = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z, 1.57, p = p, epse = 0.1, epsB = 0.005, fA = 1, fV =    1, onAxis = False) for i in range(len(t))]

print("cosmology:", equip_newtonian[0].cosmo)
print("C:", equip_newtonian[0].C())

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

# reading in Christy table 4 data
t,\
logRgeom1, logRgeom1merr, logRgeom1perr,\
logEgeom1, logEgeom1merr, logEgeom1perr, _,\
logBgeom1, logBgeom1merr, logBgeom1perr,\
logNgeom1, logNgeom1merr, logNgeom1perr,\
logngeom1, logngeom1merr, logngeom1perr,\
betageom1, betageom1merr, betageom1perr = np.loadtxt("apjad675bt4_ChatGPTprocessed_geom1.txt", dtype = str).T

logRgeom1    = np.float64(logRgeom1)
logRgeom1err = np.maximum(np.abs(np.float64(logRgeom1merr)), np.float64(logRgeom1perr))
logEgeom1    = np.float64(logEgeom1)
logEgeom1err = np.maximum(np.abs(np.float64(logEgeom1merr)), np.float64(logEgeom1perr))
logBgeom1    = np.float64(logBgeom1)
logBgeom1err = np.maximum(np.abs(np.float64(logBgeom1merr)), np.float64(logBgeom1perr))
logNgeom1    = np.float64(logNgeom1)
logNgeom1err = np.maximum(np.abs(np.float64(logNgeom1merr)), np.float64(logNgeom1perr))
logngeom1    = np.float64(logngeom1)
logngeom1err = np.maximum(np.abs(np.float64(logngeom1merr)), np.float64(logngeom1perr))
betageom1    = np.float64(betageom1)
betageom1err = np.maximum(np.abs(np.float64(betageom1merr)), np.float64(betageom1perr))

t,\
logRgeom2, logRgeom2merr, logRgeom2perr,\
logEgeom2, logEgeom2merr, logEgeom2perr, _, _, _,\
logBgeom2, logBgeom2merr, logBgeom2perr,\
logNgeom2, logNgeom2merr, logNgeom2perr,\
logngeom2, logngeom2merr, logngeom2perr,\
Gammgeom2, Gammgeom2merr, Gammgeom2perr = np.loadtxt("apjad675bt4_ChatGPTprocessed_geom2.txt", dtype = str).T

logRgeom2    = np.float64(logRgeom2)
logRgeom2err = np.maximum(np.abs(np.float64(logRgeom2merr)), np.float64(logRgeom2perr))
logEgeom2    = np.float64(logEgeom2)
logEgeom2err = np.maximum(np.abs(np.float64(logEgeom2merr)), np.float64(logEgeom2perr))
logBgeom2    = np.float64(logBgeom2)
logBgeom2err = np.maximum(np.abs(np.float64(logBgeom2merr)), np.float64(logBgeom2perr))
logNgeom2    = np.float64(logNgeom2)
logNgeom2err = np.maximum(np.abs(np.float64(logNgeom2merr)), np.float64(logNgeom2perr))
logngeom2    = np.float64(logngeom2)
logngeom2err = np.maximum(np.abs(np.float64(logngeom2merr)), np.float64(logngeom2perr))
Gammgeom2    = np.float64(Gammgeom2)
Gammgeom2err = np.maximum(np.abs(np.float64(Gammgeom2merr)), np.float64(Gammgeom2perr))

t,\
logRgeom3, logRgeom3merr, logRgeom3perr,\
logEgeom3, logEgeom3merr, logEgeom3perr, _, _, _,\
logBgeom3, logBgeom3merr, logBgeom3perr,\
logNgeom3, logNgeom3merr, logNgeom3perr,\
logngeom3, logngeom3merr, logngeom3perr,\
Gammgeom3, Gammgeom3merr, Gammgeom3perr = np.loadtxt("apjad675bt4_ChatGPTprocessed_geom3.txt", dtype = str).T

logRgeom3    = np.float64(logRgeom3)
logRgeom3err = np.maximum(np.abs(np.float64(logRgeom3merr)), np.float64(logRgeom3perr))
logEgeom3    = np.float64(logEgeom3)
logEgeom3err = np.maximum(np.abs(np.float64(logEgeom3merr)), np.float64(logEgeom3perr))
logBgeom3    = np.float64(logBgeom3)
logBgeom3err = np.maximum(np.abs(np.float64(logBgeom3merr)), np.float64(logBgeom3perr))
logNgeom3    = np.float64(logNgeom3)
logNgeom3err = np.maximum(np.abs(np.float64(logNgeom3merr)), np.float64(logNgeom3perr))
logngeom3    = np.float64(logngeom3)
logngeom3err = np.maximum(np.abs(np.float64(logngeom3merr)), np.float64(logngeom3perr))
Gammgeom3    = np.float64(Gammgeom3)
Gammgeom3err = np.maximum(np.abs(np.float64(Gammgeom3merr)), np.float64(Gammgeom3perr))

fig = plt.figure(figsize = (6, 6))

def plotfigure():
    """"""
    percentdiffsubplot(321, "", "", "$R_{ours}/R_{C24}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].Req()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].Req(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_thin[i].Req()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].Req(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_wide[i].Req()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].Req(), ddof = 1) for i in range(len(tdays))]),\
                       10**logRgeom1, np.log(10) * 10**logRgeom1 * logRgeom1err,\
                       10**logRgeom2, np.log(10) * 10**logRgeom2 * logRgeom2err,\
                       10**logRgeom3, np.log(10) * 10**logRgeom3 * logRgeom3err, legend = True, labels = ["Newtonian", r"off axis $\theta=1.05$", r"off axis $\theta=1.57$"], size = 8, loc = "lower right")
    percentdiffsubplot(322, "", "", "$E_{ours}/E_{C24}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].energyeq()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].energyeq(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_thin[i].energyeq()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].energyeq(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_wide[i].energyeq()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].energyeq(), ddof = 1) for i in range(len(tdays))]),\
                       10**logEgeom1, np.log(10) * 10**logEgeom1 * logEgeom1err,\
                       10**logEgeom2, np.log(10) * 10**logEgeom2 * logEgeom2err,\
                       10**logEgeom3, np.log(10) * 10**logEgeom3 * logEgeom3err, flipy = True)
    percentdiffsubplot(323, "", "", "$B_{ours}/B_{C24}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].magField()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].magField(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_thin[i].magField()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].magField(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_wide[i].magField()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].magField(), ddof = 1) for i in range(len(tdays))]),\
                       10**logBgeom1, np.log(10) * 10**logBgeom1 * logBgeom1err,\
                       10**logBgeom2, np.log(10) * 10**logBgeom2 * logBgeom2err,\
                       10**logBgeom3, np.log(10) * 10**logBgeom3 * logBgeom3err)
    percentdiffsubplot(324, "", "", "$N_{e,ours}/N_{e,C24}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].Ne()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].Ne(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_thin[i].Ne()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].Ne(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_wide[i].Ne()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].Ne(), ddof = 1) for i in range(len(tdays))]),\
                       10**logNgeom1, np.log(10) * 10**logNgeom1 * logNgeom1err,\
                       10**logNgeom2, np.log(10) * 10**logNgeom2 * logNgeom2err,\
                       10**logNgeom3, np.log(10) * 10**logNgeom3 * logNgeom3err, flipy = True)
    percentdiffsubplot(325, "", r"Time $t$ [days]", "$n_{ext,ours}/n_{ext,C24}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].CNMnumDens()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].CNMnumDens(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_thin[i].CNMnumDens()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].CNMnumDens(), ddof = 1) for i in range(len(tdays))]),\
                       np.array([np.mean(equip_offAxis_wide[i].CNMnumDens()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].CNMnumDens(), ddof = 1) for i in range(len(tdays))]),\
                       10**logngeom1, np.log(10) * 10**logngeom1 * logngeom1err,\
                       10**logngeom2, np.log(10) * 10**logngeom2 * logngeom2err,\
                       10**logngeom3, np.log(10) * 10**logngeom3 * logngeom3err)
    #percentdiffsubplot(326, "", "", "$\Gamma_{ours}/\Gamma_{C24}$", tdays,\
    #                   np.array([np.nan]), np.array([np.nan]),\
    #                   np.array([np.mean(equip_offAxis_thin[i].gammaBulk()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].gammaBulk(), ddof = 1) for i in range(len(tdays))]),\
    #                   np.array([np.mean(equip_offAxis_wide[i].gammaBulk()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].gammaBulk(), ddof = 1) for i in range(len(tdays))]),\
    #                   np.array([np.nan]), np.array([np.nan]),\
    #                   Gammgeom2, Gammgeom2err,\
    #                   Gammgeom3, Gammgeom3err,\
    #                   plota = False)
    #percentdiffsubplot(327, "", "time (days)", r"$\beta_{ours}/\beta_{C24}$", tdays,\
    #                   np.array([np.mean(equip_newtonian[i].betaeqN()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].betaeqN(), ddof = 1) for i in range(len(tdays))]),\
    #                   np.array([np.nan]), np.array([np.nan]),\
    #                   np.array([np.nan]), np.array([np.nan]),\
    #                   betageom1, betageom1err,\
    #                   np.array([np.nan]), np.array([np.nan]),\
    #                   np.array([np.nan]), np.array([np.nan]),\
    #                   plotb = False, plotc = False)
    
    shiftfactor = 1.1
    alpha = 1
    ax6 = plt.subplot(326)
    ax6.set_xticks([])
    ax6.set_yticks([])
    for spine in ax6.spines.values():
        spine.set_visible(False)
    
    tseries1a, tseries1aerr, tseries1b, tseries1berr, tseries1c, tseries1cerr,\
    tseries2a, tseries2aerr, tseries2b, tseries2berr, tseries2c, tseries2cerr =\
        [np.array([np.mean(equip_newtonian[i].betaeqN()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].betaeqN(), ddof = 1) for i in range(len(tdays))]),\
         np.array([np.mean(equip_offAxis_thin[i].gammaBulk()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].gammaBulk(), ddof = 1) for i in range(len(tdays))]),\
         np.array([np.mean(equip_offAxis_wide[i].gammaBulk()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].gammaBulk(), ddof = 1) for i in range(len(tdays))]),\
         betageom1, betageom1err,\
         Gammgeom2, Gammgeom2err,\
         Gammgeom3, Gammgeom3err]
    
    
    # Match x-limits with main subplots
    all_axes = plt.gcf().get_axes()
    
    shiftfactor = 1.1
    alpha = 1

    # Remove ax6 completely, just create two stacked axes manually in the bottom-right space
    # Determine bottom-right position from neighboring subplots
    pos_left_bottom = all_axes[5].get_position()  # subplot 325
    width = all_axes[4].get_position().width * 1.23 # FIXME magic number
    height = all_axes[4].get_position().height * 1.33 # FIXME magic number
    yshift = -0.045 # FIXME magic number

    # Two stacked axes in the bottom-right (roughly matching subplot 326 position)
    inner_ax1 = fig.add_axes([pos_left_bottom.x0, pos_left_bottom.y0 + 0.5*height + yshift,
                              width, 0.5*height])
    inner_ax2 = fig.add_axes([pos_left_bottom.x0, pos_left_bottom.y0 + yshift,
                              width, 0.5*height])

    # Prepare the data
    tseries1a, tseries1aerr, tseries1b, tseries1berr, tseries1c, tseries1cerr,\
    tseries2a, tseries2aerr, tseries2b, tseries2berr, tseries2c, tseries2cerr =\
        [np.array([np.mean(equip_newtonian[i].betaeqN()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].betaeqN(), ddof=1) for i in range(len(tdays))]),\
         np.array([np.mean(equip_offAxis_thin[i].gammaBulk()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_thin[i].gammaBulk(), ddof=1) for i in range(len(tdays))]),\
         np.array([np.mean(equip_offAxis_wide[i].gammaBulk()) for i in range(len(tdays))]), np.array([np.std(equip_offAxis_wide[i].gammaBulk(), ddof=1) for i in range(len(tdays))]),\
         betageom1, betageom1err,\
         Gammgeom2, Gammgeom2err,\
         Gammgeom3, Gammgeom3err]

    pdiffa, pdiffaerr = percentdifference(tseries1a, tseries1aerr, tseries2a, tseries2aerr)
    pdiffb, pdiffberr = percentdifference(tseries1b, tseries1berr, tseries2b, tseries2berr)
    pdiffc, pdiffcerr = percentdifference(tseries1c, tseries1cerr, tseries2c, tseries2cerr)

    # Plot in stacked axes
    inner_ax1.yaxis.tick_right()
    inner_ax1.yaxis.set_label_position("right")
    inner_ax1.errorbar(tdays, pdiffa, yerr=pdiffaerr, marker="s", alpha=alpha, linestyle=":", mfc="None", color=colors[0])
    inner_ax1.axhline(1, linestyle="-.", color="k")
    inner_ax1.set_ylabel(r"$\beta_{ours}/\beta_{C24}$")
    inner_ax1.set_xscale("log")
    
    inner_ax2.yaxis.tick_right()
    inner_ax2.yaxis.set_label_position("right")
    inner_ax2.errorbar(tdays*shiftfactor, pdiffb, yerr=pdiffberr, marker="d", alpha=alpha, linestyle="--", mfc="None", color=colors[1])
    inner_ax2.errorbar(tdays*shiftfactor**2, pdiffc, yerr=pdiffcerr, marker="o", alpha=alpha, linestyle="-", mfc="None", color=colors[2])
    inner_ax2.axhline(1, linestyle="-.", color="k")
    inner_ax2.set_ylabel(r"$\Gamma_{ours}/\Gamma_{C24}$")
    inner_ax2.set_xlabel(r"Time $t$ [days]")
    inner_ax2.set_xscale("log")

    # Align x-limits with left-hand subplots
    inner_ax1.set_xlim(all_axes[0].get_xlim())
    inner_ax2.set_xlim(all_axes[4].get_xlim())

    # Force figure layout to remove extra padding
    plt.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.08, hspace=0, wspace=0)

    plt.savefig("percentdiff.png", dpi=750)
    plt.savefig("percentdiff.svg")
    plt.show()

plotfigure()