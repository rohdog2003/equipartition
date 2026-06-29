# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 10:29:42 2025

@author: rohdo
"""
import numpy as np
from equipartition import Equipartition
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.constants as scont
import scipy

def printlatextable(headers, data):
    """"""
    textabular = f"l|{'c'*len(headers)}"
    texheader = " & " + " & ".join(headers) + "\\\\"
    texdata = "\\hline\n"
    for label in data:
       texdata += f"{label} & {' & '.join(map(str,data[label]))} \\\\\n"
    
    print("\\begin{tabular}{"+textabular+"}")
    print(texheader)
    print(texdata,end="")
    print("\\end{tabular}")

def printtable(equip_newtonian):
    """"""
    headers = ["$t$", "$\\log(R)$ (cm)", "$\\log(E)$ (ergs)", "$\\log(B)$ (G)", "$\\log(N_e)$", "$\\log(n_{\\text{ext}})$ ($\\text{cm}^{-3}$)", "$\\beta$"]
    
    N_T = len(equip_newtonian); print(N_T)
    
    data = dict()
    
    keys = ["Spherical", "$f_A = {:.2f}$".format(equip_newtonian[0].fA), "$f_V = {:.2f}$".format(equip_newtonian[0].fV), "$\epsilon_e = {:.3f}$".format(equip_newtonian[0].epse), "$\epsilon_B = {:.3f}$".format(equip_newtonian[0].epsB), "---", "----", "-----", "------", "-------"]; print(len(keys))
    
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
    
    #headers[-1] = "$\\Gamma$"
    
    #data = dict()
    
# =============================================================================
#     keys = ["Jet", "$\\theta_{\\text{obs}}=$" + "${:.2f}$".format(equip_offAxis_thin[0].theta), "$f_A = {:.2f}$".format(equip_offAxis_thin[0].fA), "$f_V = {:.2f}$".format(equip_offAxis_thin[0].fV), "$\epsilon_e = {:.3f}$".format(equip_offAxis_thin[0].epse), "$\epsilon_B = {:.3f}$".format(equip_offAxis_thin[0].epsB), "---"]
#     
#     for i in range(N_T):
#         k = keys[i]
#         v = ["${:.0f}".format(equip_offAxis_thin[i].tdays) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].Req())))        + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].Req()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].energyeq())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].energyeq()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].magField())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].magField()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].Ne())))         + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].Ne()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_thin[i].CNMnumDens()))) + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_thin[i].CNMnumDens()), ddof = 1)) + "$",\
#              "${:.3f}".format(np.mean(equip_offAxis_thin[i].gammaBulk()))            + "\\pm" + "{:.3f}".format(np.std(equip_offAxis_thin[i].gammaBulk(), ddof = 1)) + "$"]
#         data[k] = v
#     
#     printlatextable(headers, data)
#     
#     data = dict()
#     
#     keys = ["Jet", "$\\theta_{\\text{obs}}=$" + "${:.2f}$".format(equip_offAxis_wide[0].theta), "$f_A = {:.2f}$".format(equip_offAxis_wide[0].fA), "$f_V = {:.2f}$".format(equip_offAxis_wide[0].fV), "$\epsilon_e = {:.3f}$".format(equip_offAxis_wide[0].epse), "$\epsilon_B = {:.3f}$".format(equip_offAxis_wide[0].epsB), "---"]
#     
#     for i in range(N_T):
#         k = keys[i]
#         v = ["${:.0f}".format(equip_offAxis_wide[i].tdays) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].Req())))        + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].Req()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].energyeq())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].energyeq()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].magField())))   + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].magField()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].Ne())))         + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].Ne()), ddof = 1)) + "$",\
#              "${:.2f}".format(np.mean(np.log10(equip_offAxis_wide[i].CNMnumDens()))) + "\\pm" + "{:.2f}".format(np.std(np.log10(equip_offAxis_wide[i].CNMnumDens()), ddof = 1)) + "$",\
#              "${:.3f}".format(np.mean(equip_offAxis_wide[i].gammaBulk()))            + "\\pm" + "{:.3f}".format(np.std(equip_offAxis_wide[i].gammaBulk(), ddof = 1)) + "$"]
#         data[k] = v
#     
#     printlatextable(headers, data)
# =============================================================================

def percentdifference(tseries1, tseries1err, tseries2, tseries2err):
    """"""
    pdiff = tseries1/tseries2
    pdifferr = np.sqrt((tseries1err/tseries2)**2 + (tseries2err * tseries1/tseries2**2)**2)
    
    return pdiff, pdifferr

colors = ["#006EB2", "#E69F00", "#750000"]

def percentdiffsubplot(subplotnum, title, xlabel, ylabel, t,\
                       tseries1a, tseries1aerr,\
                       tseries2a, tseries2aerr,\
                       plota = True, plotb = True, plotc = True, top = None, bottom = None, flipy = False, legend = False, labels = [], size = 6, loc = "best", ylog = False):
    """""" # TODO replace with colorblind friendly colors
    shiftfactor = 1.1
    alpha = 1
    
    ax = plt.subplot(subplotnum)
    if flipy:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.yaxis.labelpad = 20
    
    if ylog:
        plt.yscale("log")
    
    plt.xscale("log")
    plt.xlim(45, 650)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.axhline(1, linestyle = "--", color = "gray")
    pdiffa, pdiffaerr = percentdifference(tseries1a, tseries1aerr, tseries2a, tseries2aerr)
    if plota:
        plt.errorbar(t, pdiffa, yerr = pdiffaerr, marker = "s", alpha = alpha, linestyle = ":", mfc = "None", color = colors[0])
# =============================================================================
#     pdiffb, pdiffberr = percentdifference(tseries1b, tseries1berr, tseries2b, tseries2berr)
#     if plotb:
#         plt.errorbar(t * shiftfactor, pdiffb, yerr = pdiffberr, marker = "d", alpha = alpha, linestyle = "--", mfc = "None", color = colors[1])
#     pdiffc, pdiffcerr = percentdifference(tseries1c, tseries1cerr, tseries2c, tseries2cerr)
#     if plotc:
#         plt.errorbar(t * shiftfactor**2, pdiffc, yerr = pdiffcerr, marker = "o", alpha = alpha, linestyle = "-", mfc = "None", color = colors[2])
#     
# =============================================================================
    if legend:
        pdiffa, pdiffaerr = percentdifference(tseries1a, tseries1aerr, tseries2a, tseries2aerr)
        if plota:
            plt.errorbar(t, pdiffa, yerr = pdiffaerr, marker = "s", alpha = alpha, linestyle = ":", mfc = "None", color = colors[0], label = labels[0])
# =============================================================================
#         pdiffb, pdiffberr = percentdifference(tseries1b, tseries1berr, tseries2b, tseries2berr)
#         if plotb:
#             plt.errorbar(t * shiftfactor, pdiffb, yerr = pdiffberr, marker = "d", alpha = alpha, linestyle = "--", mfc = "None", color = colors[1], label = labels[1])
#         pdiffc, pdiffcerr = percentdifference(tseries1c, tseries1cerr, tseries2c, tseries2cerr)
#         if plotc:
#             plt.errorbar(t * shiftfactor**2, pdiffc, yerr = pdiffcerr, marker = "o", alpha = alpha, linestyle = "-", mfc = "None", color = colors[2], label = labels[2])
#             
# =============================================================================
        plt.legend(prop={'size': size}, loc = loc)
        
    if top is not None:
        plt.ylim(top = top)
    if bottom is not None:
        plt.ylim(bottom = bottom)

def logerrors2linearerrors(logv, logv_errl, logv_errh):
    """"""
    return 10**logv * np.abs(logv_errl), 10**logv * np.abs(logv_errh)

# deltat, Fp, Fp_err_lo, Fp_err_hi, log(nu_p), err_lo, err_hi, log(R_eq), err_lo, err_hi, log(E_eq), err_lo, err_hi, log(B), err_lo, err_hi, log(N_e), err_lo, err_hi, log(n_ext), err_lo, err_hi, log(v_eq), err_lo, err_hi
t,\
FpmJy, FpmJy_errl, FpmJy_errh,\
lognup, lognup_errl, lognup_errh,\
logReq, logReq_errl, logReq_errh,\
logEeq, logEeq_errl, logEeq_errh,\
logB, logB_errl, logB_errh,\
logNe, logNe_errl, logNe_errh_,\
lognext, lognext_errl, lognext_errh,\
logveq, logveq_errl, logveq_errh = np.loadtxt("Cendes2021table2.txt").T

nup = 10**lognup
nup_errl, nup_errh = logerrors2linearerrors(lognup, lognup_errl, lognup_errh)
nup_err = np.maximum(nup_errl, nup_errh)

FpmJy_err = np.maximum(-FpmJy_errl, FpmJy_errh)

t = np.float64(t); print("day:", t)
FpmJy = np.float64(FpmJy); print("peak flux (mJy):", FpmJy)
FpmJy_err = np.float64(FpmJy_err); #print(FpmJy_err)

tdays = t
nup10 = nup/1e10
nup10_err = nup_err/1e10

z = 0.051
p = 2.7 # ?
p_err = 0.2

N = 5000 # TODO change to large number

equip_newtonian    = [Equipartition(scipy.stats.truncnorm.rvs(0, np.inf, FpmJy[i], FpmJy_err[i], N), scipy.stats.truncnorm.rvs(2, 3, nup10[i], nup10_err[i], N), tdays[i], z,    0, p = scipy.stats.truncnorm.rvs(0, np.inf, p, p_err, N), epse = 0.1, epsB =   0.02, fA = 1, fV = 0.36, fOmega = 4, newtonian = True, isoNewtonianNe = False, factorsFour = False, corr = True, numelectrons = True, hotprotons = True, outofequipartition = True, energysum = True) for i in range(len(t))]
#equip_offAxis_thin = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z, 1.05, p = p, epse = 0.1, epsB = 0.003, fA = 1, fV =    1, onAxis = False) for i in range(len(t))]
#equip_offAxis_wide = [Equipartition(np.random.normal(FpmJy[i], FpmJy_err[i], N), np.random.normal(nup10[i], nup10_err[i], N), tdays[i], z, 1.57, p = p, epse = 0.1, epsB = 0.005, fA = 1, fV =    1, onAxis = False) for i in range(len(t))]

print("cosmology:", equip_newtonian[0].cosmo)
print("C:", equip_newtonian[0].C())
print("gammae Newtonian:", [np.mean(equip.gammae()) for equip in equip_newtonian])
print("gammaM Newtonian:", [np.mean(equip.gammaM()) for equip in equip_newtonian])
#print("gammae offAxis theta=1.05:", [np.mean(equip.gammae()) for equip in equip_offAxis_thin])
#print("gammaM offAxis theta=1.05:", [np.mean(equip.gammaM()) for equip in equip_offAxis_thin])
#print("gammae offAxis theta=1.57:", [np.mean(equip.gammae()) for equip in equip_offAxis_wide])
#print("gammaM offAxis theta=1.57:", [np.mean(equip.gammaM()) for equip in equip_offAxis_wide])

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
        

#s +=    "off axis theta = 1.05\nt      log(R)      log(E)       log(B)       log(Ne)      log(next)      Gamma\n"
 
#for i in range(len(t)):
#    s +="{}   {:.2f}±{:.2f}  {:.2f}±{:.2f}   {:.2f}±{:.2f}    {:.2f}±{:.2f}   {:.2f}±{:.2f}      {:.1f}±{:.1f}\n".format(\
#        equip_offAxis_thin[i].tdays,\
#        np.mean(np.log10(equip_offAxis_thin[i].Req())), np.std(np.log10(equip_offAxis_thin[i].Req()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_thin[i].energyeq())), np.std(np.log10(equip_offAxis_thin[i].energyeq()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_thin[i].magField())), np.std(np.log10(equip_offAxis_thin[i].magField()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_thin[i].Ne())), np.std(np.log10(equip_offAxis_thin[i].Ne()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_thin[i].CNMnumDens())), np.std(np.log10(equip_offAxis_thin[i].CNMnumDens()), ddof = 1),\
#        np.mean(equip_offAxis_thin[i].gammaBulk()), np.std(equip_offAxis_thin[i].gammaBulk(), ddof = 1))
        
#s +=    "off axis theta = 1.57\nt      log(R)       log(E)       log(B)       log(Ne)      log(next)      Gamma\n"
#
#for i in range(len(t)):
#    s +="{}   {:.2f}±{:.2f}  {:.2f}±{:.2f}   {:.2f}±{:.2f}    {:.2f}±{:.2f}   {:.2f}±{:.2f}      {:.1f}±{:.1f}\n".format(\
#        equip_offAxis_wide[i].tdays,\
#        np.mean(np.log10(equip_offAxis_wide[i].Req())), np.std(np.log10(equip_offAxis_wide[i].Req()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_wide[i].energyeq())), np.std(np.log10(equip_offAxis_wide[i].energyeq()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_wide[i].magField())), np.std(np.log10(equip_offAxis_wide[i].magField()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_wide[i].Ne())), np.std(np.log10(equip_offAxis_wide[i].Ne()), ddof = 1),\
#        np.mean(np.log10(equip_offAxis_wide[i].CNMnumDens())), np.std(np.log10(equip_offAxis_wide[i].CNMnumDens()), ddof = 1),\
#        np.mean(equip_offAxis_wide[i].gammaBulk()), np.std(equip_offAxis_wide[i].gammaBulk()), ddof = 1)

#print(s)

with open("CendesTable2.txt", "w") as text_file:
    text_file.write(s)
    
print()
printtable(equip_newtonian)

logbetaeq = 10**logveq/scont.c
logbetaeq_errl, logbetaeq_errh = logerrors2linearerrors(logveq, logveq_errl, logveq_errh)
logbetaeq_errl /= scont.c
logbetaeq_errh /= scont.c


# redefinitions
logRgeom1, logRgeom1merr, logRgeom1perr,\
logEgeom1, logEgeom1merr, logEgeom1perr,\
logBgeom1, logBgeom1merr, logBgeom1perr,\
logNgeom1, logNgeom1merr, logNgeom1perr,\
logngeom1, logngeom1merr, logngeom1perr,\
betageom1, betageom1merr, betageom1perr =\
logReq, logReq_errl, logReq_errh,\
logEeq, logEeq_errl, logEeq_errh,\
logB, logB_errl, logB_errh,\
logNe, logNe_errl, logNe_errh_,\
lognext, lognext_errl, lognext_errh,\
logbetaeq, logbetaeq_errl, logbetaeq_errh

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

# =============================================================================
# t,\
# logRgeom2, logRgeom2merr, logRgeom2perr,\
# logEgeom2, logEgeom2merr, logEgeom2perr, _, _, _,\
# logBgeom2, logBgeom2merr, logBgeom2perr,\
# logNgeom2, logNgeom2merr, logNgeom2perr,\
# logngeom2, logngeom2merr, logngeom2perr,\
# Gammgeom2, Gammgeom2merr, Gammgeom2perr = np.loadtxt("", dtype = str).T
# 
# logRgeom2    = np.float64(logRgeom2)
# logRgeom2err = np.maximum(np.abs(np.float64(logRgeom2merr)), np.float64(logRgeom2perr))
# logEgeom2    = np.float64(logEgeom2)
# logEgeom2err = np.maximum(np.abs(np.float64(logEgeom2merr)), np.float64(logEgeom2perr))
# logBgeom2    = np.float64(logBgeom2)
# logBgeom2err = np.maximum(np.abs(np.float64(logBgeom2merr)), np.float64(logBgeom2perr))
# logNgeom2    = np.float64(logNgeom2)
# logNgeom2err = np.maximum(np.abs(np.float64(logNgeom2merr)), np.float64(logNgeom2perr))
# logngeom2    = np.float64(logngeom2)
# logngeom2err = np.maximum(np.abs(np.float64(logngeom2merr)), np.float64(logngeom2perr))
# Gammgeom2    = np.float64(Gammgeom2)
# Gammgeom2err = np.maximum(np.abs(np.float64(Gammgeom2merr)), np.float64(Gammgeom2perr))
# 
# t,\
# logRgeom3, logRgeom3merr, logRgeom3perr,\
# logEgeom3, logEgeom3merr, logEgeom3perr, _, _, _,\
# logBgeom3, logBgeom3merr, logBgeom3perr,\
# logNgeom3, logNgeom3merr, logNgeom3perr,\
# logngeom3, logngeom3merr, logngeom3perr,\
# Gammgeom3, Gammgeom3merr, Gammgeom3perr = np.loadtxt("", dtype = str).T
# 
# logRgeom3    = np.float64(logRgeom3)
# logRgeom3err = np.maximum(np.abs(np.float64(logRgeom3merr)), np.float64(logRgeom3perr))
# logEgeom3    = np.float64(logEgeom3)
# logEgeom3err = np.maximum(np.abs(np.float64(logEgeom3merr)), np.float64(logEgeom3perr))
# logBgeom3    = np.float64(logBgeom3)
# logBgeom3err = np.maximum(np.abs(np.float64(logBgeom3merr)), np.float64(logBgeom3perr))
# logNgeom3    = np.float64(logNgeom3)
# logNgeom3err = np.maximum(np.abs(np.float64(logNgeom3merr)), np.float64(logNgeom3perr))
# logngeom3    = np.float64(logngeom3)
# logngeom3err = np.maximum(np.abs(np.float64(logngeom3merr)), np.float64(logngeom3perr))
# Gammgeom3    = np.float64(Gammgeom3)
# Gammgeom3err = np.maximum(np.abs(np.float64(Gammgeom3merr)), np.float64(Gammgeom3perr))
# =============================================================================

fig = plt.figure(figsize = (6, 6))

def plotfigure():
    """"""
    percentdiffsubplot(321, "", "", r"$R_{\rm{eq,ours}}/R_{\rm{eq,BNP13}}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].Req()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].Req(), ddof = 1) for i in range(len(tdays))]),\
                       10**logRgeom1, np.log(10) * 10**logRgeom1 * logRgeom1err,\
                       legend = True, labels = ["Newtonian"], size = 8, loc = "lower right", ylog = False)
    percentdiffsubplot(322, "", "", r"$E_{\rm{eq,ours}}/E_{\rm{eq,BNP13}}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].energyeq()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].energyeq(), ddof = 1) for i in range(len(tdays))]),\
                       10**logEgeom1, np.log(10) * 10**logEgeom1 * logEgeom1err,\
                       flipy = True, ylog = False)
    percentdiffsubplot(323, "", "", r"$B_{\rm{ours}}/B_{\rm{BNP13}}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].magField()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].magField(), ddof = 1) for i in range(len(tdays))]),\
                       10**logBgeom1, np.log(10) * 10**logBgeom1 * logBgeom1err,\
                       ylog = False)
    percentdiffsubplot(324, "", "", r"$N_{\rm{e,ours}}/N_{\rm{e,BNP13}}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].Ne()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].Ne(), ddof = 1) for i in range(len(tdays))]),\
                       10**logNgeom1, np.log(10) * 10**logNgeom1 * logNgeom1err,\
                       flipy = True, ylog = False)
    percentdiffsubplot(325, "", r"Time $t$ [days]", r"$n_{\rm{ext,ours}}/n_{\rm{ext,BNP13}}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].CNMnumDens()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].CNMnumDens(), ddof = 1) for i in range(len(tdays))]),\
                       10**logngeom1, np.log(10) * 10**logngeom1 * logngeom1err,\
                       ylog = False)
    percentdiffsubplot(326, "", "Time $t$ [days]", r"$\beta_{\rm{ours}}/\beta_{\rm{BNP13}}$", tdays,\
                       np.array([np.mean(equip_newtonian[i].betaeqN()) for i in range(len(tdays))]), np.array([np.std(equip_newtonian[i].betaeqN())]),\
                       betageom1, betageom1err,\
                       flipy = True)
    
    plt.subplots_adjust(wspace=0, hspace=0)
        
    plt.savefig("percentdiffCendes.png", dpi=750)
    plt.savefig("percentdiffCendes.svg")
    plt.show()

plotfigure()