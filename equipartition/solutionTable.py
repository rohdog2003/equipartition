# -*- coding: utf-8 -*-
"""
Created on Wed May 14 13:20:11 2025

@author: rohdo
"""

import numpy as np
import scipy
import pickle

tol = 1e-15
maxiter = 500

xstep = 0.05
ystep = 0.1

t = np.arange(0, np.pi/2, xstep) # TODO should upper limit be np.pi? Physics says no math says yes.
tc = np.arange(1e-12, 0.001**(-17/24), ystep) # considers Newtonian velocity to be between 0.001c and less than 17000c

np.save("theta.npy", t)
np.save("thetac.npy", tc)

T, TC = np.meshgrid(t, tc, indexing = "ij")

def constraint_gammaBeta(u, t, tc):
    """MP23 (29) solved for zero and in terms of the four velocity"""
    if u == 1:
        return 1
    else:
        return (np.sqrt(1 + u**2) - u * np.cos(t)) * (1 - 1/(1 + u**2))**(-17/48) - tc

def solveMin(t, tc):
    """"""
    ucres = scipy.optimize.minimize(constraint_gammaBeta, 1, args = (t, tc))
    ucmin = ucres.x[0]
    
    return ucmin

@np.vectorize
def solveGammaBetaOn(t, tc, tol, maxiter):
    """"""    
    ucmin = solveMin(t, tc)
    fmin = constraint_gammaBeta(ucmin, t, tc)
        
    if (fmin < 0):
        gammaBetaOn = scipy.optimize.brentq(constraint_gammaBeta,\
                                            (1 + tc)**(-24/17), ucmin, args = (t, tc),\
                                            xtol = tol, maxiter = maxiter)
    else:
        gammaBetaOn = np.nan
    
    return gammaBetaOn
    
@np.vectorize
def solveGammaBetaOff(t, tc, tol, maxiter):
    """"""
    ucmin = solveMin(t, tc)
    fmin = constraint_gammaBeta(ucmin, t, tc)
    
    if (fmin < 0 and t != 0):
        gammaBetaOff = scipy.optimize.brentq(constraint_gammaBeta,\
                                             ucmin, tc/(1-np.cos(t)), args = (t, tc),\
                                             xtol = tol, maxiter = maxiter)
    else:
        gammaBetaOff = np.nan
        
    return gammaBetaOff

print("computing on axis solutions")
gammaBetaOn = solveGammaBetaOn(T, TC, tol, maxiter)
print("saving on axis solution")
np.save("onAxisSolutionPickle.npy", gammaBetaOn)

print("computing off axis solutions")
gammaBetaOff = solveGammaBetaOff(T, TC, tol, maxiter)
print("saving off axis solution")
np.save("offAxisSolutionPickle.npy", gammaBetaOff)
