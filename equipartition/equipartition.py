# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 09:00:20 2025

@author: rohdo
"""

import scipy.constants as scont
import scipy
import numpy as np
import cosmoCalc

class Equipartition:
    """Equipartition calculator based on Matsumoto and Piran 2023 (MP23) 
    prescription (https://doi.org/10.1093/mnras/stad1269). cgs + Gauss units
    """
    def __init__(self, FpmJy, nup10, tdays, z, theta, R17 = None, nuM10 = 1,\
                 nuA10 = 1, fA = 1, fV = 1, fOmega = 1, p = 2, onAxis = True,\
                 newtonian = False, table = False, tol = 1e-15, maxiter = 500,\
                 epse = 0.1, epsB = 0.1 * 6/11, corr = True, BDfactor = False,\
                 gammaM_newtonian = 2, hotprotons = True):
        self.table = table
        self.tol = tol
        self.maxiter = maxiter
        
        if self.table:
            thetaTable = np.load("theta.npy")
            thetacTable = np.load("thetac.npy")
            gammaBetaOnTable = np.load("onAxisSolutionPickle.npy")
            gammaBetaOffTable = np.load("offAxisSolutionPickle.npy")
            self.gammaBetaOnAxisInterp = scipy.interpolate.RectBivariateSpline(thetaTable, thetacTable, gammaBetaOnTable, kx = 1, ky = 1)
            self.gammaBetaOffAxisInterp = scipy.interpolate.RectBivariateSpline(thetaTable, thetacTable, gammaBetaOffTable, kx = 1, ky = 1)
        
        self.FpmJy = FpmJy
        self.Fp = FpmJy * 1e-26 # cgs
        self.nup10 = nup10
        self.nup = nup10 * 1e10
        self.tdays = tdays
        self.tsec = 86400 * self.tdays # s
        self.z = z
        self.theta = theta
        self.nuM10 = nuM10
        self.nuM = nuM10 * 1e10
        self.nuA10 = nuA10
        self.nuA = nuA10 * 1e10
        self.fA = fA
        self.fV = fV
        self.fOmega = fOmega
        self.p = p
        self.onAxis = onAxis;
        self.newtonian = newtonian;
        # calculated values
        self.dL = cosmoCalc.lumDistLCDM(z) * 100 # cm
        self.dL28 = self.dL/1e28
        self.epse = epse
        self.epsB = epsB
        self.xi = (1 + 1/self.epse)**hotprotons # hot proton term
        self.corr = corr
        self.BDfactor = BDfactor
        self.gammaM_newtonian = gammaM_newtonian
        
        # constants
        self.c_cgs = scont.c * 100 # cm/s
        self.m_e_cgs = scont.m_e * 1000 # g
        self.q_e_cgs = 4.80320425e-10 # statcoulombs
        self.m_p_cgs = scont.m_p * 1000 # g
        self.gamBet = self.gammaBeta() # calculate gammaBeta as object is created to reduce runtime
    
        if R17 is None: # handling default radius as equipartition radius
            self.R17 = self.Req()/1e17 # TODO make sure this factor of 1e17 should be here
        else:
            self.R17 = R17
        
        self.R = self.R17 * 1e17 # cm
    
    def eta(self):
        """MP23 (7)"""
        return self.newtonian * 1 +\
               np.logical_not(self.newtonian) * ((self.nuM10 <= self.nuA10) * 1 +\
                                                 np.logical_not((self.nuM10 <= self.nuA10)) * self.nuM10/self.nuA10)
    # forgot an equal sign on the <= and it was the worst bug ever to find
                   
    def pbar(self): # TODO defaults assume that absorption frequency is equal taken here to be greater but self.p is np.nan by default
        """"""
        return np.where(self.nuA10 < self.nuM10, 2, self.p)
    
    def eps(self):
        """"""
        return 11/(2 * (self.pbar() + 1)) * self.epsB/self.epse
    
    def deltaD(self):
        """MP23 (1)"""
        beta = np.sqrt(1 - 1/self.gammaBulk()**2)
    
        return self.newtonian * 1 +\
               np.logical_not(self.newtonian) * 1/(self.gammaBulk() * (1 - beta * np.cos(self.theta)))
    
    def gammae(self):
        """MP23 (8)"""
        return (self.C()/3) * 5.2e2 * (self.FpmJy * self.dL28**2 * self.eta()**(5/3))/\
               (self.nup10 * (1 + self.z)**3) *\
               self.gammaBulk()**2/(self.fA * self.R17**2 * self.deltaD())
    
    def Ne(self): # FIXME much too small in the newtonian case 
        """MP23 (9)"""
        return (self.C()/3)**2 * 4.1e54 * (self.FpmJy**3 * self.dL28**6 * self.eta()**(10/3))/\
               (self.nup10**5 * (1 + self.z)**8) *\
               self.gammaBulk()**4/(self.fA**2 * self.R17**4 * self.deltaD()**4)
    
    def magField(self):
        """MP23 (10)"""
        return (3/self.C())**2 * 1.3e-2 * (self.nup10**5 * (1 + self.z)**7)/\
               (self.FpmJy**2 * self.dL28**4 * self.eta()**(10/3)) *\
               (self.fA**2 * self.R17**4 * self.deltaD())/self.gammaBulk()**4
    
    def CNMnumDens(self):
        """MP23 (13)"""
        return (3 * self.Ne() * self.gammaBulk()**2)/(self.fOmega * np.pi * self.R**3)
    
    def energye(self): 
        """MP23 (14)"""
        if self.corr == False:
            return self.m_e_cgs * self.c_cgs**2 * self.gammae() * self.gammaBulk() * self.Ne()
        else:
            return (self.gammaM() / self.gammae())**(2 - self.pbar()) * self.m_e_cgs * self.c_cgs**2 * self.gammae() * self.gammaBulk() * self.Ne()
        
    def energyB(self):
        """MP23 (15)"""
        return self.magField()**2/(8 * np.pi) * self.fV * self.R**3/self.gammaBulk()**2
    
    def energytot(self):
        """MP23 (16)"""
        if self.corr == False:
            return self.energye() + self.energyB()
        else:
            return self.xi * self.energye() + self.energyB()
    
    def Req(self):
        """MP23 (17)"""
        pb = self.pbar()
        
        if self.corr:
            return self.ReqN() * self.gammaBulk() * self.deltaD()**(-(pb + 5)/(2 * pb + 13))
        
        return self.ReqN() * self.gammaBulk() * self.deltaD()**(-7/17)
    
    def ReqN(self):
        """MP23 (18)"""
        pb = self.pbar()
        
        if self.corr == False:
            return (self.C()/3)**(7/17) * 1.9e17 * (self.FpmJy**(8/17) * self.dL28**(16/17) * self.eta()**(35/51))/\
                   (self.nup10 * (1 + self.z)**(25/17)) * self.fA**(-7/17) * self.fV**(-1/17)
        else:
            #return self.Fp**((pb + 6)/(2 * pb + 13)) * self.dL**(2 * (pb + 6)/(2 * pb + 13)) * (self.xi * self.eps() * self.gammaM()**(1 - pb) *\
            #       ((pb + 1) * self.C()**(pb + 5) * self.c_cgs * self.eta()**(5/3 * (pb + 5)))/\
            #       (2**(pb + 11) * 11 * np.sqrt(3) * np.pi**(pb + 7) * self.m_e_cgs**(pb + 6) * self.nup**(2 * pb + 13) * (1 + self.z)**(3 * pb + 19) * self.fA**(pb + 5) * self.fV))**(1/(2 * pb + 13))
    
            return (self.xi**(1/(2 * pb + 13)) * self.eps()**(1/(2 * pb + 13)) * self.gammaM()**((2 - pb)/(2 * pb + 13)) * (pb + 1)**(1/(2 * pb + 13)) * self.C()**((pb + 5)/(2 * pb + 13)) * self.c_cgs**(1/(2 * pb + 13)) * self.Fp**((pb + 6)/(2 * pb + 13)) * self.dL**(2 * (pb + 6)/(2 * pb + 13)) * self.eta()**(5/3 * (pb + 5)/(2 * pb + 13)))/\
                   (2**((pb + 2)/(2 * pb + 13)) * 11**(1/(2 * pb + 13)) * np.sqrt(3)**(1/(2 * pb + 13)) * np.pi**((pb + 7)/(2 * pb + 13)) * self.m_e_cgs**((pb + 6)/(2 * pb + 13)) * self.nup * (1 + self.z)**((3 * pb + 19)/(2 * pb + 13)) * self.fA**((pb + 5)/(2 * pb + 13)) * self.fV**(1/(2 * pb + 13)))
    
    def _ReqNtilde(self):
        """"""
        pb = self.pbar()
        
        return (self.xi**(1/(2 * pb + 13)) * self.eps()**(1/(2 * pb + 13)) * self.chie()**((2 - pb)/(2 * pb + 13)) * (pb + 1)**(1/(2 * pb + 13)) * self.C()**((pb + 5)/(2 * pb + 13)) * self.c_cgs**(1/(2 * pb + 13)) * self.Fp**((pb + 6)/(2 * pb + 13)) * self.dL**(2 * (pb + 6)/(2 * pb + 13)) * self.eta()**(5/3 * (pb + 5)/(2 * pb + 13)))/\
               (2**((pb + 2)/(2 * pb + 13)) * 11**(1/(2 * pb + 13)) * np.sqrt(3)**(1/(2 * pb + 13)) * np.pi**((pb + 7)/(2 * pb + 13)) * self.m_e_cgs**((pb + 6)/(2 * pb + 13)) * self.nup * (1 + self.z)**((3 * pb + 19)/(2 * pb + 13)) * self.fA**((pb + 5)/(2 * pb + 13)) * self.fV**(1/(2 * pb + 13)))
    
    def energyeq(self):
        """MP23 (19)"""
        pb = self.pbar()
        
        if self.corr:
            return self.energyeqN() * self.gammaBulk() * self.deltaD()**(-(7 * pb + 29)/(2 * pb + 13))
        
        return self.energyeqN() * self.gammaBulk() * self.deltaD()**(-43/17)
        
    def energyeqN(self):
        """MP23 (20)"""
        pb = self.pbar()
        
        if self.corr == False:
            return (self.C()/3)**(9/17) * 6.2e49 * (self.FpmJy**(20/17) * self.dL28**(40/17) * self.eta()**(15/17))/\
                   (self.nup10 * (1 + self.z)**(37/17)) * self.fA**(-9/17) * self.fV**(6/17)
        else:
            return (self.xi**(11/(2 * pb + 13)) * self.gammaM()**(11 * (2 - pb)/(2 * pb + 13)) * np.pi**((pb + 1)/(2 * pb + 13)) * 17 * self.C()**(3 * (pb + 1)/(2 * pb + 13)) * self.c_cgs**((4 * pb + 37)/(2 * pb + 13)) * self.m_e_cgs**((pb + 12)/(2 * pb + 13)) * self.Fp**((3 * pb + 14)/(2 * pb + 13)) * self.dL**(2 * (3 * pb + 14)/(2 * pb + 13)) * self.eta()**(5 * (pb + 1)/(2 * pb + 13)) * self.fV**(2 * (pb + 1)/(2 * pb + 13)))/\
                   (2**((7 * pb - 4)/(2 * pb + 13)) * (pb + 1)**(2 * (pb + 1)/(2 * pb + 13)) * 3**(5/(2 * pb + 13)) * 11**(11/(2 * pb + 13)) * np.sqrt(3)**(1/(2 * pb + 13)) * self.q_e_cgs**2 * self.nup * (1 + self.z)**((5 * pb + 27)/(2 * pb + 13)) * self.fA**(3 * (pb + 1)/(2 * pb + 13))) *\
                   (11/17 * self.eps()**(-2 * (pb + 1)/(2 * pb + 13)) + 2 * (pb + 1)/17 * self.eps()**(11/(2 * pb + 13)))
    
    def betaeqN(self):
        """MP23 (28)"""
        return (1 + self.z) * self.ReqN()/(self.c_cgs * self.tsec)
    
    def _betaeqNtilde(self):
        """"""
        return (1 + self.z) * self._ReqNtilde()/(self.c_cgs * self.tsec)
    
    def thetac(self):
        """MP23 (37)"""
        pb = self.pbar()
        
        return self.betaeqN()**(-(2 * pb + 13)/(3 * (pb + 6))) * self.eps()**(-1/(2 * pb + 13))
    
    def _thetactilde(self):
        """"""
        pb = self.pbar()
        
        return self._betaeqNtilde()**(-(2 * pb + 13)/(3 * (pb + 6))) * self.eps()**(-1/(2 * pb + 13))
    
    def constraint_gammaBeta(u, t, tc):
        """MP23 (29) solved for zero and in terms of the four velocity"""
        return (np.sqrt(1 + u**2) - u * np.cos(t)) * (1 - 1/(1 + u**2))**(-17/48) - tc
    
    def constraint_gammaBeta_corr(u, t, tctilde, pb):
        """"""
        return (np.sqrt(1 + u**2) - u * np.cos(t)) * (1 - 1/(1 + u**2))**(-17/48) - tctilde * (np.sqrt(1 + u**2) - 1)**((pb - 2)/(3 * (pb + 6)))
    
    def solveMin(t, tc):
        """"""
        ucres = scipy.optimize.minimize(Equipartition.constraint_gammaBeta, 1, args = (t, tc))
        ucmin = ucres.x[0]
        
        return ucmin
    
    def solveMin_corr(t, tctilde, pb):
        """"""
        ucres = scipy.optimize.minimize(Equipartition.constraint_gammaBeta_corr, 1, args = (t, tctilde, pb))
        ucmin = ucres.x[0]
        
        return ucmin
    
    def gammaBulk(self):
        """MP23 (34) and (39)"""
        return np.where(self.newtonian == True, 1, np.sqrt(1 + self.gamBet**2))
    
    @np.vectorize
    def solveGammaBetaOn(t, tc, tol, maxiter):
        """"""    
        ucmin = Equipartition.solveMin(t, tc)
        fmin = Equipartition.constraint_gammaBeta(ucmin, t, tc)
            
        if (fmin < 0):
            gammaBetaOn = scipy.optimize.brentq(Equipartition.constraint_gammaBeta,\
                                                (1 + tc)**(-24/17), ucmin, args = (t, tc),\
                                                xtol = tol, maxiter = maxiter)
        else:
            gammaBetaOn = np.nan
        
        return gammaBetaOn
    
    @np.vectorize
    def solveGammaBetaOn_corr(t, tctilde, pb, tol, maxiter):
        """"""
        ucmin = Equipartition.solveMin_corr(t, tctilde, pb)
        fmin = Equipartition.constraint_gammaBeta_corr(ucmin, t, tctilde, pb)
            
        if (fmin < 0):
            gammaBetaOn = scipy.optimize.brentq(Equipartition.constraint_gammaBeta_corr,\
                                                (1 + tctilde)**(-3 * (pb + 6)/(2 * pb + 13)), ucmin, args = (t, tctilde, pb),\
                                                xtol = tol, maxiter = maxiter)
        else:
            gammaBetaOn = np.nan
        
        return gammaBetaOn
    
    @np.vectorize
    def solveGammaBetaOff(t, tc, tol, maxiter):
        """"""
        ucmin = Equipartition.solveMin(t, tc)
        fmin = Equipartition.constraint_gammaBeta(ucmin, t, tc)
        
        if (fmin < 0 and t != 0):
            gammaBetaOff = scipy.optimize.brentq(Equipartition.constraint_gammaBeta,\
                                                 ucmin, tc/(1-np.cos(t)), args = (t, tc),\
                                                 xtol = tol, maxiter = maxiter)
        else:
            gammaBetaOff = np.nan
            
        return gammaBetaOff
    
    @np.vectorize
    def solveGammaBetaOff_corr(t, tctilde, pb, tol, maxiter):
        """"""
        ucmin = Equipartition.solveMin_corr(t, tctilde, pb)
        fmin = Equipartition.constraint_gammaBeta_corr(ucmin, t, tctilde, pb)
        
        if (fmin < 0 and t != 0):
            gammaBetaOff = scipy.optimize.brentq(Equipartition.constraint_gammaBeta_corr,\
                                                 ucmin, (tctilde/(1 - np.cos(t)))**(3 * (pb + 6)/(2 * (pb + 10))), args = (t, tctilde, pb),\
                                                 xtol = tol, maxiter = maxiter)
        else:
            gammaBetaOff = np.nan
            
        return gammaBetaOff
    
    def gammaBeta(self):
        """Bulk lorentz factor times the velocity fraction i.e. the four velocity"""
        if self.corr:
            if self.newtonian: # TODO vectorize but note that there is a weird circular dependency
                return self.betaeqN()
            elif self.onAxis:
                return Equipartition.solveGammaBetaOn_corr(self.theta, self._thetactilde(), self.pbar(), self.tol, self.maxiter)
            else:
                return Equipartition.solveGammaBetaOff_corr(self.theta, self._thetactilde(), self.pbar(), self.tol, self.maxiter)
            
            # AttributeError: 'Equipartition' object has no attribute 'gamBet'
            #return np.where(self.newtonian, self.betaeqN(),\
            #                np.where(self.onAxis, Equipartition.solveGammaBetaOn_corr(self.theta, self._thetactilde(), self.pbar(), self.tol, self.maxiter),\
            #                         Equipartition.solveGammaBetaOff_corr(self.theta, self._thetactilde(), self.pbar(), self.tol, self.maxiter)))
        
        if self.table: 
            return np.where(self.newtonian, self.betaeqN(),\
                            np.where(self.onAxis, np.vectorize(self.gammaBetaOnAxisInterp)(self.theta, self.thetac()),\
                                                                                           np.vectorize(self.gammaBetaOffAxisInterp)(self.theta, self.thetac())))
        
        else: # brentq method
            return np.where(self.newtonian, self.betaeqN(),\
                            np.where(self.onAxis, Equipartition.solveGammaBetaOn(self.theta, self.thetac(), self.tol, self.maxiter),\
                                     Equipartition.solveGammaBetaOff(self.theta, self.thetac(), self.tol, self.maxiter)))
    
    def chie(self):
        """"""
        return (self.p - 2)/(self.p - 1) * self.epse * self.m_p_cgs/self.m_e_cgs
    
    def gammaM(self): 
        """"""
        if self.newtonian:
            return self.gammaM_newtonian
        
        return self.chie() * (self.gammaBulk() - 1)
    
    def C(self):
        """Shen and Zhang (A15) (https://doi.org/10.1111/j.1365-2966.2009.15212.x).
        MP23 and Barniol Duran et al. 2013 (https://doi.org/10.1088/0004-637X/772/1/78)
        take this factor to be 3.
        """
        return np.where(self.BDfactor, 3,\
                        np.where(self.nuM10 <= self.nuA10, np.sqrt(2) * (self.p + 1) * (scipy.special.gamma((3 * self.p + 22)/12) * scipy.special.gamma((3 * self.p + 2)/12))/(scipy.special.gamma((3 * self.p + 19)/12) * scipy.special.gamma((3 * self.p - 1)/12)),\
                                 (self.p + 2) * (self.p - 1/3)/(self.p + 2/3)))
        