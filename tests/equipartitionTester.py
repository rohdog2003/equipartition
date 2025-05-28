# -*- coding: utf-8 -*-
"""
Created on Sat May 10 12:01:12 2025

@author: rohdo
"""

from equipartition import Equipartition
import unittest
import pandas as pd
import numpy as np
from cosmoCalc import lumDistLCDM
import matplotlib.pyplot as plt

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

theta = np.logspace(-0.6, 0.2, 19);

T = len(tdays)
N = len(theta)

class TestEquipartionMethods(unittest.TestCase):
    
    def setUp(self):
        """"""
        self.equip_Newtonian = Equipartition(FpmJy, nup10, tdays, z, 0, newtonian = True)
        self.equip_onAxis = [Equipartition(FpmJy, nup10, tdays, z, thet, newtonian = False, onAxis = True) for thet in theta]
        self.equip_offAxis = [Equipartition(FpmJy, nup10, tdays, z, thet, newtonian = False, onAxis = False) for thet in theta]
        
        self.equip_tableNewtonian = Equipartition(FpmJy, nup10, tdays, z, 0, newtonian = True, table = True)
        self.equip_tableOnAxis = [Equipartition(FpmJy, nup10, tdays, z, thet, newtonian = False, onAxis = True, table = True) for thet in theta]
        self.equip_tableOffAxis = [Equipartition(FpmJy, nup10, tdays, z, thet, newtonian = False, onAxis = False, table = True) for thet in theta]
    
    def test_zzzDebugPrint(self):
        print()
        
        print("theta:")
        for t in theta:
            print(t)
            
        print("betaeqN:")
        for b in self.equip_Newtonian.betaeqN():
            print(b)
        
        print("lumDistLCDM:")
        for l in lumDistLCDM(z):
            print(l)
            
        #u = np.linspace(0.001, 10, 100)
        #plt.plot(u, Equipartition.constraint_gammaBeta(u, 1.5848931924611136, 0.045719442431490356**(-17/24)))
        #plt.show()
        
        self.fail("debug print")
    
    def test_ReqN(self):
        """checking that the value of the Newtonian equipartition radius is the
        same as that calculated via calculator.
        """
        expected = []
            
        result = self.equip_Newtonian.ReqN()
        
        self.fail("not yet implemented")
    
        #for i in range(T):
        #    with self.subTest():
        #        self.assertEqual(expected[i], round(result[i], ))
    
    def test_ReqNsameOnAxis(self):
        """checking that the Newtonian equipartition radius is the same 
        irregardless of on axis and angle.
        """
        expected = self.equip_Newtonian.ReqN()
        result = self.equip_onAxis[-1].ReqN()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])
    
    def test_ReqNsameOffAxis(self):
        """checking that the Newtonian equipartition radius is the same 
        irregardless of off axis and angle.
        """
        expected = self.equip_Newtonian.ReqN()
        result = self.equip_offAxis[-1].ReqN()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])
    
    def test_energyeqNsameOnAxis(self):
        """checking that the Newtonian equipartition radius is the same 
        irregardless of on axis and angle.
        """
        expected = self.equip_Newtonian.energyeqN()
        result = self.equip_onAxis[-1].energyeqN()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])
                
    def test_energyeqNsameOffAxis(self):
        """checking that the Newtonian equipartition radius is the same 
        irregardless of off axis and angle.
        """
        expected = self.equip_Newtonian.energyeqN()
        result = self.equip_offAxis[-1].energyeqN()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])    
    
    def test_betaeqN(self):
        """checking that the value of the Newtonian velocity is the same as 
        that calculated via calculator.
        """
        expected = []
        result = self.equip_Newtonian.betaeqN()
        
        self.fail("not yet implemented")
    
    def test_betaeqNsameOnAxis(self):
        """checking that the Newtonian velocity is the same irregardless of 
        on axis and angle.
        """
        expected = self.equip_Newtonian.betaeqN()
        result = self.equip_onAxis[-1].betaeqN()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])
                
    def test_betaeqNsameOffAxis(self):
        """checking that the Newtonian velocity is the same irregardless of 
        off axis and angle.
        """
        expected = self.equip_Newtonian.betaeqN()
        result = self.equip_offAxis[-1].betaeqN()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])
    
    def test_thetac(self):
        """checking that the critical angle matches values calculated via a 
        calculator.
        """
        expected = [8.89396012234,\
                    8.45126503466,\
                    7.70300186453,\
                    8.33199414175,\
                    9.98518111595,\
                    6.66942166713,\
                    8.26859082922,\
                    10.1243860666]
        
        result = self.equip_Newtonian.thetac()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(round(expected[i], 9), round(result[i], 9))
    
    def test_gammaBetaNewtonianIsBeta(self):
        """checking that for the Newtonian case that the four velocity is the 
        velocity.
        """
        expected = self.equip_Newtonian.betaeqN()
        result = self.equip_Newtonian.gammaBeta()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], result[i])
        
    def test_gammaBetaOnAxisSolvesConstraint(self):
        """checking that the on axis solution solves the implicit equation for 
        the four velocity.
        """
        u = [equip.gammaBeta() for equip in self.equip_onAxis]
        tc = [equip.thetac() for equip in self.equip_onAxis]
        
        for i in range(N):
            for j in range(T):
                with self.subTest():
                    self.assertEqual(0, round(abs(Equipartition.constraint_gammaBeta(u[i][j], theta[i], tc[i][j])), 12))
        
    def test_gammaBetaOffAxisSolvesConstraint(self):
        """checking that the off axis solution solves the implicit equation for 
        the four velocity.
        """
        u = [equip.gammaBeta() for equip in self.equip_offAxis]
        tc = [equip.thetac() for equip in self.equip_offAxis]
        
        for i in range(N):
            for j in range(T):
                with self.subTest():
                    self.assertEqual(0, round(abs(Equipartition.constraint_gammaBeta(u[i][j], theta[i], tc[i][j])), 12))
    
    def test_gammaBetaOnAxis(self):
        """checks whether the on axis solutions are correctly solved for 
        in the implicit expression by comparing values achieved via a calculator 
        for the largest angle considered. Assumes that the correct Newtonian 
        velocity is found.
        """
        expected = [0.04588,\
                    0.04933,\
                    0.05628,\
                    0.05034,\
                    0.03893,\
                    0.06913,\
                    0.05089,\
                    0.03817]
            
        result = self.equip_onAxis[-1].gammaBeta()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], round(result[i], 5))
    
    def test_gammaBetaTableOnAxis(self):
        """checks whether the on axis solutions are correctly solved for 
        in the implicit expression by comparing values achieved via a calculator 
        for the largest angle considered. Assumes that the correct Newtonian 
        velocity is found.
        """
        expected = [0.04588,\
                    0.04933,\
                    0.05628,\
                    0.05034,\
                    0.03893,\
                    0.06913,\
                    0.05089,\
                    0.03817]
            
        result = self.equip_tableOnAxis[-1].gammaBeta()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], round(result[i], 5))
    
    def test_gammaBetaOffAxis(self):
        """checks whether the off axis solutions are correctly solved for 
        in the implicit expression by comparing values achieved via a calculator 
        for the largest angle considered. Assumes that the correct Newtonian 
        velocity is found.
        """
        expected = [8.67274,\
                    8.23097,\
                    7.48285,\
                    8.11185,\
                    9.75964,\
                    6.4455,\
                    8.04851,\
                    9.89812]
            
        result = self.equip_offAxis[-1].gammaBeta()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], round(result[i], 5))
                
    def test_gammaBetaTableOffAxis(self):
        """checks whether the off axis solutions are correctly solved for 
        in the implicit expression by comparing values achieved via a calculator 
        for the largest angle considered. Assumes that the correct Newtonian 
        velocity is found.
        """
        expected = [8.67274,\
                    8.23097,\
                    7.48285,\
                    8.11185,\
                    9.75964,\
                    6.4455,\
                    8.04851,\
                    9.89812]
            
        result = self.equip_tableOffAxis[-1].gammaBeta()
        
        for i in range(T):
            with self.subTest():
                self.assertEqual(expected[i], round(result[i], 5))
    
if __name__ == '__main__':
    unittest.main()