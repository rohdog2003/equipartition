# -*- coding: utf-8 -*-
"""
Created on Wed May 14 15:51:13 2025

@author: rohdo
"""
import numpy as np
import scipy
import timeit

t = np.load("theta.npy")
tc = np.load("thetac.npy")

gammaBetaOn = np.load("onAxisSolutionPickle.npy")
gammaBetaOff = np.load("offAxisSolutionPickle.npy")

onAxisInterp = scipy.interpolate.RectBivariateSpline(t, tc, gammaBetaOn, kx = 1, ky = 1)
print("onAxis expect:", 0.17186)
print("onAxis result:", onAxisInterp(np.pi/4, np.pi)[0][0])

offAxisInterp = scipy.interpolate.RectBivariateSpline(t, tc, gammaBetaOff, kx = 1, ky = 1)
print("offAxis expect:", 10.53027)
print("offAxis result:", offAxisInterp(np.pi/4, np.pi)[0][0])

print("timing for a large call")
start = timeit.timeit()
a = np.vectorize(offAxisInterp)(np.pi/4, tc)
end = timeit.timeit()
print(a)
print("elapsed time:", end - start)