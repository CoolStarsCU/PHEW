#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from astropy import log

"""
Calculate the line depth of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam

Args:
----------
xmin,xmax - the specified interval of the spectrum to plot 
excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
n - the number of Monte Carlo iterations 

Returns:
-------
- the mean and standard deviation of the line depth measured n times
- the spectrum plotted with the pseudo-continuum (yellow)
- a histogram of the line depth distribution
"""

def measure_line_depth(filename,xmin,xmax,exclude_min,exclude_max,n):
    sp = p.Spectrum(filename)
    sp.xarr.units = 'micron'
    sp.xarr.xtype = 'wavelength'
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars',color='grey')
    sp.baseline(xmin=xmin, xmax=xmax,exclude=[exclude_min, exclude_max],subtract=False,
                reset_selection=False,hightlight_fitregion=False,order=0)
    F = sp.data
    F = np.asarray(F)
    F_lambda = F.min()
    F_zero = sp.baseline.baselinepars
    line_depth = 1 - F_lambda/F_zero
    line_depth = np.float(line_depth) 
    print line_depth   
        
    sp2 = sp.copy()
    depth = []

    for w in range(n):
        sp2.data = sp.data + np.random.randn(sp.data.size)*sp.error
        sp2.baseline(xmin=xmin, xmax=xmax,exclude=[exclude_min, exclude_max],subtract=False,
                     reset_selection=False,hightlight_fitregion=False,order=0)
        F = sp2.data
        F = np.asarray(F)
        F_lambda = F.min()
        F_zero = sp2.baseline.baselinepars
        dist = 1 - F_lambda/F_zero
        dist = np.float(dist)
        depth.append(dist)

    depth = np.array(depth)

    plt.figure()
    mu,sigma = norm.fit(depth)
    print mu, sigma

    n,bins,patches = plt.hist(depth,10,normed=True,facecolor='lightblue')
    y = mlab.normpdf(bins,mu,sigma)
    plt.plot(bins,y,'r--',linewidth=2)
    plt.grid(True)
    plt.ylabel('Probability')
    plt.xlabel('Line Depth')
    plt.show()
    
if __name__=="__main__":
    xmin,xmax,exclude_min,exclude_max,n