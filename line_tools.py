#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from astropy import log

"""Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
Args:
----------
xmin,xmax - the specified interval in wavelength space  
excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
n - the number of Monte Carlo iterations 
Returns:
-------
- the mean and standard deviation of the equivalent width measured n times
- the spectrum plotted with the Voigt profile line fit (blue), the pseudo-continuum (yellow), and the approximated rectangle (green) 
- a histogram of the EqW distribution    
"""

def equivalent_width(filename,xmin,xmax,exclude_min,exclude_max,n):
    sp = p.Spectrum(filename)
    sp.xarr.units ='microns'                                                           
    sp.xarr.xtype = 'wavelength'
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', color='grey')                                         
    sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], 
                subtract=False, highlight_fitregion=False, 
                selectregion=True, order=0)    
    sp.specfit(plot=True, fittype='voigt', color='blue',
               guesses='moments', vheight=True)
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, components=False, 
                   annotate=True, loc='lower left', xmin=None, xmax=None)

    sp2 = sp.copy()
    EQWs = []              

    for w in range(n):
        sp2.data = sp.data + np.random.randn(sp.data.size)*sp.error
        sp2.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], 
                     subtract=False, highlight_fitregion=False, 
                     selectregion=True, order=0)           
        sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
        dist = sp2.specfit.EQW(plotcolor='g', fitted=False, components=False, 
                               annotate=True, loc='lower left', xmin=None, xmax=None)

        EQWs.append(dist)
    EQWs = np.array(EQWs)
    EQWs = EQWs*10000

    plt.figure()
    mu,sigma = norm.fit(EQWs) 
    print mu, sigma

    n,bins,patches = plt.hist(EQWs, 10, normed=True, facecolor='green', histtype='stepfilled')      
    y = mlab.normpdf(bins,mu,sigma)
    plt.plot(bins,y,'r--',linewidth=2)
    plt.grid(True)
    plt.ylabel('Probability')
    plt.xlabel('EQW')           
    plt.show()