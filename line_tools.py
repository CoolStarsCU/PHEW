#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from astropy import log

<<<<<<< HEAD
"""Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
=======
"""
Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
>>>>>>> a2e51d542e166aa0c2b16458ede8699373d797b4
Args:
----------
xmin,xmax - the specified interval of the spectrum to plot
excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
n - the number of Monte Carlo iterations 
<<<<<<< HEAD
=======

>>>>>>> a2e51d542e166aa0c2b16458ede8699373d797b4
Returns:
-------
- the mean and standard deviation of the equivalent width measured n times
- the spectrum plotted with the Voigt profile line fit (blue), the pseudo-continuum (yellow), 
  and the approximated rectangle (green) 
- a histogram of the EqW distribution
"""

<<<<<<< HEAD
def equivalent_width(filename,xmin,xmax,exclude_min,exclude_max,n):
    sp = p.Spectrum(filename)
    sp.xarr.units ='microns'                                                           
    sp.xarr.xtype = 'wavelength'
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', color='grey')                                         
    sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], 
=======
def measure_equivalent_width(filename,xmin,xmax,exclude_min,exclude_max,n):
    sp = p.Spectrum(filename)
    sp.xarr.units = 'micron' 
    sp.xarr.xtype = 'wavelength'
    print sp.xarr
    print sp.data

    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', 
               color='grey')
    sp.baseline(xmin=xmin, xmax=xmax,  
                exclude=[exclude_min,exclude_max], 
>>>>>>> a2e51d542e166aa0c2b16458ede8699373d797b4
                subtract=False, highlight_fitregion=False, 
                selectregion=True, order=0)    
    sp.specfit(plot=True, fittype='voigt', color='blue',
               guesses='moments', vheight=True)
<<<<<<< HEAD
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
=======
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, 
                   components=False, annotate=True, 
                   loc='lower left', xmin=None, xmax=None)

    sp2 = sp.copy()
    EQWs = np.zeros(n)

    for w in range(n):
        sp2.data = sp.data + np.random.randn(sp.data.size)*sp.error
        sp2.baseline(xmin=xmin, xmax=xmax, 
                     exclude=[exclude_min,exclude_max],  
                     subtract=False, highlight_fitregion=False, 
                     selectregion=True, order=0)
        sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
        dist = sp2.specfit.EQW(plotcolor='g', fitted=False, 
                               components=False, annotate=True,  
                               loc='lower left', xmin=None, xmax=None)

        # Convert from microns to angstroms
        EQWs[w] = dist*10000.0
>>>>>>> a2e51d542e166aa0c2b16458ede8699373d797b4

    plt.figure()
    mu,sigma = norm.fit(EQWs) 
    print mu, sigma

<<<<<<< HEAD
    n,bins,patches = plt.hist(EQWs, 10, normed=True, facecolor='green', histtype='stepfilled')      
    y = mlab.normpdf(bins,mu,sigma)
    plt.plot(bins,y,'r--',linewidth=2)
    plt.grid(True)
    plt.ylabel('Probability')
    plt.xlabel('EQW')           
    plt.show()
=======
    num_bins = np.floor(np.log10(n)) * 10

    # Light blue histogram for contrast and for R/G colorblind folks
    n,bins,patches = plt.hist(EQWs, num_bins, normed=True, 
                              facecolor='lightblue',  
                              histtype='stepfilled')      
    y = mlab.normpdf(bins,mu,sigma)
    plt.plot(bins,y,'r--',linewidth=2)

    # Plot the median and 68th percentile values for comparison
    perc = np.percentile(EQWs, [16, 50, 84])
    ax = plt.gca()
    for p_value in perc:
        ax.axvline(p_value,linestyle=":",color="k",lw=1.5)

    plt.ylabel('Probability')
    plt.xlabel('EQW')           
    plt.show()

if __name__=="__main__":
    # xmin,xmax,exclude_min,exclude_max,n
    equivalent_width('2m1821_61_08jun11.txt',1.250, 1.255, 1.2514, 
                     1.2538, 100)
>>>>>>> a2e51d542e166aa0c2b16458ede8699373d797b4
