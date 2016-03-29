#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from astropy import log

"""
Calculate the full width at half maximum (FWHM) of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam

Args:
----------
xmin,xmax - the specified interval of the spectrum to plot 
excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
n - the number of Monte Carlo iterations 

Returns:
-------
- the mean and standard deviation of the FWHM measured n times
- the spectrum plotted with the the Voigt profile line fit (blue), the pseudo-continuum (yellow), 
and the FWHM (blue)
- a histogram of the FWHM distribution 
"""

def measure_fwhm(filename,xmin,xmax,exclude_min,exclude_max,n):
    vf = p.spectrum.models.inherited_voigtfitter.voigt_fitter()
    
    #read & plot in spectrum 
    sp = p.Spectrum(filename)
    sp.xarr.units = 'micron'
    sp.xarr.xtype = 'wavelength'
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars',color='grey')  
    
    #fit pseudocontinuum
    sp.baseline(xmin=xmin, xmax=xmax,exclude=[exclude_min, exclude_max],subtract=False,
                reset_selection=False,hightlight_fitregion=False,order=0)

    #fit voigt profile
    sp.specfit(plot=True, fittype='voigt', color='blue', guesses='moments', 
               vheight=True)

    #measure FWHM
    fwhm = sp.specfit.measure_approximate_fwhm(threshold='error', emission=False, 
                                               interpolate_factor=1024, plot=True, 
                                               grow_threshold=1)
    sp.plotter.refresh()
    fwhm = "{0:0.06f}".format(fwhm)
    fwhm = fwhm[:-7]
    print fwhm
    xarr_fit_units = 'microns'
    plt.ylabel('Normalized Flux')
    plt.xlabel('Wavelength ($\mu m$)')
    
 
    #copy flux array for Monte Carlo
    sp.specfit(guesses=sp.specfit.parinfo.values)
    sp2 = sp.copy()
    FWHM = []
    
    #Monte Carlo iterations to estimate uncertainty
    for w in range(n):
        sp2.data = sp.data + np.random.randn(sp.data.size)*sp.error
        sp2.baseline(xmin=xmin, xmax=xmax,exclude=[exclude_min, exclude_max], subtract=False,
                     reset_selection=False, highlight_fitregion=False, order=0)
        sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
        dist = sp2.specfit.measure_approximate_fwhm(threshold='error', emission=False,interpolate_factor=1024,grow_threshold=1)
        dist = "{0:0.07f}".format(dist)
        dist = dist[:-7]
        dist = float(dist)
        print dist
        FWHM.append(dist)
        
    FWHM = np.array(FWHM)
    FWHM = FWHM*10000
    
    #plot FWHM probability distribution
    plt.figure()
    mu,sigma = norm.fit(FWHM)
    print mu, sigma

    n,bins,patches = plt.hist(FWHM,10,normed=True,facecolor='lightblue')
    y = mlab.normpdf(bins,mu,sigma)
    plt.plot(bins,y,'r--',linewidth=2)
    plt.grid(True)
    plt.ylabel('Probability')
    plt.xlabel('FWHM ($\AA$)')
    plt.show()
    
if __name__=="__main__":
    xmin,xmax,exclude_min,exclude_max,n