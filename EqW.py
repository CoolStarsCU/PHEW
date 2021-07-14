#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from os import path
import pdb


def equivalent_width(filename, xmin, xmax, exclude_min, exclude_max, n, fldr=None, name=None):
    """"
    This is the main function to be invoked by the user. This function gathers the spectrum and metadata, sets up the output PDF figure, and calls the measure_equivalent_width() function to perform the calculation.

    Args:
    ----------
    filename - String, spectrum file name with full path
    xmin,xmax - Integers, the specified interval in wavelength space, which defines the region of interest
    excludemin, excludemax - Integers, the specified interval in wavelength space of the spectral feature, which binds the edges of the spectral feature itself
    n - Integer, the number of times the EqW measurement is repeated in the MCMC step
    fldr - String, location where output figure is desired; if None, figure is saved in current folder
    name - String, if not None, it uses it to label the object

    Returns:
    -------
    - the mean and standard deviation of the equivalent width measured n times
    - A figure with a plot of the full spectrum; a plot of the spectral feature fit, with the Voigt profile line fit (blue), the pseudo-continuum (orange), and the approximated rectangle (green); and a histogram with the MCMC results of the EqW distribution
    """

    # Check that inputs are valid ---------------------------------------------
    # fldr
    if fldr is None:
        fldr = './'
    else:
        if not path.exists(fldr):
            fldr = './'

    # Invoke pyspeckit to perform equivalent width measurement ----------------
    mu, sigma = measure_equivalent_width(filename, xmin, xmax, exclude_min, exclude_max, n)

    # Save figure -------------------------------------------------------------
    fig = plt.gcf()
    plt.savefig(fldr + name + '_EWfit.pdf')

    return np.array([mu, sigma]) 


def measure_equivalent_width(filename, xmin, xmax, exclude_min, exclude_max, n):
    """Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
    
    Args:
    =======
    xmin,xmax - the specified interval of the spectrum to plot
    excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
    n - the number of Monte Carlo iterations 
    
    Returns:
    =======
    - the mean and standard deviation of the equivalent width measured n times
    - the spectrum plotted with the Voigt profile line fit (blue), the pseudo-continuum (yellow), 
    and the approximated rectangle (green) 
    - a histogram of the EqW distribution
    """
    sp = p.Spectrum(filename) # This has trouble loading flux error dimension
    sp.unit = 'Angstrom' 
    sp.xarr.xtype = 'wavelength'
    
    # Get an estimate of flux error from noise in continuum
    icont = np.where(((sp.xarr.value >= xmin) & (sp.xarr.value < exclude_min)) |  \
                    ((sp.xarr.value > exclude_max) & (sp.xarr.value <= xmax)))[0]
    tmpcont = sp.flux[icont]
    tmperr = np.std(tmpcont)
    
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', 
               color='grey')
    sp.baseline(xmin=xmin, xmax=xmax,  
                exclude=[exclude_min,exclude_max], 
                subtract=False, highlight_fitregion=False, 
                selectregion=True, order=0)    
    sp.specfit(plot=True, fittype='voigt', color='blue',
               guesses='moments', vheight=True)
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, components=False, 
                   annotate=True, loc='lower left', xmin=None, xmax=None)

    sp2 = sp.copy()
    EQWs = []              

    for w in range(n):
        sp2.data = sp.data + np.random.randn(sp.data.size)*tmperr
        sp2.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], 
                     subtract=False, highlight_fitregion=False, 
                     selectregion=True, order=0)
        sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
        dist = sp2.specfit.EQW(plotcolor='g', fitted=False, components=False, 
                               annotate=True, loc='lower left', xmin=None, xmax=None)

        EQWs.append(dist)
    EQWs = np.array(EQWs)
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, 
                   components=False, annotate=True, 
                   loc='lower left', xmin=None, xmax=None)

    # sp2 = sp.copy()
    # EQWs = np.zeros(n)

    # for w in range(n):
    #     sp2.data = sp.data + np.random.randn(sp.data.size)*sp.error
    #     sp2.baseline(xmin=xmin, xmax=xmax, 
    #                  exclude=[exclude_min,exclude_max],  
    #                  subtract=False, highlight_fitregion=False,
    #                  selectregion=True, order=0)
    #     sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
    #     dist = sp2.specfit.EQW(plotcolor='g', fitted=False, 
    #                            components=False, annotate=True,  
    #                            loc='lower left', xmin=None, xmax=None)


    fig = plt.figure()
    mu,sigma = norm.fit(EQWs) 
    print(mu, sigma)
    ax = fig.add_subplot(111)
    
    num,bins,patches = ax.hist(EQWs, 10, facecolor='green', density=True, histtype='stepfilled')      
    y = norm.pdf(bins,mu,sigma)
    ax.plot(bins,y,'r--',linewidth=2)
    ax.grid(True)
    ax.set_ylabel('Probability')
    ax.set_xlabel('EQW')           
    plt.show()
    num_bins = np.floor(np.log10(n)) * 10
    pdb.set_trace()
    # Light blue histogram for contrast and for R/G colorblind folks
    n,bins,patches = plt.hist(EQWs, int(num_bins), density=True,
                              facecolor='lightblue',  
                              histtype='stepfilled')      
    y = norm.pdf(bins,mu,sigma)
    ax.plot(bins,y,'r--',linewidth=2)

    # Plot the median and 68th percentile values for comparison
    perc = np.percentile(EQWs, [16, 50, 84])
    ax = plt.gca()
    for p_value in perc:
        ax.axvline(p_value,linestyle=":",color="k",lw=1.5)

    plt.ylabel('Probability')
    plt.xlabel('EQW')           
    plt.show()

    return mu, sigma

if __name__=="__main__":
    # xmin,xmax,exclude_min,exclude_max,n
    equivalent_width('2m1821_61_08jun11.txt',1.250, 1.255, 1.2514, 
                     1.2538, 1000)
                     
def final_plots(filename,xmin,xmax,exclude_min,exclude_max):
    vf = p.spectrum.models.inherited_voigtfitter.voigt_fitter()
    
    sp = p.Spectrum(filename)
    sp.xarr.units = 'micron'
    sp.xarr.xtype = 'wavelength'
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='fill',color='grey')  
    sp.baseline(xmin=xmin, xmax=xmax,exclude=[exclude_min,exclude_max],subtract=False,
                reset_selection=False,hightlight_fitregion=False,order=0)
    sp.specfit(plot=True, fittype='voigt', color='magenta', guesses='moments', 
               vheight=True)
    fwhm = sp.specfit.measure_approximate_fwhm(threshold='error', emission=False, 
                                               interpolate_factor=1024, plot=True, 
                                               grow_threshold=1, color='magenta')     
    sp.plotter.refresh()
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, continuum=0.5, components=False, annotate=True, loc='lower left', xmin=None,
    xmax=None)
    sp.plotter.refresh()
    xarr_fit_units = 'microns'
    plt.ylabel('Normalized Flux')
    plt.xlabel('Wavelength ($\mu m$)')