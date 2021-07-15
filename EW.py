#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
import astropy.io.fits as pf
from os import path
import pdb


def equivalent_width(spec, xmin, xmax, exclude_min, exclude_max, n, fldr=None, name=None):
    """"
    This is the main function to be invoked by the user. This function gathers the spectrum and metadata, sets up the output PDF figure, and calls the measure_equivalent_width() function to perform the calculation.

    Args:
    ----------
    spec - Numpy Array, List or String; if array or list, it should contain the spectrum with wavelength in Angstrom, flux, and optional flux errors as columns; if string, it should contain the filename with full path of the spectrum file
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
    # Is spec a filename?
    if isinstance(spec, str):
        filename = spec
    else:
        filename = None
        # Is spec an array?
        specarr = None
        if isinstance(spec, (list, np.ndarray)):
            specarr = list(spec)
        else:
            # spec is neither a filename nor a list with spectrum
            if filename is None:
                print('ERROR: spec parameter has the wrong format.\n')
                return
        # Is the array an appropriate one?
        if len(specarr) not in (2,3):
            print('ERROR: spec parameter has the wrong format.\n')
            return
    
    # Is fldr included?
    if fldr is None:
        fldr = './'
    else:
        if not path.exists(fldr):
            fldr = './'

    sp = None

    # Attempt to load spectrum using input filename ---------------------------
    if filename is not None:
        sp = readspec(filename)
        if sp == 1:
            return
    
    # Get spectrum data from input array --------------------------------------
    # (Note: the empty header pf.Header() is added to avoid a pyspeckit warning message)
    if sp is None:
        if len(specarr) < 3:
            # Only wavelength and flux
            sp = p.Spectrum(data=specarr[1], xarr=specarr[0], \
                            xarrkwargs={'unit':'Angstrom'}, \
                            maskdata=True, header=pf.Header()) 
        else:
            if np.isnan(specarr[2]).all():
                # Only wavelength and flux
                sp = p.Spectrum(data=specarr[1], xarr=specarr[0], \
                                xarrkwargs={'unit':'Angstrom'}, \
                                maskdata=True, header=pf.Header())                 
            else:
                # Wavelength, flux, and flux errors
                sp = p.Spectrum(data=specarr[1], xarr=specarr[0], error=specarr[2], \
                                xarrkwargs={'unit':'Angstrom'}, \
                                maskdata=True, header=pf.Header())
        sp.xarr.xtype = 'wavelength'
    
    # Mask negative flux values
    msk = sp.data.data < 0
    sp.data.mask = sp.data.mask | msk


    # Invoke pyspeckit to perform equivalent width measurement ----------------
    result = measure_equivalent_width(sp, xmin, xmax, exclude_min, exclude_max, n)
    if result is None:
        return

    # Save figure -------------------------------------------------------------
    fig = plt.gcf()
    plt.savefig(fldr + name + '_EWfit.pdf')

    return result 


def measure_equivalent_width(sp, xmin, xmax, exclude_min, exclude_max, n):
    """Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
    
    Args:
    =======
    sp - pyspeckit instance of class Spectrum
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

    fig = plt.figure()
    mu,sigma = norm.fit(EQWs) 
    print(mu, sigma)
    ax = fig.add_subplot(111)
    
    num,bins,patches = ax.hist(EQWs, 10, facecolor='green', density=True, histtype='stepfilled')      
    y = norm.pdf(bins,mu,sigma)
    ax.plot(bins,y,'r--',linewidth=2)
    ax.grid(True)
    ax.set_ylabel('Probability')
    ax.set_xlabel('EW')           
    plt.show()
    num_bins = np.floor(np.log10(n)) * 10
    
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
    plt.xlabel('EW')           
    plt.show()

    return [mu, sigma]

def readspec(fname):
    """ Uses pyspeckit to attempt to read fits file of spectrum and load it into an instance of class Spectrum."""

    try:
        sp = p.Spectrum(fname, maskdata=True)
    except FileNotFoundError:
        print('ERROR: Spectrum file not found.\n')
        return 1
    except TypeError:
        print('ERROR: Spectrum file could not be loaded.\n')
        return 1
    
    # Fix wavelength units
    sp.xarr.xtype = 'wavelength' 
    if sp.xarr.unit != 'Angstrom':
        sp.xarr.convert_to_unit('Angstrom')

    # Warn the user if no flux errors present; p.Spectrum has trouble loading flux errors sometimes
    ierrs = np.where(sp.error.data != 0)[0]
    if len(ierrs) == 0:
        print('WARNING: No flux errors were loaded.\n')

    return sp

 
if __name__=="__main__":
    # xmin,xmax,exclude_min,exclude_max,n
    equivalent_width('2m1821_61_08jun11.txt',1.250, 1.255, 1.2514, 
                     1.2538, 1000)
                     