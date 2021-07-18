#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import astropy.io.fits as pf
from os import path

"""
Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit.
"""

def equivalent_width(spec, xmin, xmax, exclude_min, exclude_max, n=1000, fldr=None, name=None, bandloc=6563, speclims=None):
    """"
    This is the main function to be invoked by the user. This function gathers the spectrum and metadata, sets up the output PDF figure, and calls the measure_equivalent_width() function to perform the calculation.

    Args:
    ----------
    spec - Numpy Array, List or String; if array or list, it should contain the spectrum with wavelength in Angstrom, flux, and optional flux errors as columns; if string, it should contain the filename with full path of the spectrum file
    xmin,xmax - Integers, the specified interval in wavelength space, which defines the region of interest
    excludemin, excludemax - Integers, the specified interval in wavelength space of the spectral feature, which binds the edges of the spectral feature itself
    n - Integer, the number of times the EqW measurement is repeated in the MCMC routine
    fldr - String, location where output figure is desired; if None, figure is saved in current folder
    name - String, if not None, it uses it to label the object
    bandloc - Float or Integer, the central location of the spectral line in Angstrom
    speclims - Numpy Array or List, the minimum and maximum wavelength values in Angstrom to be plotted; the fitting routine ignores this input

    Returns:
    -------
    - The mean and standard deviation of the distribution of n equivalent width measurements; if the MCMC routine is not run, then only the preliminary equivalent width measured is returned
    - A figure with three panels: one with the full spectrum, one with the fit to the spectral line, and one with a histogram from the MCMC results
    """

    # Define some constants ---------------------------------------------------
    global L_GRAY, GRAY, BLUE, GREEN, ORANGE, BLORDER
    L_GRAY = '#f2f2f2'
    GRAY = '#999999'
    BLUE = '#2166ac'
    GREEN = '#006633'
    ORANGE = '#FF9933'
    BLORDER = 1 # order for baseline fitting 

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
            # spec is neither a filename nor an array with spectrum
            if filename is None:
                print('ERROR: spec parameter has the wrong format.\n')
                return
        # Is the array an appropriate one?
        if len(specarr) not in (2,3):
            print('ERROR: spec parameter has the wrong format.\n')
            return
    
    # Are xmin, xmax, exclude_min, exclude_max reasonable?
    if xmin > xmax or exclude_min > exclude_max or xmin > exclude_min or xmax < exclude_max:
        print('ERROR: region of interest and/or edges of spectral line are incorrect.')
        return

    # Is fldr included?
    if fldr is None:
        fldr = './'
    else:
        if not path.exists(fldr):
            fldr = './'
    
    # Is name provided?
    if name is not None:
        objname = name
    else:
        objname = 'Object' # This could be something more meaningful...

    # Are speclims in the right format?
    if speclims is not None:
        if not isinstance(speclims, (list, np.ndarray)):
            print('ERROR: speclims parameter has the wrong format.\n')
            return
        elif len(speclims) != 2:
            print('ERROR: speclims parameter has the wrong format.\n')
            return

    # Attempt to load spectrum using input filename ---------------------------
    sp = None
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
    msk = sp.flux.data < 0
    sp.flux.mask = sp.flux.mask | msk

    # Set up figure -----------------------------------------------------------
    plt.close()
    plt.rc('font', size=7)
    plt.ion()
    fig = plt.figure(1, figsize=(3.30709*1.3,3.30709*1.6))
    ax1 = fig.add_subplot(311) # TOP PANEL To plot the full spectrum
    ax2 = fig.add_subplot(312) # MID PANEL To plot the fit
    ax3 = fig.add_subplot(313) # BOTTOM PANEL To plot the histogram
    plt.subplots_adjust(hspace=0.3, top=0.97, bottom=0.06, right=0.97, left=0.12)
    
    # Plot full spectrum (TOP PANEL) ------------------------------------------
    if speclims is None:
        tmpmin = sp.xarr.value[1]
        tmpmax = sp.xarr.value[-2]
    else:
        tmpwarning = False
        tmpmin = speclims[0]
        if tmpmin < sp.xarr.value[0]:
            tmpmin = sp.xarr.value[1]
            tmpwarning = True
        tmpmax = speclims[1]
        if tmpmax > sp.xarr.value[-1]:
            tmpmax = sp.xarr.value[-2]
            tmpwarning = True
        if tmpwarning:
            print('WARNING: Spectrum limits provided are beyond the range of the spectrum.\n')
    irange = np.where((sp.xarr.value >= tmpmin) & (sp.xarr.value <= tmpmax))[0]
    ax1.plot(sp.xarr.value[irange], sp.flux.data[irange], drawstyle='steps-mid', linewidth=0.8, color='k')
    # Beautify full spectrum plot
    ax1.set_xlim(xmin=tmpmin, xmax=tmpmax)
    ax1.set_title(objname, y=0.97)
    ax1.set_xlabel(r'Wavelength $(\AA)$', labelpad=0)
    ax1.set_ylabel('Flux', labelpad=-1)
    ax1.axvline(x=bandloc, linestyle=':', color=GRAY, linewidth=0.8)

    # Check the validity of xmin, xmax, exclude_min, exclude_max
    if xmin < sp.xarr.value[0] or xmax > sp.xarr.value[-1]:
        print('ERROR: region of interest defined by (xmin, xmax) is not possible.')
        return

    # Invoke pyspeckit to plot fit (MID PANEL) --------------------------------
    # (Note: the EW result here is only preliminary and is only provided to the user if they choose not to execute the next step )
    resultprel = measure_equivalent_width(sp, xmin, xmax, exclude_min, exclude_max, n, \
                                          bandloc, fig=fig, ax=ax2)
    if resultprel is None:
        return

    # Perform MCM routine to calculate EW and e_EW ----------------------------
    plt.show()
    tmpCommand = input('Continue with MCMC step? ([y]/n) ')
    if tmpCommand.upper() == 'N':
        # Eliminate the bottom panel
        ax3.set_visible(False)
        plt.savefig(fldr + objname + '_EWfit.pdf')
        return resultprel

    sp2 = sp.copy()
    EQWs = np.zeros(n)
    print('Begin MCMC iterations...')
    for w in range(n):
        if w % 100 == 0 or w == n - 1:
            print(w, end='  ', flush=True)
        # Determine how to add noise to spectrum
        if np.sum(sp.error) == 0:
            # If spectrum has no flux uncertainties, then get noise from noise in continuum
            icont = np.where(((sp.xarr.value >= xmin) & (sp.xarr.value < exclude_min)) |  \
                            ((sp.xarr.value > exclude_max) & (sp.xarr.value <= xmax)))[0]
            tmpcont = sp.data[icont]
            tmperr = np.std(tmpcont)
        else:
            # Otherwise, get noise from flux uncertainties
            tmperr = sp.error
        # Generate a new version of the flux array by adding random noise
        sp2.data = sp.data + np.random.randn(sp.data.size) * tmperr

        # Invoke pyspeckit to do equivalent width measurement
        mcmcresult = measure_equivalent_width(sp2, xmin, xmax, exclude_min, exclude_max, n, bandloc)
        EQWs[w] = mcmcresult

    # Calculate stats of MCMC results
    mu, sigma = norm.fit(EQWs)
    result = [mu, sigma]
    percentiles = [16, 50, 84]
    perc = np.percentile(EQWs, percentiles)
    
    # Plot histogram (BOTTOM PANEL)
    nums,bins,ptchs = ax3.hist(EQWs, 10, density=True, facecolor=GREEN, histtype='stepfilled')
    y = norm.pdf(bins, mu, sigma)
    ax3.plot(bins, y, color=ORANGE, linestyle='--', linewidth=2)    
    
    for ip,p_value in enumerate(perc):
        ax3.axvline(p_value, linestyle=':', color='k', linewidth=1.5)
        perctxt = ' ' + str(percentiles[ip]) + r'$^{th}$'
        # We want the x coords to be data coords, and the y coords to span 0-1 in axes coords:
        trans = ax3.get_xaxis_transform()
        ax3.text(p_value, 0.01, perctxt, fontsize=6, transform=trans)
    
    # Beautify histogram plot
    ax3.set_ylabel('PDF')
    ax3.set_xlabel(r'EW ($\AA$)', labelpad=0)
    resulttxt = 'EW = ' + format(mu,'.3f') + r' $\pm$ ' + format(sigma,'.3f')
    ax3.text(0.02,0.9, resulttxt, transform=ax3.transAxes, fontsize=7, \
             bbox=dict(facecolor=L_GRAY, ec='none', pad=0.3, boxstyle='round'))
    
    # Eliminate the preliminary EW result from mid panel
    ax2.get_children()[7].set_visible(False)

    # Save figure -------------------------------------------------------------
    plt.savefig(fldr + objname + '_EWfit.pdf')

    return result


def measure_equivalent_width(sp, xmin, xmax, exclude_min, exclude_max, n, bandloc, fig=None, ax=None):
    """Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
    
    Args:
    =======
    sp - pyspeckit instance of class Spectrum
    xmin,xmax - the specified interval of the spectrum to plot
    excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
    n - the number of Monte Carlo iterations
    bandloc - Integer or Float, wavelength location in Angstrom of the spectral feature
    fig - Matplotlib pyplot figure containing the axis where fit will be drawn
    ax - Matplotlib pyplot axis where the fit will be drawn; if None, function assumes that it is being invoked by MCMC step

    Returns:
    =======
    - the mean and standard deviation of the equivalent width measured n times
    - the spectrum plotted with the Voigt profile line fit (blue), the pseudo-continuum (yellow), 
    and the approximated rectangle (green) 
    - a histogram of the EqW distribution
    """
 

    # Get an estimate of flux error from noise in continuum
    # icont = np.where(((sp.xarr.value >= xmin) & (sp.xarr.value < exclude_min)) |  \
    #                 ((sp.xarr.value > exclude_max) & (sp.xarr.value <= xmax)))[0]
    # tmpcont = sp.flux[icont]
    # tmperr = np.std(tmpcont)
    
    # Set up plotter (when required) and fit baseline
    if ax is None:
        sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], \
                    subtract=False, reset_selection=False, order=BLORDER)
    else:
        sp.plotter(axis=ax, clear=False, xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', \
                   color=GRAY)
        sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], \
                    subtract=False, highlight_fitregion=False, reset_selection=False, \
                    selectregion=True, order=BLORDER)
        sp.baseline.annotate(loc='upper left', fontsize=6)

    # Fit Voigt profile to spectral feature
    if ax is None:
        sp.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
        ew = sp.specfit.EQW(fitted=True, xmin=None, xmax=None)
    else:
        sp.specfit(fittype='voigt', color=BLUE, guesses='moments')
        ew = sp.specfit.EQW(plot=True, plotcolor=GREEN, fitted=True, xmin=None, xmax=None)

        # Annotate results into figure
        sp.specfit.annotate(loc='lower left', fontsize=6, frameon=True)
        txt = 'EW = ' + format(ew,'.2f') + r' $\AA$'
        ax.text(0.995,0.02, txt, transform=ax.transAxes, fontsize=6, ha='right')

        # Beautify plot
        sp.plotter.axis.set_xlabel(r'Wavelength $(\AA)$', labelpad=0)
        sp.plotter.axis.set_ylabel('Flux', labelpad=-1)
        sp.plotter.axis.axvline(x=exclude_min, linestyle=':', color=GRAY, linewidth=0.8)
        sp.plotter.axis.axvline(x=exclude_max, linestyle=':', color=GRAY, linewidth=0.8)
        # We don't want the fit parameters box to block the bottom edge of exclude_* indicators:
        sp.specfit.fitleg.set_bbox_to_anchor((0,0.03), transform=ax.transAxes)
        lgd = ax.get_legend()
        plt.setp(lgd.get_frame(), color=L_GRAY)
        sp.specfit.fitleg.set_alpha(0.3)
        
        fig.canvas.draw() # This refreshes the figure so that everything shows
    
    return ew

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
    
    # Fix some parameters
    sp.specname = '' # This way, it is not printed as a title for the fit panel
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
                     