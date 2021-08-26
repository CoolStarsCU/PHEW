#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import astropy.io.fits as pf
from os import path
import pdb

"""
Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit.
"""

def equivalent_width(spec, bandloc, xmin, xmax, exclude_min, exclude_max, mc=True, n=1000, fldr=None, name=None, speclims=None, outname=None, blorder=1, interactive=True, clobber=True):
    """"
    This is the main function to be invoked by the user. This function gathers the spectrum and metadata, sets up the output PDF figure, and calls the __measure_equivalent_width() function to perform the calculation.

    Args:
    ----------
    spec - Numpy Array, List or String; if array or list, it should contain the spectrum with wavelength in Angstrom, flux, and optional flux errors as columns; if string, it should contain the filename with full path of the spectrum file
    bandloc - Float or Integer, the central location of the spectral line in Angstrom
    xmin,xmax - Integers, the specified interval in wavelength space, which defines the region of interest
    excludemin, excludemax - Integers, the specified interval in wavelength space of the spectral feature, which binds the edges of the spectral feature itself; the pseudo-continuum will be defined using the spectral ranges (xmin,exclude_min) and (exclude_max, xmax)
    mc - Boolean, whether to perform the MC iteration to estimate EW uncertainty; this parameter is relevant only in non-interactive mode (interactive=False)
    n - Integer, the number of times the EW measurement is repeated in the MC iteration
    fldr - String, location where output figure is desired; if None, figure is saved in current folder
    name - String, if not None, it uses it to label the object
    speclims - Numpy Array or List, the minimum and maximum wavelength values in Angstrom to be plotted; the fitting routine ignores this input
    outname - String, specifies the output figure filename, if user wants something different from the default name convention used here (i.e., name + _EWfit)
    blorder - Integer, the order of the polynomial for fitting the pseudo-continuum
    interactive - Boolean, whether to ask the user to confirm performing the MC iteration; if True, the mc parameter is overridden by the user interactively
    clobber - Boolean, whether to overwrite existing figure file

    Returns:
    -------
    - A list with the EW and its error e_EW, adopted from the mean and standard deviation of the distribution of n equivalent width measurements; if the MC iteration is not run, then only the preliminary equivalent width measured is returned as a float
    - A figure with three panels: one with the full spectrum, one with the fit to the spectral line, and if the MC routine is run, one additional panel with a histogram from the MC results
    """

    # Define some constants ---------------------------------------------------
    global L_GRAY, GRAY, BLUE, GREEN, ORANGE
    L_GRAY = '#f2f2f2'
    GRAY = '#999999'
    BLUE = '#2166ac'
    GREEN = '#006633'
    ORANGE = '#FF9933'

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
                print('ERROR: spec parameter has the wrong format.')
                return
        # Is the array an appropriate one?
        # (there is probably a more elegant way to check the acceptable array shape and size)
        if len(specarr) not in (2,3):
            print('ERROR: spec parameter has the wrong format.')
            return
        if len(specarr[0]) < 10:
            print('ERROR: spec parameter has the wrong format.')
            return
    
    # Are xmin, xmax, exclude_min, exclude_max reasonable?
    if xmin > xmax or exclude_min > exclude_max or xmin > exclude_min or xmax < exclude_max:
        print('ERROR: region of interest and/or edges of spectral line are not appropriate.')
        return

    # is bandloc reasonable?
    if isinstance(bandloc, (list, np.ndarray, str)):
        print('ERROR: bandloc parameter has the wrong format.')
        return
    if bandloc < xmin or bandloc > xmax:
        print('ERROR: bandloc parameter is outside the region of interest (xmin, xmax).')
        return
    
    # Is mc Boolean?
    if not isinstance(mc, bool):
        print('ERROR: mc parameter must be a boolean.')
        return
    
    # Is n a number?
    if not isinstance(n, (int, float)):
        print('ERROR: n parameter must be an integer.')
        return
    n = int(n)
    # We are NOT checking whether n is too small.

    # Is fldr included?
    if fldr is None:
        fldr = './'
    else:
        if fldr[-1] !='/':
            fldr = fldr + '/'
        if not path.exists(fldr):
            print('WARNING: fldr parameter invalid. Output figure will be saved in currenty directory.')
            fldr = './'
    
    # Is name provided?
    if name is not None:
        objname = name
    else:
        objname = 'Object' # This could be something more meaningful...

    # Are speclims in the right format?
    if speclims is not None:
        if not isinstance(speclims, (list, np.ndarray)):
            print('ERROR: speclims parameter has the wrong format.')
            return
        elif len(speclims) != 2:
            print('ERROR: speclims parameter has the wrong format.')
            return

    # Is outname reasonable?
    if outname is not None and not isinstance(outname, str):
        print('ERROR: outname parameter should be a string.')
        return
    if outname is None:
        fname = objname
    else:
        fname = outname

    # Is blorder reasonable?
    if not isinstance(blorder, (int, float)):
        print('ERROR: blorder parameter should be an integer number.')
        return
    blorder = int(blorder)

    # Is interactive Boolean?
    if not isinstance(interactive, bool):
        print('ERROR: interactive parameter must be a boolean.')
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

    # Normalize spectrum and determine relevant ranges ------------------------
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
            print('WARNING: Spectrum limits provided are beyond the range of the spectrum.')
    # Define the range to be plotted in the top panel
    irange = np.where((sp.xarr.value >= tmpmin) & (sp.xarr.value <= tmpmax))[0]
    
    # Determine if spectral line to be measured is within the wavelength range of the spectrum
    if xmin < sp.xarr.value[0] or xmax > sp.xarr.value[-1]:
        print('ERROR: spectral line to be measured is beyond the range of the spectrum.')
        return

    # Normalize spectrum by the flux at the spectral line location
    inorm = np.where(sp.xarr.value >= bandloc)[0][0]
    normval = sp.flux.data[inorm]
    fluxnorm = np.ma.core.MaskedArray(data=sp.flux.data / normval * 10,
                                      mask=sp.flux.mask, fill_value=sp.flux.fill_value)
    sp.flux = fluxnorm
    if np.sum(sp.error) != 0:
        efluxnorm = np.ma.core.MaskedArray(data=sp.error.data / normval * 10,
                                           mask=sp.error.mask, fill_value=sp.error.fill_value)  
        sp.error = efluxnorm

    # Set up figure -----------------------------------------------------------
    plt.close()
    plt.rc('font', size=7)
    plt.ion()
    fig = plt.figure(1, figsize=(3.30709*1.3,3.30709*1.6))
    ax1 = fig.add_subplot(311) # TOP PANEL To plot the full spectrum
    ax2 = fig.add_subplot(312) # MID PANEL To plot the fit
    ax3 = fig.add_subplot(313) # BOTTOM PANEL To plot the histogram
    plt.subplots_adjust(hspace=0.3, top=0.97, bottom=0.06, right=0.97, left=0.1)
    ax3.set_visible(False)

    # Plot full spectrum (TOP PANEL) ------------------------------------------
    ax1.plot(sp.xarr.value[irange], sp.flux.data[irange], drawstyle='steps-mid', \
             linewidth=0.8, color='k')
    ax1.set_xlim(xmin=tmpmin, xmax=tmpmax)
    ax1.set_title(objname, y=0.97)
    ax1.set_xlabel(r'Wavelength $(\AA)$', labelpad=0)
    ax1.set_ylabel('Normalized Flux', labelpad=-1)
    ax1.axvline(x=bandloc, linestyle=':', color=GRAY, linewidth=0.8)

    # Check the validity of xmin, xmax, exclude_min, exclude_max
    if xmin < sp.xarr.value[0] or xmax > sp.xarr.value[-1]:
        print('ERROR: region of interest defined by (xmin, xmax) is not possible.')
        return

    # Invoke pyspeckit to plot fit (MID PANEL) --------------------------------
    # (Note: the EW result here is only preliminary and is only provided to the user if they choose not to execute the next step)
    resultprel, guesses = __measure_equivalent_width(sp, xmin, xmax, exclude_min, exclude_max, \
                                                   blorder, ax=ax2)
    if resultprel is None:
        print('ERROR: Pyspeckit could not fit the spectral feature.')
        return

    # Determine whether the MC iteration should be run ----------------------
    if interactive:
        # In interative mode, ask the user if the MC iteration should be run
        fig.canvas.draw() # This refreshes the figure so that everything drawn so far shows
        plt.show()
        tmpCommand = None
        while tmpCommand not in ('y', 'n', ''):
            tmpCommand = input('Continue with MC step? ([y]/n) ')
            if tmpCommand == 'n':
                fig.text(0.5, 0.18, '(MC iteration skipped)', ha='center', fontstyle='italic')
                __savefig(fldr, fname, clobber)
                plt.close()
                return resultprel
            elif tmpCommand == '' or tmpCommand == 'y':
                mc = True
    else:
        # In non-interactive mode, use the mc parameter to determine whether to run the MC iteration 
        if not mc:
            fig.text(0.5, 0.18, '(mc iteration skipped)', ha='center', fontstyle='italic')
            __savefig(fldr, fname, clobber)
            plt.close()
            return resultprel

    # Perform MC iteration to calculate EW and e_EW -------------------------
    if mc:
        sp2 = sp.copy()
        EQWs = np.zeros(n)
        print('Begin MC iterations...')
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

            mcresult = __measure_equivalent_width(sp2, xmin, xmax, exclude_min, exclude_max, \
                                                  blorder, guesses)
            EQWs[w] = mcresult

        # Calculate stats of MC results
        mu, sigma = norm.fit(EQWs)
        result = [mu, sigma]
        percentiles = [16, 50, 84]
        perc = np.percentile(EQWs, percentiles)
        
        # Plot histogram (BOTTOM PANEL)
        ax3.set_visible(True)
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
        __savefig(fldr, fname, clobber)
        plt.close()

        return result


def __measure_equivalent_width(sp, xmin, xmax, exclude_min, exclude_max, blorder, guesses=None, ax=None):
    """Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit.
    
    Args:
    ----------
    sp - pyspeckit instance of class Spectrum
    xmin,xmax - the specified interval of the spectrum to plot
    excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature
    blorder - Integer, the order of the polynomial used to fit the pseudo-continuum
    guesses - List, the parameter guesses for the Voigt profile to be fitted; it is only used when the function is invoked by the MC iteration
    ax - Matplotlib pyplot axis where the fit will be drawn; if None, it is assumed that function is being invoked by the MC iteration

    Returns:
    ----------
    - The equivalent width measurement as a float
    """
    
    # Set up plotter (when required) and fit baseline
    if ax is None:
        sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], \
                    subtract=False, reset_selection=False, order=blorder)
    else:
        sp.plotter(axis=ax, clear=False, xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', \
                   color=GRAY)
        sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], \
                    subtract=False, highlight_fitregion=False, reset_selection=False, \
                    selectregion=True, order=blorder)
        sp.baseline.annotate(loc='upper left', fontsize=6)

    # Fit Voigt profile to spectral feature
    if ax is None:
        # Do the fit, but don't plot anything
        # Use the parameter values from the preliminary fit to guess the parameters here
        # (Ignore a numpy divide-by-zero error coming from within Pyspeckit)
        with np.errstate(divide='ignore'):
            sp.specfit(fittype='voigt', guesses=guesses)
            ew = sp.specfit.EQW(fitted=True, xmin=None, xmax=None)

            return ew
    else:
        # Do the fit and plot it
        # Use the moments to guess the parameters
        # (Ignore a numpy divide-by-zero error coming from within Pyspeckit)
        with np.errstate(divide='ignore'):
            sp.specfit(fittype='voigt', color=BLUE, guesses='moments')
            ew = sp.specfit.EQW(plot=True, plotcolor=GREEN, fitted=True, xmin=None, xmax=None)

        # Annotate results into figure
        sp.specfit.annotate(loc='upper right', fontsize=6, frameon=True)
        txt = 'EW = ' + format(ew,'.2f') + r' $\AA$'
        ax.text(0.995, 0.02, txt, transform=ax.transAxes, fontsize=6, ha='right')

        # Beautify plot
        sp.plotter.axis.set_xlabel(r'Wavelength $(\AA)$', labelpad=0)
        sp.plotter.axis.set_ylabel('Normalized Flux', labelpad=-1)
        sp.plotter.axis.axvline(x=exclude_min, linestyle=':', color=GRAY, linewidth=0.8)
        sp.plotter.axis.axvline(x=exclude_max, linestyle=':', color=GRAY, linewidth=0.8)
        # Put box around baseline fit description
        ax.artists[0].set_frame_on(True)
        plt.setp(ax.artists[0].get_frame(), color=L_GRAY)
        # Put box around spectral fit description
        plt.setp(ax.artists[1].get_frame(), color=L_GRAY)
        ax.artists[1].set_alpha = 0.1        
        sp.specfit.fitleg.set_bbox_to_anchor((1,1), transform=ax.transAxes) # This line does nothing if it's (1,1). We leave it here in case we want to move the box around...
        
        # In addition to the EW result, return the Voigt fit parameters from this preliminary fit, so that they can be used as guesses in the MC iteration
        return ew, sp.specfit.parinfo.values

def readspec(fname):
    """ Uses pyspeckit to attempt to read fits file of spectrum and load it into an instance of class Spectrum."""

    try:
        sp = p.Spectrum(fname, maskdata=True)
    except FileNotFoundError:
        print('ERROR: Spectrum file not found.')
        return 1
    except TypeError:
        print('ERROR: Spectrum file could not be loaded.')
        return 1
    
    # Fix some parameters
    sp.specname = '' # This way, it is not printed as a title for the fit panel
    sp.xarr.xtype = 'wavelength' 
    if sp.xarr.unit != 'Angstrom':
        sp.xarr.convert_to_unit('Angstrom')

    # Warn the user if no flux errors present; p.Spectrum has trouble loading flux errors sometimes
    ierrs = np.where(sp.error.data != 0)[0]
    if len(ierrs) == 0:
        print('WARNING: No flux errors were loaded.')

    return sp


def __savefig(fldr, fname, clobber):
    """ Saves output figure. Determines output filename to avoid over-writing, if requested."""

    if clobber:
        plt.savefig(fldr + fname + '_EWfit.pdf')
    else:
        exists = True
        cntnum = 0
        while exists:
            if cntnum == 0:
                cnt = ''
            else:
                cnt = '_' + str(cntnum)
            exists = path.exists(fldr + fname + '_EWfit' + cnt + '.pdf')
            if exists:
                cntnum = cntnum + 1
        plt.savefig(fldr + fname + '_EWfit' + cnt + '.pdf')
    return
                     