#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt

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