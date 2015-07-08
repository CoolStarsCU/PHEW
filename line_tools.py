#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm


"""Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam

PARAMETERS:
----------
xmin,xmax - the specified interval in wavelength space  
excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature 
n - the number of times the EqW measurement is repeated to 


RETURNS:
-------
- the mean and standard deviation of the equivalent width measured n times
- the spectrum plotted with the Voigt profile line fit (blue), the pseudo-continuum (yellow), and the approximated rectangle (green) 
- a histogram of the EqW distribution    
"""

# high-res NIRSPEC 
D = {'2M1821':'/Users/munazzaalam/Desktop/blueseq/2m1821_61_08jun11.txt'}   

# optical, low-res SALT 
# D = {'RX1132':'/Users/munazzaalam/saltspec/RX1132-2651A_000_2014-02-27.txt'}      

def equivalent_width(xmin,xmax,exclude_min,exclude_max,n):
    for obj in D.keys():
        obj = unicode(obj)

        # for data in D[obj].keys():
        #     data = unicode(data)

            # sp = p.Spectrum(D[obj][data])
        sp = p.Spectrum(D[obj])
        sp.xarr.units ='microns'                                                           
        sp.xarr.xtype = 'wavelength'
        sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars', color='grey')     
        sp.specfit(fittype='voigt', color='blue')                                     
        sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], subtract=False, reset_selection=True, highlight_fitregion=False, order=0)    
        sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, continuum=0.5, components=False, annotate=True, loc='lower left', xmin=None, xmax=None)

        sp2 = sp.copy()
        EQWs = []              

        for w in range(n):
            sp2.data = sp.data + np.random.randn(sp.data.size)*sp.error
            sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
            sp2.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], subtract=False, reset_selection=True, highlight_fitregion=False, order=0)
            dist = sp2.specfit.EQW(plotcolor='g', fitted=False, continuum=0.5, components=False, annotate=True, loc='lower left', xmin=None, xmax=None)
            EQWs.append(dist)
        EQWs = np.array(EQWs)
        EQWs = EQWs*10000

        plt.figure()
        mu,sigma = norm.fit(EQWs)
        print obj, mu, sigma

        n,bins,patches = plt.hist(EQWs,10,normed=True,facecolor='green',histtype='stepfilled')      
        y = mlab.normpdf(bins,mu,sigma)
        plt.plot(bins,y,'r--',linewidth=2)
        plt.grid(True)
        plt.ylabel('Probability')
        plt.xlabel('EQW ($\AA$)')           
        plt.show()