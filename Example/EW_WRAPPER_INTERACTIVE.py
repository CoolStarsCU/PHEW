#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 11:37:13 2021
@author: stanislavdelaurentiis

Instructions:
1. from EW_WRAPPER_INTERACTIVE import *
2. Run the function fit_runner(), with the absolute path to the directory with all fits file spectra as input. With fit_runner(), user will decide which spectra have a measurable spectral line, and which wavelength parameters to apply to each spectrum. The output will be in a csv file called "EW_FITS.csv".
3. Run the function mc_run(), with the absolute path to the directory with all fits file spectra as input. The function mc_run() uses the parameters defined in the previous step to run the MC routine in the EW.equivalent_width() function to calculate the final EW and uncertainty measurements for each spectrum. The output will be in a csv file called "EW_VALS_ALL.csv".

"""

import readspec as rs
import csv
import os
import subprocess
import shutil
import EW as ew


def fit_runner(path):
    '''
    (Created By Stanislav DeLaurentiis, 2021)
    
    Will run through all the fits files in the input directory and
    apply EW.equivalent_width() (with mc=False and interactive=False) to each file, ensuring that
    the curve is fit well to the spectral feature.
    Will then organize each fits file and its output pdf figure into a directory
    correlating with its parameters, while recording specname, parameters, and 
    type in EW_FITS.csv. 
    
    Parameters
    ----------
    path : string
        The absolute path of the directory that contains the desired fits files.

    Returns
    -------
    None.
    
    Outputs
    --------
    EW_FITS.csv

    '''
    
    f=open(path+'/EW_FITS.csv', 'w')
    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
    writer.writeheader()
    f.close()
    genparamcounter=0
    amount_in=0
    for spec in os.listdir(path):
        if spec.endswith('.fits'):
                amount_in=amount_in+1
    if amount_in==0:
        print('ERROR: There are no .fits files in the directory to fit.')
        return(None)
    determined=[]
    while amount_in>0:
        genparamcounter=genparamcounter+1
        print('There are '+str(amount_in)+ ' .fits files that need to be fit')
        print('Please input the general paramaters you wish to apply to these spectra')
        bandloc=int(input('Band Loc: '))
        xmin=int(input('XMIN: '))
        xmax=int(input('XMAX: '))
        excludemin=int(input('EXCLUDE MIN: '))
        excludemax=int(input('EXCLUDE MAX: '))
        c=0
        for spec in os.listdir(path):
            if str(spec).endswith('.fits'):
                if spec in determined:
                    continue
                c=c+1
                print('\n'+spec)
                print(c)
                print('\n')
                tmploc = spec.find('.fit')
                specname= spec[:tmploc]
                try:
                    specinfo=rs.read_spec(path+"/"+spec) #reading the fits file into array
                    ew.equivalent_width(specinfo, bandloc, xmin, xmax, excludemin, excludemax, name=specname, fldr=path, mc=False, interactive=False)
                except(AttributeError, ValueError):
                    print('ERROR: Unable to read this spectrum.')
                    f=open(path+'/EW_FITS.csv', 'a')
                    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                    writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'UNABLE TO READ', 'COMMENTS':''})
                    f.close()
                    amount_in=amount_in-1
                    determined.append(spec)
                    continue
                subprocess.call(['open', path+'/'+specname+'_EWfit.pdf'])
                print('Save this fit using these parameters?')
                first_yn=input('y or (n): ')
                if first_yn=='':
                    first_yn='n'
                if first_yn=='y':
                    print('File Saved.')
                    print('PARAMS: '+str([xmin, xmax, excludemin, excludemax]))
                    comment=input('Any comments? [Press enter if none]\n')
                    f=open(path+'/EW_FITS.csv', 'a')
                    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                    writer.writerow({'SPEC':specname, 'BANDLOC':bandloc, 'XMIN':xmin, 'XMAX':xmax, 'EXCLUDEMIN':excludemin, 'EXCLUDEMAX':excludemax, 'TYPE':'FITTED', 'COMMENTS':comment})
                    f.close()
                    determined.append(spec)
                    amount_in=amount_in-1
                    continue
                print('Re-apply parameters for this individual spectrum?')
                smallparam_yn=input('y or (n): ')
                if smallparam_yn=='':
                    smallparam_yn='n'
                if smallparam_yn=='n':
                    print('Is the spectrum flat and EW effectively zero?')
                    zero_yn=input('y or (n): ')
                    if zero_yn=='':
                        zero_yn='n'
                    if zero_yn=='y':
                        print('This spectrum will now be filed as having ZERO EW')
                        comment=input('Any comments? [Press enter if none]\n')
                        f=open(path+'/EW_FITS.csv', 'a')
                        writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                        writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'ZERO', 'COMMENTS':comment})
                        f.close()
                        determined.append(spec)
                        amount_in=amount_in-1
                        continue
                    if zero_yn=='n':
                        print('Attempt to fit this spectrum again?')
                        hope_yn=input('y or (n): ')
                        if hope_yn=='':
                            hope_yn='n'
                        if hope_yn=='n':
                            print('This spectrum will be filed as undetermined')
                            comment=input('Any comments? [Press enter if none]\n')
                            f=open(path+'/EW_FITS.csv', 'a')
                            writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                            writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'UNDETERMINED', 'COMMENTS':comment})
                            f.close()
                            determined.append(spec)
                            amount_in=amount_in-1
                            continue
                        if hope_yn=='y':
                            print('The spectrum will remain in the pool for further review')
                            continue
                while smallparam_yn=='y':
                    print('Enter your choice of unique parameters')
                    small_xmin=int(input('XMIN: '))
                    small_xmax=int(input('XMAX: '))
                    small_excludemin=int(input('EXCLUDE MIN: '))
                    small_excludemax=int(input('EXCLUDE MAX: '))
                    ew.equivalent_width(specinfo, bandloc, small_xmin, small_xmax, small_excludemin, small_excludemax, name=specname, fldr=path, mc=False, interactive=False)
                    subprocess.call(['open', path+'/'+specname+'_EWfit.pdf'])
                    print('Re-apply parameters for this individual spectrum?')
                    smallparam_yn=input('y or (n): ')
                    if smallparam_yn=='':
                        smallparam_yn='n'
                print('Save this spectrum according to these parameters?')
                unique_yn=input('y or (n): ')
                if unique_yn=='y':
                    print('This spectrum and its output will now be filed as fitted')
                    print('File Saved.')
                    print('PARAMS: '+str([small_xmin, small_xmax, small_excludemin, small_excludemax]))
                    comment=input('Any comments? [Press enter if none]\n')
                    f=open(path+'/EW_FITS.csv', 'a')
                    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                    writer.writerow({'SPEC':specname, 'BANDLOC':bandloc, 'XMIN':small_xmin, 'XMAX':small_xmax, 'EXCLUDEMIN':small_excludemin, 'EXCLUDEMAX':small_excludemax, 'TYPE':'FITTED', 'COMMENTS':comment})
                    f.close()
                    amount_in=amount_in-1
                    determined.append(spec)
                    continue
                if unique_yn=='':
                    unique_yn='n'
                if unique_yn=='n':
                    print('Is the spectrum flat and EW effectively zero?')
                    zero_yn=input('y or (n): ')
                    if zero_yn=='':
                        zero_yn='n'
                    if zero_yn=='y':
                        print('This spectrum will now be filed as having ZERO EW')
                        comment=input('Any comments? [Press enter if none]\n')
                        f=open(path+'/EW_FITS.csv', 'a')
                        writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                        writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'ZERO', 'COMMENTS':comment})
                        f.close()
                        amount_in=amount_in-1
                        determined.append(spec)
                        continue
                    if zero_yn=='n':
                        print('Attempt to fit this spectrum again?')
                        hope_yn=input('y or (n): ')
                        if hope_yn=='':
                            hope_yn='n'
                        if hope_yn=='n':
                            print('This spectrum will be filed as undetermined')
                            comment=input('Any comments? [Press enter if none]\n')
                            f=open(path+'/EW_FITS.csv', 'a')
                            writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                            writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'UNDETERMINED', 'COMMENTS':comment})
                            f.close()
                            amount_in=amount_in-1
                            determined.append(spec)
                            continue
                        if hope_yn=='y':
                            print('The spectrum will remain in the pool for further review')
                            continue

    return()

def mc_run(path):
    '''
    (Created By Stanislav DeLaurentiis 2020)
    
    An automated loop that applies the MC routine using EW.equivalent_width (with mc=True and interactive=False) to all fits files that had accurate fits.
    Writes a csv file containing the EW values from the routine.
    
    Parameters
    ----------
    path : string
        The absolute path of the directory containing all fits files.

    Returns
    -------
    data1: dict
        A dictionary containing the spec name, and corresponding ew, and ew error.
    
    Outputs
    -------
    EW_VALS_ALL.csv

    '''
    data1={}
    data0={}
    
    f=open(path+'/EW_VALS_ALL.csv', 'w')
    writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EW_MU', 'EW_SIG'])
    writer.writeheader()
    f.close()
    
    specparams={}
    
    f=open(path+'/EW_FITS.csv', 'rU')
    reader=csv.DictReader(f)
    for line in reader:
        specname=line['SPEC']
        typee=line['TYPE']
        if typee=='FITTED':
            bandloc=int(line['BANDLOC'])
            xmin=int(line['XMIN'])
            xmax=int(line['XMAX'])
            excludemin=int(line['EXCLUDEMIN'])
            excludemax=int(line['EXCLUDEMAX'])
            specparams[specname]=[bandloc, xmin, xmax, excludemin, excludemax]
        if typee=='ZERO':
            data0[specname]=[0, 0]
            data1[specname]=[0, 0]
            
    f.close()       
            
    for specname in specparams:
            bandloc=specparams[specname][0]
            xmin=specparams[specname][1]
            xmax=specparams[specname][2]
            excludemin=specparams[specname][3]
            excludemax=specparams[specname][4]
            specinfo=rs.read_spec(path+'/'+specname+'.fits')
            histdata=ew.equivalent_width(specinfo, bandloc, xmin, xmax, excludemin, excludemax, name=specname, mc=True, interactive=False, clobber=True,fldr=path, n=10)
            eqwmu=float(histdata[0])
            eqwsig=float(histdata[1])
            print('\n')
            print('EW: '+str(eqwmu)+', EW_SIG: '+str(eqwsig))
            print('\n')
            f=open(path+'/EW_VALS_ALL.csv', 'a')
            writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EW_MU', 'EW_SIG'])
            writer.writerow({'SPECNAME':specname, 'EW_MU':eqwmu, 'EW_SIG':eqwsig})
            f.close()
            data1[specname]=[eqwmu, eqwsig]
    
    for line in data0:
        spec=line
        eqwmu=float(data0[line][0])
        eqwsig=float(data0[line][1])
        f=open(path+'/EW_VALS_ALL.csv', 'a')
        writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EW_MU', 'EW_SIG'])
        writer.writerow({'SPECNAME':spec, 'EW_MU':eqwmu, 'EW_SIG':eqwsig})
        f.close()
        
    return(data1)




