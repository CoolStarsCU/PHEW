#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 11:37:13 2021

@author: stanislavdelaurentiis
"""

import antools3 as at
import csv
import os
import subprocess
import shutil
import EW as ew


def fit_runner(path):
    '''
    (Created By Stanislav DeLaurentiis, 2020)
    
    Will run through all the fits files in the inputted directory and
    apply antools.equivalent_width (with errors=False) to each file, ensuring that
    the curve is fit well to the spectral feature.
    Will then organize each fits file and its output pdf figure into a directory
    correlating with its parameters, while recording specname, paramters, and 
    type in 'EQW_FITS.csv. 
    

    Parameters
    ----------
    path : string
        The absolute path of the directory that contains the deisred fits files and is ok'd to be drained.

    Returns
    -------
    None.
    
    Outputs
    --------
    EQW_FITS.csv

    '''
    
    f=open('EQW_FITS.csv', 'w')
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
                specname=spec.split('.')[0]
                try:
                    specinfo=at.read_spec(path+"/"+spec)[0] #reading the fits file into array
                    ew.equivalent_width(specinfo, bandloc, xmin, xmax, excludemin, excludemax, name=specname, fldr=path, mcmc=False, interactive=False)
                except(AttributeError, ValueError):
                    print('ERROR: Unable to read this spectrum.')
                    f=open('EQW_FITS.csv', 'a')
                    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                    writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'UNABLE TO READ', 'COMMENTS':''})
                    f.close()
                    amount_in=amount_in-1
                    determined.append(spec)
                    continue
                subprocess.call(['open', path+'/'+specname+'_EWfit.pdf'])
                print('Would you like to save this fit using these parameters?')
                first_yn=input('y or (n): ')
                if first_yn=='':
                    first_yn='n'
                if first_yn=='y':
                    print('File Saved.')
                    print('PARAMS: '+str([xmin, xmax, excludemin, excludemax]))
                    comment=input('Any comments? [Press enter if none]\n')
                    f=open('EQW_FITS.csv', 'a')
                    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                    writer.writerow({'SPEC':specname, 'BANDLOC':bandloc, 'XMIN':xmin, 'XMAX':xmax, 'EXCLUDEMIN':excludemin, 'EXCLUDEMAX':excludemax, 'TYPE':'FITTED', 'COMMENTS':comment})
                    f.close()
                    determined.append(spec)
                    amount_in=amount_in-1
                    continue
                print('Would you like to re-apply parameters for this individual spectrum?')
                smallparam_yn=input('y or (n): ')
                if smallparam_yn=='':
                    smallparam_yn='n'
                if smallparam_yn=='n':
                    print('Is the spectrum flat and ew effectively zero?')
                    zero_yn=input('y or (n): ')
                    if zero_yn=='':
                        zero_yn='n'
                    if zero_yn=='y':
                        print('This spectrum will now be filed as having ZERO ew')
                        comment=input('Any comments? [Press enter if none]\n')
                        f=open('EQW_FITS.csv', 'a')
                        writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                        writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'ZERO', 'COMMENTS':comment})
                        f.close()
                        determined.append(spec)
                        amount_in=amount_in-1
                        continue
                    if zero_yn=='n':
                        print('Would you like to continue attempting to fit this spectrum?')
                        hope_yn=input('y or (n): ')
                        if hope_yn=='':
                            hope_yn='n'
                        if hope_yn=='n':
                            print('This spectrum will be filed as undetermined')
                            comment=input('Any comments? [Press enter if none]\n')
                            f=open('EQW_FITS.csv', 'a')
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
                    print('Please enter your choice of unique parameters')
                    small_xmin=int(input('XMIN: '))
                    small_xmax=int(input('XMAX: '))
                    small_excludemin=int(input('EXCLUDE MIN: '))
                    small_excludemax=int(input('EXCLUDE MAX: '))
                    ew.equivalent_width(specinfo, bandloc, small_xmin, small_xmax, small_excludemin, small_excludemax, name=specname, fldr=path, mcmc=False, interactive=False)
                    subprocess.call(['open', path+'/'+specname+'_EWfit.pdf'])
                    print('Would you like to re-apply parameters for this individual spectrum?')
                    smallparam_yn=input('y or (n): ')
                    if smallparam_yn=='':
                        smallparam_yn='n'
                print('Would you like to save this spectrum according to these parameters?')
                unique_yn=input('y or (n): ')
                if unique_yn=='y':
                    print('This spectrum and its output will now be filed as fitted')
                    print('File Saved.')
                    print('PARAMS: '+str([small_xmin, small_xmax, small_excludemin, small_excludemax]))
                    comment=input('Any comments? [Press enter if none]\n')
                    f=open('EQW_FITS.csv', 'a')
                    writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                    writer.writerow({'SPEC':specname, 'BANDLOC':bandloc, 'XMIN':small_xmin, 'XMAX':small_xmax, 'EXCLUDEMIN':small_excludemin, 'EXCLUDEMAX':small_excludemax, 'TYPE':'FITTED', 'COMMENTS':comment})
                    f.close()
                    amount_in=amount_in-1
                    determined.append(spec)
                    continue
                if unique_yn=='':
                    unique_yn='n'
                if unique_yn=='n':
                    print('Is the spectrum flat and ew effectively zero?')
                    zero_yn=input('y or (n): ')
                    if zero_yn=='':
                        zero_yn='n'
                    if zero_yn=='y':
                        print('This spectrum will now be filed as having ZERO ew')
                        comment=input('Any comments? [Press enter if none]\n')
                        f=open('EQW_FITS.csv', 'a')
                        writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX', 'TYPE', 'COMMENTS'])
                        writer.writerow({'SPEC':specname, 'BANDLOC':'', 'XMIN':'', 'XMAX':'', 'EXCLUDEMIN':'', 'EXCLUDEMAX':'', 'TYPE':'ZERO', 'COMMENTS':comment})
                        f.close()
                        amount_in=amount_in-1
                        determined.append(spec)
                        continue
                    if zero_yn=='n':
                        print('Would you like to further attempt to fit this spectrum?')
                        hope_yn=input('y or (n): ')
                        if hope_yn=='':
                            hope_yn='n'
                        if hope_yn=='n':
                            print('This spectrum will be filed as undetermined')
                            comment=input('Any comments? [Press enter if none]\n')
                            f=open('EQW_FITS.csv', 'a')
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

def mcmc_run(path):
    '''
    (Created By Stanislav DeLaurentiis 2020)
    
    An automated loop that applies the MCMC routine to all fits files
    that had accurate fits.
    Writes a csv file contianing the eqw values from the routine.
    
    (Needs to be used with antools2.py)

    Parameters
    ----------
    path : string
        The absolute path of the directory containing all fits files.
    steps : int
        How many steps each MCMC routine runs through (recommended=1000).

    Returns
    -------
    None.
    
    Outputs
    -------
    EQW_VALS_ALL.csv

    '''
    data1={}
    data0={}

    f=open('EQW_VALS_ALL.csv', 'w')
    writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EQW_MU', 'EQW_SIG'])
    writer.writeheader()
    f.close()
    
    specparams={}
    
    f=open('EQW_FITS.csv', 'rU')
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
            specinfo=at.read_spec(path+'/'+specname+'.fits')[0]
            histdata=ew.equivalent_width(specinfo, bandloc, xmin, xmax, excludemin, excludemax, name=specname, mcmc=True, interactive=False, clobber=True,fldr=path, n=10)
            eqwmu=float(histdata[0])
            eqwsig=float(histdata[1])
            print('\n')
            print('EQW: '+str(eqwmu)+', EQW_SIG: '+str(eqwsig))
            print('\n')
            f=open('EQW_VALS_ALL.csv', 'a')
            writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EQW_MU', 'EQW_SIG'])
            writer.writerow({'SPECNAME':specname, 'EQW_MU':eqwmu, 'EQW_SIG':eqwsig})
            f.close()
            data1[specname]=[eqwmu, eqwsig]
    
    for line in data0:
        spec=line
        eqwmu=float(data0[line][0])
        eqwsig=float(data0[line][1])
        f=open('EQW_VALS_ALL.csv', 'a')
        writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EQW_MU', 'EQW_SIG'])
        writer.writerow({'SPECNAME':spec, 'EQW_MU':eqwmu, 'EQW_SIG':eqwsig})
        f.close()
        
    return(data1)




