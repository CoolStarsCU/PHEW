#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 17:14:05 2021

@author: stanislavdelaurentiis
"""

import EW as ew
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import antools3 as at

path='/Users/stanislavdelaurentiis/Desktop/TEST_PHEW/'

spec_count=0
determined=[]
for spec in os.listdir(path):
    spec_count=spec_count+1

f=open('EQW_FITS.csv', 'w')
writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX'])
writer.writeheader()
f.close()
    
while spec_count>0:
    bandloc=int(input('BANDLOC: '))
    xmin=int(input('XMIN: '))
    xmax=int(input('XMAX: '))
    exclude_min=int(input('EXCLUDE MIN: '))
    exclude_max=int(input('EXCLUDE MAX: '))
    print('\n')
    for spec in os.listdir(path):
        if spec in determined:
            continue
        specname=spec.split('.')[0]
        try:
            specinfo=at.read_spec(path+spec)[0]
            ew.equivalent_width(specinfo, bandloc, xmin, xmax, exclude_min, exclude_max, mcmc=False, interactive=False, name=specname)
        except(AttributeError, ValueError):
            print(str(specname)+' is a FAULTED SPECTRUM')
            spec_count=spec_count-1
            determined.append(spec)
            continue
        print('\n')
        print(specname)
        verdict=input('KEEP FIT (y or (n))? ')
        if verdict=='' or verdict=='n':
            pass
        if verdict=='y':
            f=open('EQW_FITS.csv', 'a')
            writer=csv.DictWriter(f, fieldnames=['SPEC', 'BANDLOC', 'XMIN', 'XMAX', 'EXCLUDEMIN', 'EXCLUDEMAX'])
            writer.writerow({'SPEC':specname, 'BANDLOC':bandloc, 'XMIN':xmin, 'XMAX':xmax, 'EXCLUDEMIN':exclude_min, 'EXCLUDEMAX':exclude_max})
            f.close()
            spec_count=spec_count-1
            determined.append(spec)

f=open('EQW_VALS_ALL.csv', 'w')
writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EQW_MU', 'EQW_SIG'])
writer.writeheader()
f.close()
    


specparams={}
    
f=open('EQW_FITS.csv', 'rU')
reader=csv.DictReader(f)
for line in reader:
    specname=line['SPEC']
    bandloc=int(line['BANDLOC'])
    xmin=int(line['XMIN'])
    xmax=int(line['XMAX'])
    excludemin=int(line['EXCLUDEMIN'])
    excludemax=int(line['EXCLUDEMAX'])
    specparams[specname]=[bandloc, xmin, xmax, excludemin, excludemax]

f.close()



for specname in specparams:
    bandloc=specparams[specname][0]
    xmin=specparams[specname][1]
    xmax=specparams[specname][2]
    excludemin=specparams[specname][3]
    excludemax=specparams[specname][4]
    histdata=ew.equivalent_width(path+'/'+specname+'.fits', bandloc, xmin, xmax, excludemin, excludemax, name=specname, mcmc=True, interactive=False, clobber=True)
    eqwmu=float(histdata[0])
    eqwsig=float(histdata[1])
    print('EQW: '+str(eqwmu)+', EQW_SIG: '+str(eqwsig))
    f=open('EQW_VALS_ALL.csv', 'a')
    writer=csv.DictWriter(f, fieldnames=['SPECNAME', 'EQW_MU', 'EQW_SIG'])
    writer.writerow({'SPECNAME':specname, 'EQW_MU':eqwmu, 'EQW_SIG':eqwsig})
    f.close()











            