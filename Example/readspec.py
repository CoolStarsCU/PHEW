import numpy as np
from astropy.io import ascii
import astropy.io.fits as pf

def read_spec(specFiles, errors=True, negtonan=False, linear=False, ext=0, verbose=True, header=False):
    '''
    (by Alejandro N |uacute| |ntilde| ez, Jocelyn Ferrara)

    Read spectral data from fits or ascii files. It returns a list of numpy arrays with wavelength in position 0, flux in position 1 and error values (if requested) in position 2. More than one file name can be provided simultaneously.

    **Limitations**: Due to a lack of set framework for ascii file headers, this function assumes ascii files to have wavelength in column 1, flux in column 2, and (optional) error in column 3. Ascii spectra are assumed to be linear, so the kwarg *linear* is disabled for ascii files. Fits files that have multiple spectral orders will not be interpreted correctly with this function.

    *specFiles*
      String with fits file name (with full path); it can also be a python list of file names.
    *errors*
      Boolean, whether to return error values for the flux data; return nans if unavailable.
    *negtonan*
      Boolean, whether to set negative flux values equal to zero.
    *linear*
      Boolean, whether to return spectrum only if it is linear. If it cannot verify linearity, it will assume linearity.
    *ext*
      Integer, the extension number of the fits file where the spectrum is.
    *verbose*
      Boolean, whether to print warning messages.
    *header*
      Boolean, whether to also return the fits file header.
    '''

    # 1. Convert specFiles into a list type if it is only one file name
    if isinstance(specFiles, str):
        specFiles = [specFiles,]

    try:
        specFiles[0]
    except TypeError:
        print('ERROR: File name(s) in invalid format.')
        return

    # 2. Initialize array to store spectra
    specData = [None] * len(specFiles)

    # 3. Loop through each file name:
    for spFileIdx,spFile in enumerate(specFiles):
        if spFile is None: continue
        # 3.1 Determine the type of file it is
        isFits = False
        extension = spFile[-4:].lower()
        if extension == 'fits' or extension == '.fit':
            isFits = True
        
        # 3.2. Get data from file
        if isFits:
            isSDSS = False
            isHARPS = False
            # Determine table index to extract the data
            if ext != 0:
                tmpext = ext
                fitsData = pf.getdata(spFile, ext=tmpext)
            elif ext == 0:
                # Even if specified by user, try to determine the relevant extension
                try:
                    tmpHead = pf.getheader(spFile, ext=0)
                    
                    # Telescope exceptions
                    try:
                        tmptelescope = tmpHead['TELESCOP'].upper()
                    except KeyError:
                        tmptelescope = ''
                    if tmptelescope.find('SDSS') != -1:
                        isSDSS = True
                        tmpext = 1
                    if tmptelescope.find('LAMOST') != -1:
                        isLAMOST = True
                    
                    if not isSDSS:
                        if tmpHead['NAXIS'] == 0:
                            try:
                                if tmpHead['NAXIS1'] < 100:
                                    tmpext = 2
                                else:
                                    tmpext = 1
                            except KeyError:
                                tmpext = 1
                        else:
                            tmpext = 0
                    fitsData = pf.getdata(spFile, ext=tmpext)
                except IOError:
                    print('Could not open ' + str(spFile) + '.')
                    continue
            # Re-shape SDSS data array to make it compatible with the rest of this code
            if isSDSS:
                fitsData = np.array(fitsData.tolist()).T
            
            # Now determine the table index to extract header info with wavelength solution
            tmpHead = pf.getheader(spFile, ext=tmpext)
            if isSDSS:
                fitsHeader = pf.getheader(spFile, ext=0) 
            else:
                fitsHeader = tmpHead.copy()

            # Determine if spectrum is from ESO Observatory
            try:
                tmpHead = pf.getheader(spFile, ext=0)
                tmpfind = tmpHead['TELESCOP'].find('ESO')
            except KeyError:
                tmpfind = -1
            if tmpfind != -1:
                isHARPS = True

        # Assume ascii file otherwise
        else:
            try:
                aData = ascii.read(spFile)
                specData[spFileIdx] = [aData[0].tonumpy(), aData[1].tonumpy()]
                if len(aData) >= 3 and errors:
                    specData[spFileIdx].append(aData[2].tonumpy())
            except IOError:
                print('ERROR: Could not open ' + str(spFile) + '.')
                continue

        # 3.3. Check if data in fits file is linear
        if isFits:
            KEY_TYPE = ['CTYPE1']
            setType  = set(KEY_TYPE).intersection(set(fitsHeader.keys()))
            if len(setType) == 0:
                if verbose:
                    print('WARNING: Data in ' + spFile + ' assumed to be linear.')
                isLinear = True
            else:
                valType = fitsHeader[setType.pop()]
                if valType.strip().upper() == 'LINEAR':
                    isLinear = True
                else:
                    isLinear = False
            if linear and not isLinear:
                if verbose:
                    print('ERROR: Data in ' + spFile + ' is not linear.')
                return

        # 3.4. Get wl, flux & error data from fits file
        #      (returns wl in pos. 0, flux in pos. 1, error values in pos. 2)
        if isFits:
            specData[spFileIdx] = __get_spec(fitsData, fitsHeader, spFile, errors, \
                                             isHARPS=isHARPS, verb=verbose)
            if specData[spFileIdx] is None:
                continue

            # Generate wl axis when needed
            if specData[spFileIdx][0] is None:
                specData[spFileIdx][0] = __create_waxis(fitsHeader, \
                                         len(specData[spFileIdx][1]), spFile, \
                                         verb=verbose)
            # If no wl axis generated, then clear out all retrieved data for object
            if specData[spFileIdx][0] is None:
                specData[spFileIdx] = None
                continue

        # 3.5. Set negative flux values equal to zero (next step sets them to nans)
        if negtonan:
            negIdx = np.where(specData[spFileIdx][1] < 0)
            if len(negIdx[0]) > 0:
                specData[spFileIdx][1][negIdx] = 0
                if verbose:
                    print('%i negative data points found in %s.' \
                            % (len(negIdx[0]), spFile))

        # 3.6. Set zero flux values as nans (do this always)
        zeros = np.where(specData[spFileIdx][1] == 0)
        if len(zeros[0]) > 0:
            specData[spFileIdx][1][zeros] = np.nan

    # 4. Simplify output if only one spectrum retrieved
    if len(specData) == 1:
        specData = specData[0]

    if header:
        return specData, fitsHeader
    else:
        return specData

def __create_waxis(fitsHeader, lenData, fileName, verb=True):
    # Generates a wavelength (wl) axis using header data from fits file.

    # Define key names in
    KEY_MIN  = ['COEFF0','CRVAL1']         # Min wl
    KEY_DELT = ['COEFF1','CDELT1','CD1_1'] # Delta of wl
    KEY_OFF  = ['LTV1']                    # Offset in wl to subsection start

    # Find key names for minimum wl, delta, and wl offset in fits header
    setMin  = set(KEY_MIN).intersection(set(fitsHeader.keys()))
    setDelt = set(KEY_DELT).intersection(set(fitsHeader.keys()))
    setOff  = set(KEY_OFF).intersection(set(fitsHeader.keys()))

    # Get the values for minimum wl, delta, and wl offset, and generate axis
    if len(setMin) >= 1 and len (setDelt) >= 1:
        nameMin = setMin.pop()
        valMin  = fitsHeader[nameMin]

        nameDelt = setDelt.pop()
        valDelt  = fitsHeader[nameDelt]

        if len(setOff) == 0:
            valOff = 0
        else:
            nameOff = setOff.pop()
            valOff  = fitsHeader[nameOff]

        # generate wl axis
        if nameMin == 'COEFF0':
            wAxis = 10 ** (np.arange(lenData) * valDelt + valMin)
        else:
            wAxis = (np.arange(lenData) * valDelt) + valMin - (valOff * valDelt)
    else:
        wAxis = None
        if verb:
            print('ERROR: Could not re-create wavelength axis for ' + fileName + '.')

    return wAxis

def __get_spec(fitsData, fitsHeader, fileName, errorVals, isHARPS, verb=True):
    # Interprets spectral data from fits file.
    # Returns wavelength (wl) data in pos. 0, flux data in pos. 1, and if requested, error values in pos. 2.
    
    if errorVals:
        validData = [None] * 3
    else:
        validData = [None] * 2

    # Fetch spectral data from ESO spectra
    if isHARPS:
        validData = [fitsData['WAVE'][0], fitsData['FLUX'][0], fitsData['ERR'][0]]
        return validData

    # Identify number of data sets in fits file
    dimNum = len(fitsData)

    fluxIdx  = None
    waveIdx  = None
    sigmaIdx = None
    isSDSS = False
    try:
        if fitsHeader['TELESCOP'].upper().find('LAMOST') != -1:
            isLAMOST = True
        else:
            isLAMOST = False
    except KeyError:
        isLAMOST = False
    
    # Identify data sets in fits file
    if dimNum == 1:
        fluxIdx = 0
    elif dimNum == 2:
        if len(fitsData[0]) == 1:
            sampleData = fitsData[0][0][20]
        else:
            sampleData = fitsData[0][20]
        if sampleData < 0.0001:
            # 0-flux, 1-unknown
            fluxIdx  = 0
        else:
            waveIdx = 0
            fluxIdx = 1
    elif dimNum == 3:
        waveIdx  = 0
        fluxIdx  = 1
        sigmaIdx = 2
    elif dimNum == 4:
    # 0-flux clean, 1-flux raw, 2-background, 3-sigma clean
        fluxIdx  = 0
        sigmaIdx = 3
    elif dimNum == 5:
        if isLAMOST:
            # 0-flux, 1-inv.var, 2-wl, 3-andmask, 4-ormask
            fluxIdx = 0
            sigmaIdx = 1 # This column is actually 1/sigma^2
            waveIdx = 2
        else:
            # 0-flux, 1-continuum substracted flux, 2-sigma, 3-mask array, 4-unknown
            fluxIdx  = 0
            sigmaIdx = 2
    elif dimNum == 8:
        # SDSS spectra
        fluxIdx = 0
        waveIdx = 1 # This column is actually log10(wl)
        sigmaIdx = 2 # This column is actually 1/sigma^2
        isSDSS = True
    elif dimNum > 10:
    # Implies that only one data set in fits file: flux
        fluxIdx = -1
        if np.isscalar(fitsData[0]):
            fluxIdx = -1
        elif len(fitsData[0]) in range(2,4):
        # Data comes in a xxxx by 2 matrix (ascii origin)
            tmpWave = []
            tmpFlux = []
            if len(fitsData[0]) == 3:
                tmpeFlux = []
            for pair in fitsData:
                tmpWave.append(pair[0])
                tmpFlux.append(pair[1])
                if len(fitsData[0]) == 3:
                    tmpeFlux.append(pair[2])
            if len(fitsData[0]) == 3:
                fitsData = [tmpWave, tmpFlux, tmpeFlux]
            else:
                fitsData = [tmpWave,tmpFlux]
            fitsData = np.array(fitsData)

            waveIdx = 0
            fluxIdx = 1
        else:
        # Indicates that data is structured in an unrecognized way
            fluxIdx = None
    else:
        fluxIdx = None

    # Fetch wave data set from fits file
    if fluxIdx is None:
    # No interpretation known for fits file data sets
        validData = None
        if verb:
            print('ERROR: Unable to interpret data in ' + fileName + '.')
        return validData
    else:
        if waveIdx is not None:
            if len(fitsData[waveIdx]) == 1:
            # Data set may be a 1-item list
                validData[0] = fitsData[waveIdx][0]
            else:
                if isSDSS:
                    validData[0] = 10**fitsData[waveIdx]
                else:
                    validData[0] = fitsData[waveIdx]
                # Convert from vacuum wl to air wl
                if isSDSS or isLAMOST:
                    validData[0] = validData[0] / (1.0 +  5.792105E-2/(238.0185 \
                                   - (1E4/validData[0])**2) + 1.67917E-3/(57.362 \
                                   - (1E4/validData[0])**2))

    # Fetch flux data set from fits file
    if fluxIdx == -1:
        validData[1] = fitsData
    else:
        if len(fitsData[fluxIdx]) == 1:
            validData[1] = fitsData[fluxIdx][0]
        else:
            validData[1] = fitsData[fluxIdx]
            if isSDSS:
                validData[1] = validData[1] * 1E-17

    # Fetch sigma data set from fits file, if requested
    if errorVals:
        if sigmaIdx is None:
            validData[2] = np.array([np.nan] * len(validData[1]))
        else:
            if len(fitsData[sigmaIdx]) == 1:
                validData[2] = fitsData[sigmaIdx][0]
            else:
                if isSDSS or isLAMOST:
                    validData[2] = 1 / np.sqrt(fitsData[sigmaIdx])
                else:
                    validData[2] = fitsData[sigmaIdx]
                if isSDSS:
                    validData[2] = validData[2] * 1E-17

        # If all sigma values have the same value, replace them with nans
        if np.nanmin(validData[2]) == np.nanmax(validData[2]):
            validData[2] = np.array([np.nan] * len(validData[1]))
    
    # Check ascending order of spectrum using wavelength axis
    if validData[0] is not None:
        if validData[0][0] > validData[0][-1]:
            for i in range(len(validData)):
                if validData[i] is not None:
                    validData[i] = validData[i][::-1]

    return validData
