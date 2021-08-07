import numpy as np
from astropy.io import ascii
import types
import scipy.interpolate as spi
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import pdb

def read_table(name, delimiter='\t', comment='#', fmt=None, ds=1):
    '''
    Reads ascii tables and converts them cleanly into numpy arrays.
    '''
    if fmt is not None:
        datanp = ascii.read(name, guess=False, delimiter=delimiter, \
                            comment=comment, header_start=0, \
                            data_start=ds, format=fmt)
    else:
        datanp = ascii.read(name, guess=False, delimiter=delimiter, \
                            comment=comment, header_start=0, \
                            data_start=ds)

    return datanp

def angsep(ra1deg, dec1deg, ra2deg, dec2deg, angle=False):
    '''
    Determine separation in degrees between celestial objects.
       ra1deg, dec1deg - primary point(s); can be arrays
       ra2deg, dec2deg - secondary point(s); can be arrays or scalars
       angle - if True, it will calculate angle E of N.
       All arguments are in decimal degrees.
       Returns distance in arcdegrees, angles between -180 and 180 degrees.
    '''
    ra1rad = ra1deg * np.pi / 180
    dec1rad = dec1deg * np.pi / 180
    ra2rad = ra2deg * np.pi / 180
    dec2rad = dec2deg * np.pi / 180

    # calculate scalar product for determination of angular separation
    x = np.cos(ra1rad) * np.cos(dec1rad) * np.cos(ra2rad) * np.cos(dec2rad)
    y = np.sin(ra1rad) * np.cos(dec1rad) * np.sin(ra2rad) * np.cos(dec2rad)
    z = np.sin(dec1rad) * np.sin(dec2rad)

    rad = np.arccos(x + y + z) # Sometimes gives warnings when coords match

    # use Pythagoras approximation if rad < 1 arcsec
    sep = np.choose(rad<0.000004848, (np.sqrt((np.cos(dec1rad) * (ra1rad-ra2rad))**2 \
                    + (dec1rad - dec2rad)**2), rad))

    # Angular separation
    sep = sep * 180 / np.pi

    if angle:
        deltaDEC = dec1rad - dec2rad
        deltaRA = ra1rad - ra2rad
        angledeg = np.arctan2(-deltaRA, -deltaDEC) * 180 / np.pi

        return sep, angledeg
    else:
        return sep


def deg2sex(ras, decs):
    ''' Converts RA and DEC from decimal to sexagesimal. Returns string.
        Arguments:
          ras - string(s) of RA in degrees
          decs - string(s) of DEC in degrees
    '''
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    if type(ras) == list or type(ras) == np.ndarray:
        new_coords = []
        for irow in range(0,len(ras)):
            c = SkyCoord(float(ras[irow]), float(decs[irow]), \
                         frame='icrs', unit='deg')
            new_coords.append(c.to_string('hmsdms'))

    else:
        c = SkyCoord(float(ras), float(decs), frame='icrs', unit='deg')
        new_coords = c.to_string('hmsdms')

    return new_coords

def sex2deg(ras, decs):
    ''' Converts RA and DEC from sexagesimal to decimal.
        Arguments:
          ras - string(s) of RA in sexagesimal degrees (HH MM SS.SS)
          decs - string(s) of DEC in sexagesimal degrees (+-DD MM SS.SS)
    '''
    if type(ras) == list or type(ras) == np.ndarray:
        new_ras = []
        new_decs = []
        for irow in range(0,len(ras)):
            parts_ra = ras[irow].rsplit(' ')
            if len(parts_ra) == 1:
                parts_ra = ras[irow].rsplit(':')
            parts_dec = decs[irow].rsplit(' ')
            if len(parts_dec) == 1:
                parts_dec = decss[irow].rsplit(':')
            ra_deg = float(parts_ra[0]) * 15. + float(parts_ra[1]) / 4. + float(parts_ra[2]) / 240.
            dec_deg = float(parts_dec[0]) + float(parts_dec[1]) / 60. + float(parts_dec[2]) / 3600.
            new_ras.append(ra_deg)
            new_decs.append(dec_deg)
        new_ras = np.array(new_ras)
        new_decs = np.array(new_decs)

        return new_ras, new_decs

    else:
        parts_ra = ras.rsplit(' ')
        if len(parts_ra) == 1:
            parts_ra = ras.rsplit(':')
        parts_dec = decs.rsplit(' ')
        if len(parts_dec) == 1:
            parts_dec = decs.rsplit(':')
        ra_deg = float(parts_ra[0]) * 15. + float(parts_ra[1]) / 4. + float(parts_ra[2]) / 240.
        dec_deg = float(parts_dec[0]) + float(parts_dec[1]) / 60. + float(parts_dec[2]) / 3600.

        return ra_deg, dec_deg

def matchsorted(ra, dec, ra1, dec1, tol, angle=False, closest=True):
    ''' Find closest ra,dec within tol to a target in an ra-sorted list of ra,dec.
        Arguments:
          ra - Right Ascension decimal degrees (numpy sorted in ascending order)
          dec - Declination decimal degrees (numpy array)
          ra1 - RA to match (scalar, decimal degrees)
          dec1 - Dec to match (scalar, decimal degrees)
          tol - Matching tolerance in arcseconds.
          angle - Boolean, whether to return angle formed by matched sources.
          closest - Boolean, whether to return only the closest match.
        Returns:
          ibest - index of the (best) match(es) within tol; -1 if no match within tol
          sep - separation (defaults to tol if no match within tol)
          angle - angle (defaults to 0 if no match within tol)
    '''
    tol = tol / 3600.
    if isinstance(tol, float):
        # Case for one tolerance radius for all objects
        i1 = np.searchsorted(ra, ra1 - tol) - 5
        i2 = np.searchsorted(ra, ra1 + tol) + 5
    else:
        # Case for one tolerance radius for each object
        i1 = np.searchsorted(ra + tol, ra1) - 5
        i2 = np.searchsorted(ra - tol, ra1) + 5
    if i1 < 0:
        i1 = 0
    if angle:
        sep, ang = angsep(ra[i1:i2], dec[i1:i2], ra1, dec1, angle=angle)
    else:
        sep = angsep(ra[i1:i2], dec[i1:i2], ra1, dec1, angle=angle)
    
    if isinstance(tol, float):
        imatches = np.where(sep < tol)[0]
    else:
        imatches = np.where(sep < tol[i1:i2])[0]
    if len(imatches) == 0:
        if angle:
            return [-1], [tol * 3600.], [0]
        else:
            return [-1], [tol * 3600.]

    ibest = np.argmin(sep[imatches])
    #indices = np.argsort(sep)
    #if sep[indices[0]] > tol:
    #    if angle:
    #        return -1, tol * 3600., 0
    #    else:
    #        return -1, tol * 3600.
    #ibest = indices[0] + i1
    #imult = indices[np.where(sep[indices] < tol)[0]] + i1
    #imult = np.where(sep < tol)[0]

    if angle:
        if closest:
            return [imatches[ibest] + i1], [sep[imatches][ibest] * 3600.], \
                   [ang[imatches[ibest]]]
        else:
            return imatches + i1, sep[imatches] * 3600., ang[imatches]
    else:
        if closest:
            return [imatches[ibest] + i1], [sep[imatches][ibest] * 3600.]
        else:
            return imatches + i1, sep[imatches] * 3600.

def smooth(x,window_len=11,window='hanning'):
    """
    smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1],x,x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    
    return y[int(window_len / 2 - 1):-int(window_len / 2)]


def mean_comb(spectra, weights=None, mask=None, robust=None, forcesimple=False, extremes=False, renormalize=False):
    '''
    (by Alejandro N |uacute| |ntilde| ez)

    Combine spectra using a (weighted) mean. The output is a python list with mask wavelength in position 0, mean flux in position 1, and variance in position 2. If flux uncertainties are given, then mean is a weighted mean, and variance is the "variance of the mean" (|sigma|  :sub:`mean`  :sup:`2`). If no flux uncertainties are given, then mean is a straight mean (<x>), and variance is the square of the standard error of the mean (|sigma| :sup:`2`/n). If no mask is given, the wavelength array of the first spectrum will be used as mask.

    This function mimics IDL mc_meancomb (by Mike Cushing), with some restrictions.

    *spectra*
        Python list of spectra, where each spectrum is an array having wavelength in position 0, flux in position 1, and optional uncertainties in position 2.
    *weights*
      List of weights corresponding to each spectrum (must add up to one). If none, then each spectrum is assumed to have the same weight. THIS ONLY WORKS IF I GIVE IT TWO SPECTRA.
    *mask*
      Array of wavelengths to be used as mask for all spectra. If none, then the wavelength array of the first spectrum is used as mask.
    *robust*
      Float, the sigma threshold to throw bad flux data points out. If none given, then all flux data points will be used.
    *forcesimple*
      Boolean, whether to calculate a straight mean and variance even if weights are available.
    *extremes*
      Boolean, whether to include the min and max flux values at each masked pixel.
    *renormalize*
      Boolean, whether to re-normalized the spectra agains the calculated combined spectrum, in which case the spectra will be returned in a list, with masked values.

    '''
    # Check inputs
    try:
        spectra[0]
    except TypeError:
        print('Spectra invalid.')
        return
    if mask is not None:
        try:
            mask[0]
        except TypeError:
            print('Mask invalid.')
            return
    if robust is not None:
        try:
            float(robust)
        except TypeError:
            print('Robust invalid.')
            return

    # 1. Generate mask using the first spectrum given
    if mask is None:
        # Use x-axis (i.e. wl) values of first spectrum as mask for all others
        wl_mask = spectra[0][0]
    else:
        wl_mask = mask
    numPoints = len(wl_mask)
    numSpec = len(spectra)

    # 2. Check if uncertainties were given
    uncsGiven = True
    if forcesimple:
        uncsGiven = False
    for spec in spectra:
        if uncsGiven:
            try:
                uncs = spec[2]
            except IndexError:
                uncsGiven = False
                continue
            nanIdx = np.where(np.isfinite(uncs))
            if len(nanIdx[0]) == 0:
                uncsGiven = False

    # 3D-array that will hold interpolated spectra
    # (it omits wavelength dimension, since all spectra have the same one)
    if uncsGiven:
        dims = 2
    else:
        dims = 1
    ip_spectra = np.zeros((numPoints, dims, numSpec)) * np.nan

    # 3. Interpolate spectra using mask
    for spIdx, spec in enumerate(spectra):
        wl = spec[0]
        fluxRaw= spec[1]
        if uncsGiven:
            unc = spec[2]

        # Eliminate outliers if requested
        if robust is not None:
            flux = clean_outliers(fluxRaw, robust)
        else:
            flux = fluxRaw

        if spIdx == 0:
            # No need to interpolate first spectrum
            flux_new = flux
            if uncsGiven:
                unc_new = unc
        else:
            ip_func_flux = spi.interp1d(wl, flux, bounds_error=False)
            flux_new = ip_func_flux(wl_mask.tolist())
            if uncsGiven:
                ip_func_unc = spi.interp1d(wl, unc, bounds_error=False)
                unc_new = ip_func_unc(wl_mask.tolist())

        ip_spectra[:,0,spIdx] = flux_new
        if uncsGiven:
            ip_spectra[:,1,spIdx] = unc_new

    # 4. Calculate mean and variance of flux values
    if weights is None:
        wgts = np.ones(len(spectra))
    else:
        wgts = weights
    if uncsGiven:
        mvarraw = 1. / np.nansum(1. / (wgts * ip_spectra[:,1,:]), axis=1) # 1/Sum(1/sigma_i^2)
        wmean = np.nansum(wgts * ip_spectra[:,0,:] / ip_spectra[:,1,:], axis=1) # Sum(x_i/sigma_i^2)
        mean = wmean * mvarraw
        mvar = mvarraw
        # Correct weighted sample variance for small sample
        #meantile = np.tile(mean, (numSpec,1)).T
        #V1 = 1 / mvarraw
        #V2 = np.nansum(ip_spectra[:,1,:]**2, axis=1)
        #mvar = V1 / (V1**2 - V2) * \
        #       np.nansum((ip_spectra[:,0,:] - meantile)**2 / ip_spectra[:,1,:], axis=1)
    else:
        mvar = np.nanstd(ip_spectra[:,0,:], axis=1) ** 2 # /numSpec -- I think I dont need this
        mean = np.nanmean(ip_spectra[:,0,:], axis=1)

    # 5. Calculate extreme flux values if requested
    if extremes:
        minF = np.nanmin(ip_spectra[:,0,:], axis=1)
        maxF = np.nanmax(ip_spectra[:,0,:], axis=1)

    # 5. Create the combined spectrum
    if extremes:
        specComb = [wl_mask, mean, mvar, minF, maxF]
    else:
        specComb = [wl_mask, mean, mvar]

    # 6. Re-normalize spectra to calculated combined spectrum, if requested
    if renormalize:
        renorm_spectra = []
        for ispec in range(0, numSpec):
            tmpflux = ip_spectra[:,0,ispec]
            renormfac = np.median(tmpflux / mean) # mean is the flux of the combined spectrum
            if uncsGiven:
                tmpunc = ip_spectra[:,1,ispec]
                renorm_spectra.append([wl_mask, tmpflux / renormfac, tmpunc / renormfac])
            else:
                renorm_spectra.append([wl_mask, tmpflux / renormfac])

        return specComb, renorm_spectra
    else:
        return specComb

def norm_spec(specData, limits, flag=False):
    '''
    (by Alejandro N |uacute| |ntilde| ez)

    Normalize a spectrum using a band (i.e. a portion) of the spectrum specified by *limits*.

    *specData*
      Spectrum as a Python list with wavelength in position 0, flux in position 1, and (optional) error values in position 2. More than one spectrum can be provided simultaneously, in which case *specData* shall be a list of lists.
    *limits*
      Python list with lower limit in position 0 and upper limit in position 1. If more than one spectrum provided, these limits will be applied to all spectra.
    *flag*
      Boolean, whether to warn if normalization limits were shrinked in the case when they fall outside spectrum. If set to *True*, *norm_spec* returns the normalized spectra AND a boolean flag.
    '''

    # Convert specData to list or spectra if it consists only of one
    if len(specData) <= 3 and len(specData[0]) > 10:
        specData = [specData]

    # Initialize objects
    finalData = [None] * len(specData)

    # Check that given limits are reasonable
    if limits[0] >= limits[1]:
        print('norm_spec: the Min and Max values specified are not reasonable.')
        return None

    # Re-define normalizing band (specified in limits) for each spectrum in case
    # the limits fall outside of the spectrum range
    all_lims = [None] * len(specData)
    flagged = False
    for spIdx, spData in enumerate(specData):
        smallest = limits[0]
        largest  = limits[1]
        if spData is None:
            continue

        tmpNans = np.where(np.isfinite(spData[1]))
        if len(tmpNans[0]) != 0:
            if spData[0][tmpNans[0][0]] > smallest:
                smallest = spData[0][tmpNans[0][0]]
                flagged = True
            if spData[0][tmpNans[0][-1]] < largest:
                largest = spData[0][tmpNans[0][-1]]
                flagged = True

        all_lims[spIdx] = [smallest, largest]
    lims = [smallest, largest]

    # Loop through each spectral data set
    for spIdx, spData in enumerate(specData):

        # 1) Skip if data is missing
        if spData is None:
            continue

        # 2) Determine if spectra come with error values
        if len(spData) == 3:
            errors = True
        else:
            errors = False

        # 3) Determine minimum wavelength value for band
        smallIdx = np.where(spData[0] < all_lims[spIdx][0])

        # If lower limit < all values in spectrum wavelength points, then
        # make band's minimum value = first data point in spectrum
        try:
            smallIdx[0]
        except IndexError:
            minIdx = 0
            smallIdx = [None]

        # If lower limit > all values in spectrum wavelength points, then
        # no band can be selected
        if smallIdx != [None]:
            if len(smallIdx[0]) == len(spData[0]):
                print('norm_spec: the wavelength data for object is outside limits.' )
                continue
            else:
                minIdx = smallIdx[0][-1] + 1

        # 4) Determine maximum wavelength value for band
        largeIdx = np.where(spData[0] > all_lims[spIdx][1])

        # If upper limit > all values in spectrum wavelength points, then
        # make band's maximum value = last data point in spectrum
        try:
            largeIdx[0]
        except IndexError:
            maxIdx = len(spData[0])
            largeIdx = [None]

        # If upper limit < all values in spectrum wavelength points, then
        # no band can be selected
        if largeIdx != [None]:
            if len(largeIdx[0]) == len(spData[0]):
                print('norm_spec: the wavelength data for object is outside limits.')
                continue
            else:
                maxIdx = largeIdx[0][0]

        # 5) Check for consistency in the computed band limits
        if maxIdx - minIdx < 2:
            print('norm_spec: The Min and Max values specified yield no band.')
            continue

        # 6) Select flux band from spectrum
        fluxSelect = spData[1][minIdx:maxIdx]
        fluxSelect = np.array(fluxSelect)

        # 7) Select error value band from spectrum
        if errors is True:
            errorSelect = spData[2][minIdx:maxIdx]
            errorSelect = np.array(errorSelect)

        # 8) Normalize spectrum using arithmetic mean
        notNans = np.where(np.isfinite(fluxSelect))
        avgFlux = np.mean(fluxSelect[notNans])
        finalFlux = spData[1] / avgFlux

        finalData[spIdx] = [spData[0], finalFlux]

        if errors is True:
            #notNans  = np.where(np.isfinite(errorSelect))
            #avgError = np.mean(errorSelect[notNans])
            finalErrors = spData[2] / avgFlux

            finalData[spIdx] = [spData[0], finalFlux, finalErrors]

    if flag:
        return finalData, flagged
    else:
        return finalData

def read_spec(specFiles, errors=True, atomicron=False, negtonan=False, plot=False, linear=False, templ=False, verbose=True, header=False):
    '''
    (by Alejandro N |uacute| |ntilde| ez, Jocelyn Ferrara)

    Read spectral data from fits or ascii files. It returns a list of numpy arrays with wavelength in position 0, flux in position 1 and error values (if requested) in position 2. More than one file name can be provided simultaneously.

    **Limitations**: Due to a lack of set framework for ascii file headers, this function assumes ascii files to have wavelength in column 1, flux in column 2, and (optional) error in column 3. Ascii spectra are assumed to be linear, so the kwarg *linear* is disabled for ascii files. Fits files that have multiple spectral orders will not be interpreted correctly with this function.

    *specFiles*
      String with fits file name (with full path); it can also be a python list of file names.
    *errors*
      Boolean, whether to return error values for the flux data; return nans if unavailable.
    *atomicron*
      Boolean, if wavelength units are in Angstrom, whether to convert them to microns.
    *negtonan*
      Boolean, whether to set negative flux values equal to zero.
    *plot*
      Boolean, whether to plot the spectral data, including error bars when available.
    *linear*
      Boolean, whether to return spectrum only if it is linear. If it cannot verify linearity, it will assume linearity.
    *templ*
      Boolean, whether data to extract is of a template spectrum, which means it includes avg flux, flux variance, min and max flux at each wavelength.
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
        print('File name(s) in invalid format.')
        return

    # 2. Initialize array to store spectra
    specData = [None] * len(specFiles)

    # 3. Loop through each file name:
    for spFileIdx,spFile in enumerate(specFiles):
        if spFile is None: continue
        # 3.1 Determine the type of file it is
        isFits = False
        ext = spFile[-4:].lower()
        if ext == 'fits' or ext == '.fit':
            isFits = True
        
        # 3.2. Get data from file
        if isFits:
            isSDSS = False
            isLAMOST = False
            try:
                # Determine table index to extract the data
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

        # Assume ascii file otherwise
        else:
            try:
                aData = ascii.read(spFile)
                specData[spFileIdx] = [aData[0].tonumpy(), aData[1].tonumpy()]
                if len(aData) >= 3 and errors:
                    specData[spFileIdx].append(aData[2].tonumpy())
            except IOError:
                print('Could not open ' + str(spFile) + '.')
                continue

        # 3.3. Check if data in fits file is linear
        if isFits:
            KEY_TYPE = ['CTYPE1']
            setType  = set(KEY_TYPE).intersection(set(fitsHeader.keys()))
            if len(setType) == 0:
                if verbose:
                    print('Data in ' + spFile + ' assumed to be linear.')
                isLinear = True
            else:
                valType = fitsHeader[setType.pop()]
                if valType.strip().upper() == 'LINEAR':
                    isLinear = True
                else:
                    isLinear = False
            if linear and not isLinear:
                if verbose:
                    print('Data in ' + spFile + ' is not linear.')
                return

        # 3.4. Get wl, flux & error data from fits file
        #      (returns wl in pos. 0, flux in pos. 1, error values in pos. 2)
        #      (If template spec: min flux in pos. 3, max flux in pos. 4)
        if isFits:
            specData[spFileIdx] = __get_spec(fitsData, fitsHeader, spFile, errors, \
                                             templ=templ, verb=verbose)
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

        # 3.5. Convert units in wl-axis from Angstrom into microns if desired
        if atomicron:
            if specData[spFileIdx][0][-1] > 8000:
                specData[spFileIdx][0] = specData[spFileIdx][0] / 10000

        # 3.6. Set negative flux values equal to zero (next step sets them to nans)
        if negtonan:
            negIdx = np.where(specData[spFileIdx][1] < 0)
            if len(negIdx[0]) > 0:
                specData[spFileIdx][1][negIdx] = 0
                if verbose:
                    print('%i negative data points found in %s.' \
                            % (len(negIdx[0]), spFile))

        # 3.7. Set zero flux values as nans (do this always)
        zeros = np.where(specData[spFileIdx][1] == 0)
        if len(zeros[0]) > 0:
            specData[spFileIdx][1][zeros] = np.nan


    # 4. Plot the spectra if desired
    if plot:
        plot_spec(specData, ploterrors=True)

    # 5. Clear up memory
    fitsData  = ''

    if header:
        return specData, fitsHeader
    else:
        return specData

def snr(spec, rng=None):
    '''
    (by Alejandro N |uacute| |ntilde| ez)

    Calculate signal-to-noise in a spectrum.

    *spec*
      Spectrum as a Python list with wavelength in position 0, flux in position 1, and error values in position 2. It can also be a list of spectra. If no errors available, then it calculates SNR based on this: http://www.stecf.org/software/ASTROsoft/DER_SNR/der_snr.py. 
    *rng*
      list, indicating in wavelength space the range of interest. If None, it computes signal-to-noise for the whole spectrum.
    '''

    # Convert spec into a list type if it is only one spectrum
    if len(spec[0]) > 3:
        spec = [spec,]
    
    snr = np.array([np.nan] * len(spec))
    for js,s in enumerate(spec):
        i = np.where((s[1] != 0.0) & (np.isfinite(s[1])))[0]
        flux = np.array(s[1][i])
        wl = np.array(s[0][i])
        try:
            e_flux = np.array(s[2][i])
            i = np.where(np.isfinite(e_flux))[0]
            if len(i) > 0:
                errors = True
            else:
                errors = False
        except IndexError:
            errors = False
        
        if errors:
            if rng is None:
                snr[js] = np.median(flux / e_flux)
            else:
                if rng[0] >= rng[1]:
                    print('Wavelength range incorrect.')
                    return
                else:
                    i = np.where((wl > rng[0]) & (wl < rng[1]))[0]
                    if len(i) == 0:
                        print('No flux data within specified range.')
                        return
                    else:
                        snr[js] = np.median(flux[i] / e_flux[i])
        else:
            if rng is None:
                n = len(flux)
                flx = flux.copy()
            else:
                if rng[0] >= rng[1]:
                    print('Wavelength range incorrect.')
                    return
                else:
                    i = np.where((wl > rng[0]) & (wl < rng[1]))[0]
                    n = len(i)
                    flx = flux[i]
            if n < 4:
                print('At least 4 flux data points are needed for this calculation.')
                return
            else:
                signal = np.median(flx)
                noise  = 0.6052697 * np.median(np.abs(2.0 * flx[2:n-2] - flx[0:n-4] - flx[4:n]))
                snr[js] = signal / noise

    return snr

def plot_spec(specData, ploterrors=False):
    '''
    (by Alejandro N |uacute| |ntilde| ez)

    Plot a spectrum. If more than one spectrum is provided simultaneously, it will plot all spectra on top of one another.

    This is a quick and dirty tool to visualize a set of spectra. It is not meant to be a paper-ready format. You can use it, however, as a starting point.

    *specData*
      Spectrum as a Python list with wavelength in position 0, flux in position 1, and (optional) error values in position 2. More than one spectrum can be provided simultaneously, in which case *specData* shall be a list of lists.
    *ploterrors*
      Boolean, whether to include flux error bars when available. This will work only if all spectra have error values.

    '''

    # Check that there is data to plot
    allNone = True
    for spData in specData:
        if spData is not None:
            allNone = False
            break
    if allNone:
        return

    # Fix specData list dimensions when necessary
    if len(specData) == 2 or len(specData) == 3:
        if len(specData[0]) > 3:
            specData = [specData]

    # Initialize figure
    plt.close()
    fig = plt.figure(1)
    fig.clf()

    # Set plot titles
    TITLE   = 'SPECTRAL DATA'
    X_LABEL = 'Wavelength'
    Y_LABEL = 'Flux'

    # Initialize plot within figure
    subPlot = fig.add_subplot(1,1,1)
    subPlot.set_title(TITLE)
    subPlot.set_xlabel(X_LABEL)
    subPlot.set_ylabel(Y_LABEL)

    # Check if all spectra have error values
    errorsOK = True
    for spData in specData:
        if len(spData) != 3:
            errorsOK = False

    # Plot spectra
    for spData in specData:
        if spData is not None:
            if errorsOK and ploterrors:
                subPlot.errorbar(spData[0], spData[1], spData[2], \
                          capsize=2, drawstyle='steps-mid')
            else:
                subPlot.plot(spData[0], spData[1], drawstyle='steps-mid')

    return fig

def edit_header(fitsfiles, keyword, val, hdu=0):
    """
    Edit a card on the fits file header using the parameters provided.

    Args:
    ----------
    fitsfile - String, the full path of the fits file; if only a filename is provided, it will look for the file in the current directory. It can also be a python list of names.
    keyword - String, the name of the keyword to edit.
    val - String, the value that the keyword will have.
    hdu - Int, the index of the hdu to be edited.
    Returns:
    ----------
    - None.
    """

    import datetime

    # Convert fitsfiles into a list type if it is only one file name
    if isinstance(fitsfiles, str):
        fitsfiles = [fitsfiles,]

    for fitsfl in fitsfiles:

        # Read fits file data
        FitsHDU = pf.open(fitsfl, 'update')
        try:
            tmp = FitsHDU[hdu].data.shape
        except IndexError:
            print('hdu index does not exist for ' + fitsfl)
            print('Skipping this file.')
            continue
        try:
            tmp = FitsHDU[hdu].header[keyword]
        except KeyError:
            print('Keyword does not exist for ' + fitsfl)
            print('Skipping this file.')
            continue
        
        # Replace keyword value with new one
        FitsHDU[hdu].header[keyword] = val
        today = datetime.datetime.now().strftime('%Y-%m-%d')
        origcomment = FitsHDU[hdu].header.comments[keyword]
        FitsHDU[hdu].header.comments[keyword] = origcomment + ' ---Updated on ' + today + ' by antools.py.'

        FitsHDU.flush()

    return 

def crop_fits(fitsfile, xsize, ysize, croploc='center', suffix=None):
    """
    Crop a fits image using the parameters provided. If file has more than one image, it only considers the first one.

    Args:
    ----------
    fitsfile - String, the full path of the fits file; if only a filename is provided, it will look for the file in the current directory.
    xsize - Int, the desired X size (columns) in pixels.
    ysize - Int, the desired Y size (rows) in pixels.
    croploc - ['center'(default), 'upper right', 'upper left', 'lower left', 'lower right'], set location around which to crop image. If 'center', then it crops image centered in the image center. If 'upper right', then it crops image to size [xsize,ysize] anchored in the upper right corner. And so on...
    suffix - String, suffix to add to new fits file. If it is None, then the original fits file is overwritten with the new one.
    Returns:
    ----------
    - the new fits HDU, including the original header information.
    - It also saves a copy of the newly created fits file in the same folder as the original file, with an added suffix to its name, if "suffix" is specified.
    """

    import os

    # Get file path, if provided, and filename
    filepath = fitsfile.rsplit('/',1)[0]
    if filepath == fitsfile:
        filepath = ''
        filename = fitsfile.rsplit('.',1)[0]
    else:
        filepath = filepath + '/'
        filename = fitsfile.rsplit('/',1)[1].rsplit('.',1)[0]

    # Read fits file data
    FitsHDU = pf.open(fitsfile)
    Im = FitsHDU[0].data
    FitsHeader = FitsHDU[0].header
    xsizeorig = FitsHeader['NAXIS1']
    ysizeorig = FitsHeader['NAXIS2']

    # Determine pixel limits for cropping
    if croploc == 'center':
        center = [int(xsizeorig/2), int(ysizeorig/2)]
        xstart = center[0] - int(xsize/2) + 1
        xstop = center[0] + int(xsize/2) + 1
        ystart = center[1] - int(ysize/2)
        ystop = center[1] + int(ysize/2)
    elif croploc == 'upper right':
        xstart = xsizeorig - xsize + 1
        xstop = xsizeorig + 1
        ystart = ysizeorig - ysize
        ystop = ysizeorig + 1
    elif croploc == 'upper left':
        xstart = 1
        xstop = xsize + 1
        ystart = ysizeorig - ysize + 1
        ystop = ysizeorig + 1
    elif croploc == 'lower left':
        xstart = 1
        xstop = xsize + 1
        ystart = 1
        ystop = ysize + 1
    elif croploc == 'lower right':
        xstart = xsizeorig - xsize + 1
        xstop = xsizeorig + 1
        ystart = 1
        ystop = ysize + 1
    else:
        print('croploc not recognized.')
        return None

    # Check that cropping dimensions are OK
    if any((xstart<1, xstop<1, ystart<1,ystop<1)):
        print('xsize/ysize dimensions are too large.')
        return None
    if any((xstart>xsizeorig+1, xstop>xsizeorig+1)):
        print('xsize dimensions are too large.')
        return None
    if any((ystart>ysizeorig+1, ystop>ysizeorig+1)):
        print('ysize dimensions are too large.')
        return None

    #Crop the image
    Im = Im[ystart:ystop, xstart-1:xstop]
    FitsHDU[0].data=Im

    #Write it to a new file
    if suffix is not None:
        suffix = '_' + suffix
    else:
        suffix = ''
    OutFile = filepath + filename + suffix + '.fits'
    if os.path.exists(OutFile) : os.remove(OutFile)
    FitsHDU.writeto(OutFile)

    return FitsHDU


def __create_waxis(fitsHeader, lenData, fileName, verb=True):
    # Function used by read_spec only
    # (by Alejo)
    # Generates a wavelength (wl) axis using header data from fits file.

    # Define key names in
    KEY_MIN  = ['COEFF0','CRVAL1']         # Min wl
    KEY_DELT = ['COEFF1','CDELT1','CD1_1'] # Delta of wl
    KEY_OFF  = ['LTV1']            # Offset in wl to subsection start

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
            print('Could not re-create wavelength axis for ' + fileName + '.')

    return wAxis


def __get_spec(fitsData, fitsHeader, fileName, errorVals, templ=False, verb=True):
    # Function used by read_spec only
    # (by Alejo)
    # Interprets spectral data from fits file.
    # Returns wavelength (wl) data in pos. 0, flux data in pos. 1, and if requested, error values in pos. 2.
    # If templ, also returns min flux in pos. 3 and max flux in pos. 4
    
    if templ:
        validData = [None] * 5
    elif errorVals:
        validData = [None] * 3
    else:
        validData = [None] * 2

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
        if templ:
            # 0-wl, 1-avg flux, 2-flux variance, 3-min flux, 4-max flux
            waveIdx = 0
            fluxIdx = 1
            sigmaIdx = 2
            minIdx = 3
            maxIdx = 4
        else:
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
        elif len(fitsData[0]) == 2:
        # Data comes in a xxxx by 2 matrix (ascii origin)
            tmpWave = []
            tmpFlux = []
            for pair in fitsData:
                tmpWave.append(pair[0])
                tmpFlux.append(pair[1])
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
            print('Unable to interpret data in ' + fileName + '.')
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

    # Fetch template data when relevant
    if templ:
        validData[3] = fitsData[minIdx]
        # validData[4] = fitsData[maxIdx]
    
    # Check ascending order of spectrum using wavelength axis
    if validData[0] is not None:
        if validData[0][0] > validData[0][-1]:
            for i in range(len(validData)):
                if validData[i] is not None:
                    validData[i] = validData[i][::-1]

    return validData


def equivalent_width(spec, xmin, xmax, exclude_min, exclude_max, n, fldr, name=None, errors=True, head=None, normalize=True, band='Halpha', fitted=True, multi=False):
    """Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam

    Args:
    ----------
    spec - String, fits filename
    xmin,xmax - Integers, the specified interval in wavelength space, which defines the region of interest
    excludemin, excludemax - Integers, the specified interval in wavelength space of the spectral feature, which binds the edges of the spectral feature itself
    n - Integer, the number of times the EqW measurement is repeated in the MCMC step
    fldr - String, location where output figure is desired
    name - String, if not None, it uses it to label the star
    errors - Boolean, whether to perform the MCMC routine to calculate 1-sigma errors for EqW
    head - Fits header of fits file with data
    normalize - Boolean, whether to normalize the flux values
    fitted - Boolean, whether EqW is calculated using fitted model; if False, it uses instead the data points
    multi - Boolean, whether two spectral features are fitted simultaneously instead of one (NII or SII features)
    Returns:
    -------
    - the mean and standard deviation of the equivalent width measured n times
    - A figure with a plot of the full spectrum; a plot of the spectral feature fit, with the Voigt profile line fit (blue), the pseudo-continuum (orange), and the approximated rectangle (green); and a histogram with the MCMC results of the EqW distribution
    """

    import pyspeckit as p
    import matplotlib.pyplot as plt, matplotlib.mlab as mlab
    import numpy as np
    from scipy.stats import norm
    import astropy.io.fits as pyfits

    GRAY = '#999999'
    BLORDER = 1 # order for baseline fitting

    # Set band parameters
    if band == 'Halpha':
        normwl = 6555.
        bandloc = 6563.
    elif band == 'NII':
        if multi:
            normwl = 6545.
            bandloc = 6548.
            bandloc2 = 6584.
        else:
            normwl = 5750.
            bandloc = 5755.
    elif band == 'SII':
        normwl = 6711.
        bandloc = 6717.
        bandloc2 = 6731.

    # Get data from fits file
    if spec.endswith('.fits') or spec.endswith('.fit'):
        if head is None:
            data, head = read_spec(spec, header=True)
            data = data[0] # This is needed bc read_spec() returns a list, even if it's just one file
        else:
            data = read_spec(spec)[0]
    else:
        tb = read_table(spec, delimiter=' ', ds=0)
        data = np.vstack([tb.columns[0], tb.columns[1]])
    
    if data is None: return None

    # Get object name
    # (Don't get name from fits header, bc sometimes it's wrong)
    if name is not None:
        objnm = name
    else:
        objnm = spec.split('/')[-1].split('.')[-2] # This is just the filename

    # Set up figure
    plt.rc('font', size=8)
    fig = plt.figure(1, figsize=(6*1.2,6))
    plt.subplots_adjust(top=0.96, bottom=0.07, right=0.98, left=0.08)

    if multi:
        numcols = 3
    else:
        numcols = 2
    ax1 = plt.subplot2grid((2,numcols), (0,0), 1, numcols) # To plot the fit
    ax2 = plt.subplot2grid((2,numcols), (1,0), 1, 1) # To plot the full spectrum
    ax3 = plt.subplot2grid((2,numcols), (1,1), 1, 1) # To plot the histogram
    if multi:
        ax4 = plt.subplot2grid((2,numcols), (1,2), 1, 1) # To plot the second histogram

    # Plot plain spectrum
    tmpmin = data[0][1]
    if tmpmin < 4500:
        tmpmin = 4500.
    tmpmax = data[0][-2]
    irange = np.where((data[0]>=tmpmin) & (data[0]<=tmpmax))[0]
    mean = np.nanmean(data[1][irange])
    std = np.nanstd(data[1][irange])
    iclip = np.where((data[1][irange] > mean - 3*std) & \
                        (data[1][irange] < mean + 3*std))[0] # Clip 3sigma flux outliers
    inorm = np.where(data[0] >= normwl)[0][0] # Always normalize spectrum flux in this figure
    normval = data[1][inorm]
    ax2.plot(data[0][irange][iclip], data[1][irange][iclip] / normval, \
             drawstyle='steps-mid', linewidth=0.8)
    ax2.set_xlim(xmin=tmpmin, xmax=tmpmax)
    ax2.set_xlabel(r'Wavelength $(\AA)$')
    ax2.set_ylabel(r'Flux / Flux(' + format(int(normwl)) + ' $\AA$)')
    ax2.axvline(x=bandloc, linestyle='--', color=GRAY, linewidth=0.8)
    if band == 'NII':
        ax2.axvline(x=bandloc2, linestyle='--', color=GRAY, linewidth=0.8)

    # Normalize flux values (bring them close to unity to make aid the fitting routine)
    if normalize:
        data[1] = data[1] / normval * 10
        if len(data) == 3:
            data[2] = data[2] / normval * 10
    
    # Load spectrum onto PySpecKit class
    if len(data) < 3:
        # Only wavelength and flux arrays
        sp = p.Spectrum(data=data[1], xarr=data[0], header=head, \
                        xarrkwargs={'unit':'angstroms'}) 
    else:
        # Only wavelength and flux arrays
        if np.all(np.isnan(data[2])):
            sp = p.Spectrum(data=data[1], xarr=data[0], header=head, \
                            xarrkwargs={'unit':'angstroms'})            
        # Wavelength, flux, and e_flux arrays
        else: 
            sp = p.Spectrum(data=data[1], xarr=data[0], error=data[2], header=head, \
                            xarrkwargs={'unit':'angstroms'})
    sp.xarr.xtype = 'wavelength'
    if name is not None or sp.specname == '':
        sp.specname = objnm
    if normalize:
        sp.unit = 'Normalized flux'

    # Set up plotter and fit baseline
    sp.plotter(axis=ax1, clear=False, xmin=xmin, xmax=xmax, ymin=0, \
               errstyle='bars', color='grey')
    sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], \
                subtract=False, reset_selection=False, highlight_fitregion=False, \
                order=BLORDER)
    sp.baseline.annotate(loc='upper right', fontsize=8)

    # Fit Voigt profile to spectral feature
    if multi:
        if band == 'NII':
            tmpguess = [20,6548.,0.8,0.5,50,6584.,0.8,0.5] # amp, delX, sigma, gamma
        elif band == 'SII':
            tmpguess = [50,6717,0.8,0.5,50,6731.,0.8,0.5] # amp, delX, sigma, gamma
        sp.specfit(fittype='voigt', color='blue', loc='center right', multifit=multi, \
                   guesses=tmpguess)
        # Calculate equivalent width using the fit above
        ew = sp.specfit.EQW(plot=True, plotcolor='g', fitted=fitted, components=multi, \
                            annotate=True, loc='lower left', xmin=None, xmax=None)
    else:
        tmpguess =  [1., bandloc, 1., 1.] # None
        sp.specfit(fittype='voigt', color='blue', loc='center right', \
                   guesses=tmpguess)
        # Calculate equivalent width using the fit above
        ew = sp.specfit.EQW(plot=True, plotcolor='g', fitted=fitted, xmin=None, xmax=None)
        sp.specfit.annotate(loc='center right', fontsize=8)
        txt = 'EqW = ' + format(ew,'.2f') + r' $\AA$'
        ax1.text(0.86,0.02, txt, transform=ax1.transAxes)
    
    # Beautify plot and save it
    sp.specfit.fitleg.set_bbox_to_anchor((1,0.3),transform=sp.plotter.axis.transAxes)
    sp.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
    ylbl = sp.plotter.axis.get_ylabel()
    sp.plotter.axis.axvline(x=exclude_min, linestyle=':', color=GRAY)
    sp.plotter.axis.axvline(x=exclude_max, linestyle=':', color=GRAY)
    if multi:
        sp.plotter.axis.axvline(x=bandloc, linestyle='--', color='green')
        sp.plotter.axis.axvline(x=bandloc2, linestyle='--', color='green')
        tmplgd = sp.plotter.axis.get_legend()
        tmplgd.set_bbox_to_anchor((0.98,0.3), transform=sp.plotter.axis.transAxes)
        tmplgd.set_frame_on(True)

    # Print figure to allow user to determine if fit is acceptable
    plt.savefig(fldr + objnm + '_EqWfit.pdf') 

    # Do MCMC using the original result as starting point for param values
    if errors:
        sp2 = sp.copy()
        EQWs = []
        for w in range(n):
            print(w)
            if np.all(np.isnan(data[2])):
                # Get error from noise in continuum
                icont = np.where(((sp.xarr.value >= xmin) & (sp.xarr.value < exclude_min)) |  \
                                ((sp.xarr.value > exclude_max) & (sp.xarr.value <= xmax)))[0]
                tmpcont = sp.data[icont]
                tmperr = np.std(tmpcont)
            else:
                # Get error from flux uncertainties
                tmperr = sp.error

            sp2.data = sp.data + np.random.randn(sp.data.size) * tmperr
            sp2.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max], \
                         subtract=False, reset_selection=False, order=BLORDER)
            if multi:
                sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values, \
                            multifit=multi)
                dist = sp2.specfit.EQW(fitted=fitted, components=multi, \
                                       annotate=True, xmin=None, xmax=None)
            else:
                sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values)
                dist = sp2.specfit.EQW(fitted=fitted, \
                                       annotate=True, xmin=None, xmax=None)
            EQWs.append(dist)
        EQWs = np.array(EQWs)

        # Calculate stats of MCMC array and make histogram with results
        if multi:
            mu, sigma = norm.fit(EQWs[:,0])
            mu2, sigma2 = norm.fit(EQWs[:,1])

            n,bins,ptchs = ax3.hist(EQWs[:,0], 10, normed=True, facecolor='green', \
                                   histtype='stepfilled')
            n,bins2,ptchs = ax4.hist(EQWs[:,1], 10, normed=True, facecolor='green', \
                                    histtype='stepfilled')
        else:
            mu, sigma = norm.fit(EQWs)
            n,bins,ptchs = ax3.hist(EQWs, 10, normed=True, facecolor='green', \
                                   histtype='stepfilled')

        # Beautify histogram plot
        y = mlab.normpdf(bins, mu, sigma)
        ax3.plot(bins,y,'r--',linewidth=2)
        ax3.grid(True)
        ax3.set_ylabel('Count')
        ax3.set_xlabel(r'EQW ($\AA$)')
        txt = r'$\mu=$' + format(mu,'.3f') + r', $\sigma=$' + format(sigma,'.3f')
        ax3.text(0.02,0.94, txt, transform=ax3.transAxes, fontsize=8, color='white', \
                 bbox=dict(facecolor='green', ec='none', pad=0.3, boxstyle='round'))
        if multi:
            y = mlab.normpdf(bins2, mu2, sigma2)
            ax4.plot(bins2,y,'r--',linewidth=2)
            ax4.grid(True)
            ax4.set_ylabel('Count')
            ax4.set_xlabel(r'EqW ($\AA$)')
            txt = r'$\mu=$' + format(mu2,'.3f') + r', $\sigma=$' + format(sigma2,'.3f')
            ax4.text(0.02,0.94, txt, transform=ax4.transAxes, fontsize=8, color='white', \
                     bbox=dict(facecolor='green', ec='none', pad=0.3, boxstyle='round'))

        plt.savefig(fldr + objnm + '_EqWfit.pdf')

        if multi:
            return np.array([mu, sigma, mu2, sigma2])
        else:
            return np.array([mu, sigma])
    else:
        plt.savefig(fldr + objnm + '_EqWfit.pdf')
        if multi:
            return np.array(ew, 0.)
        else:
            return np.array([ew, 0.])
