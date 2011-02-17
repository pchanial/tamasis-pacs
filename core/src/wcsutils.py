import numpy as np
import pyfits
import tamasisfortran as tmf

from . import var

__all__ = [ 
    'angle_lonlat',
    'barycenter_lonlat',
    'create_fitsheader',
    'minmax_degrees',
    'mean_degrees',
]

def angle_lonlat(lon1, lat1, lon2=None, lat2=None):
    """
    Returns the angle between vectors on the celestial sphere in degrees.

    Parameters
    ----------
    lon1, lon2 : array of numbers
        longitude in degrees
    lat1, lat2 : array of numbers
        latitude in degrees

    Example
    -------
    >>> angle_lonlat((ra1,dec1), (ra2,dec2))
    >>> angle_lonlat(ra1, dec1, ra2, dec2)
    """

    if lon2 is None and lat2 is None:
        lon2, lat2 = lat1
        lon1, lat1 = lon1
    lon1 = np.array(lon1, dtype=var.FLOAT_DTYPE, ndmin=1, copy=False).ravel()
    lat1 = np.array(lat1, dtype=var.FLOAT_DTYPE, ndmin=1, copy=False).ravel()
    lon2 = np.array(lon2, dtype=var.FLOAT_DTYPE, ndmin=1, copy=False).ravel()
    lat2 = np.array(lat2, dtype=var.FLOAT_DTYPE, ndmin=1, copy=False).ravel()
    angle = tmf.angle_lonlat(lon1, lat1, lon2, lat2)
    if angle.size == 1:
        angle = float(angle)
    return angle


#-------------------------------------------------------------------------------


def barycenter_lonlat(lon, lat):
    """
    Returns the barycenter of vectors on the celestial sphere.

    Parameters
    ----------
    lon : array of numbers
        longitude in degrees
    lat : array of numbers
        latitude in degrees
    """
    lon = np.asarray(lon, dtype=var.FLOAT_DTYPE).ravel()
    lat = np.asarray(lat, dtype=var.FLOAT_DTYPE).ravel()
    return tmf.barycenter_lonlat(lon, lat)


#-------------------------------------------------------------------------------


def create_fitsheader(array, extname=None, crval=(0.,0.), crpix=None,
                      ctype=('RA---TAN','DEC--TAN'), cunit='deg', cd=None,
                      cdelt=None, naxis=None):
    """
    Return a FITS header

    Parameters
    ----------
    array : array_like
        An array from which the dimensions will be extracted. Note that
        by FITS convention, the dimension along X is the second value 
        of the array shape and that the dimension along the Y axis is 
        the first one. If None is specified, naxis keyword must be set
    extname : None or string
        if a string is specified ('' can be used), the returned header
        type will be an Image HDU (otherwise a Primary HDU)
    crval : 2 element array, optional
        Reference pixel values (FITS convention)
    crpix : 2 element array, optional
        Reference pixel (FITS convention)
    ctype : 2 element string array, optional
        Projection types
    cunit : string or 2 element string array
        Units of the CD matrix (default is degrees/pixel)
    cd : 2 x 2 array
        Astrometry parameters
            CD1_1 CD1_2
            CD2_1 CD2_2
    cdelt : 2 element array
        Physical increment at the reference pixel
    naxis : 2 element array
        (NAXIS1,NAXIS2) tuple, to be specified only if array argument is None

    Examples
    --------
    >>> map = Map.ones((10,100), unit='Jy/pixel')
    >>> map.header = create_fitsheader(map, cd=[[-1,0],[0,1]])
    """

    if array is None:
        if naxis is None:
            raise ValueError('An array argument or naxis keyword should be sp' \
                             'ecified.')
        typename = 'float64'
    else:
        if not isinstance(array, np.ndarray):
            raise TypeError('The input is not an ndarray.')
        naxis = tuple(reversed(array.shape))
        if array.dtype.itemsize == 1:
            typename = 'uint8'
        elif array.dtype.names is not None:
            typename = None
        else:
            typename = array.dtype.name

    if type(naxis) not in (list, tuple):
        naxis = (naxis,)
    numaxis = len(naxis)

    if extname is None:
        card = pyfits.createCard('simple', True)
    else:
        card = pyfits.createCard('xtension', 'IMAGE', 'Image extension')
    header = pyfits.Header([card])
    if typename is not None:
        header.update('bitpix', pyfits.PrimaryHDU.ImgCode[typename],
                      'array data type')
    header.update('naxis', numaxis, 'number of array dimensions')
    for dim in range(numaxis):
        header.update('naxis'+str(dim+1), naxis[dim])
    if extname is None:
        header.update('extend', True)
    else:
        header.update('pcount', 0, 'number of parameters')
        header.update('gcount', 1, 'number of groups')
        header.update('extname', extname)

    if cd is not None:
        cd = np.asarray(cd, dtype=np.float64)
        if cd.shape != (2,2):
            raise ValueError('The CD matrix is not a 2x2 matrix.')
    else:
        if cdelt is None:
            return header
        if not isinstance(cdelt, (list, tuple, np.ndarray)):
            cdelt = (-cdelt, cdelt)
        cd = np.array(((cdelt[0], 0), (0,cdelt[1])))

    crval = np.asarray(crval, dtype=np.float64)
    if crval.size != 2:
        raise ValueError('CRVAL does not have two elements.')

    if crpix is None:
        crpix = (np.array(naxis) + 1) / 2.
    else:
        crpix = np.asarray(crpix, dtype=np.float64)
    if crpix.size != 2:
        raise ValueError('CRPIX does not have two elements.')

    ctype = np.asarray(ctype, dtype=np.string_)
    if ctype.size != 2:
        raise ValueError('CTYPE does not have two elements.')

    if not isinstance(cunit, (list, tuple, np.ndarray)):
        cunit = (cunit, cunit)
    cunit = np.asarray(cunit, dtype=np.string_)
    if cunit.size != 2:
        raise ValueError('CUNIT does not have two elements.')

    header.update('crval1', crval[0])
    header.update('crval2', crval[1])
    header.update('crpix1', crpix[0])
    header.update('crpix2', crpix[1])
    header.update('cd1_1' , cd[0,0])
    header.update('cd2_1' , cd[1,0])
    header.update('cd1_2' , cd[0,1])
    header.update('cd2_2' , cd[1,1])
    header.update('ctype1', ctype[0])
    header.update('ctype2', ctype[1])
    header.update('cunit1', cunit[0])
    header.update('cunit2', cunit[1])

    return header


#-------------------------------------------------------------------------------


def get_cdelt_crota2(header):
    try:
        cd = np.array([[header['CD1_1'], header['CD1_2']],
                       [header['CD2_1'], header['CD2_2']]], float)
    except KeyError:
        if any([not k in header for k in ('CDELT1', 'CDELT2', 'CROTA2')]):
            raise KeyError('Header has no astrometry.')
        return np.array([header['CDELT1'], header['CDELT2']]), header['CROTA2']

    det = np.linalg.det(cd)
    sgn = det / np.abs(det)
    cdelt = np.array([sgn * np.sqrt(cd[0,0]**2 + cd[0,1]**2),
                      np.sqrt(cd[1,1]**2 + cd[1,0]**2)])
    rot = np.arctan2(-cd[1,0], cd[1,1])
    return cdelt, np.rad2deg(rot)
    

#-------------------------------------------------------------------------------


def mean_degrees(array):
    """
    Returns the mean value of an array of values in degrees, by taking into 
    account the discrepancy at 0 degree
    """
    return tmf.mean_degrees(np.asarray(array, dtype=var.FLOAT_DTYPE).ravel())


#-------------------------------------------------------------------------------


def minmax_degrees(array):
    """
    Returns the minimum and maximum value of an array of values in degrees, 
    by taking into account the discrepancy at 0 degree.
    """
    return tmf.minmax_degrees(np.asarray(array, dtype=var.FLOAT_DTYPE).ravel())
