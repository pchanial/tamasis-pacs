import numpy
import os
import pyfits
import re
import tamasisfortran as tmf
import time

__all__ = [ 'any_neq', 'create_fitsheader' ]

def create_fitsheader(array, crval=(0.,0.), crpix=None, ctype=('RA---TAN','DEC--TAN'), cunit='deg', cd=None, cdelt=None):
    """
    Return a FITS header

    Parameters
    ----------
    arra : array_like
        An array from which the dimensions will be extracted. Note that
        by FITS convention, the dimension along X is the second value 
        of the array shape and that the dimension along the Y axis is 
        the first one.
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

    Examples
    --------
    >>> map = Map.ones((10,100), unit='Jy/pixel')
    >>> map.header = create_fitsheader(map, cd=[[-1,0],[0,1]])
    """

    if not isinstance(array, numpy.ndarray):
        raise TypeError('The input is not an ndarray.')

    axisn = tuple(reversed(array.shape))

    naxis = len(axisn)
    header = pyfits.Header()
    header.update('simple', True)
    header.update('bitpix', pyfits.PrimaryHDU.ImgCode[array.dtype.name])
    header.update('extend', True)
    header.update('naxis', naxis)
    for dim in range(naxis):
        header.update('naxis'+str(dim+1), axisn[naxis-dim-1])

    if cd is not None:
        cd = numpy.asarray(cd, dtype=numpy.float64)
        if cd.shape != (2,2):
            raise ValueError('The CD matrix is not a 2x2 matrix.')
    else:
        if cdelt is None:
            return header
        if _my_isscalar(cdelt):
            cdelt = (-cdelt, cdelt)
        cd = numpy.array(((cdelt[0], 0), (0,cdelt[1])))

    crval = numpy.asarray(crval, dtype=numpy.float64)
    if crval.size != 2:
        raise ValueError('CRVAL does not have two elements.')

    if crpix is None:
        crpix = numpy.array(axisn) / 2 + 1
    else:
        crpix = numpy.asarray(crpix, dtype=numpy.float64)
    if crpix.size != 2:
        raise ValueError('CRPIX does not have two elements.')

    ctype = numpy.asarray(ctype, dtype=numpy.string_)
    if ctype.size != 2:
        raise ValueError('CTYPE does not have two elements.')

    cunit = numpy.asarray(cunit, dtype=numpy.string_)
    if _my_isscalar(cunit):
        cunit = (cunit, cunit)
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


# should use numpy's allclose instead
def any_neq(a,b, precision):
    mask = numpy.isnan(a)
    if numpy.any(mask != numpy.isnan(b)):
        return True
    return numpy.any((abs(a-b) > 10.**(-precision) * abs(a))[numpy.isnan(a).__neg__()])
#-------------------------------------------------------------------------------


def get_type(data):
    """Returns input's data type."""
    data_ = numpy.asarray(data)
    type_ = data_.dtype.type.__name__
    if type_[-1] == '_':
        type_ = type_[0:-1]
    if type_ != 'object':
        return type_
    return type(data).__name__

   
#-------------------------------------------------------------------------------


def _my_issctype(dtype):
    """Hack around numpy.issctype bug"""
    return numpy.issctype(dtype) and str(dtype)[0:2] != '|S'
   
   
#-------------------------------------------------------------------------------


def _my_isscalar(data):
    """Hack around numpy.isscalar bug"""
    return numpy.rank(data) == 0 if isinstance(data, numpy.ndarray) else numpy.isscalar(data)
