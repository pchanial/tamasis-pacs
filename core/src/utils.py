import config
import kapteyn
import numpy
import os
import pyfits
import re
import tamasisfortran as tmf
import time
from matplotlib import pyplot
from numpyutils import _my_isscalar

__all__ = [ 'create_fitsheader', 'MaskPolicy', 'mean_degrees', 'minmax_degrees', 'pointing', 'plot_scan', 'airy_disk', 'aperture_circular', 'distance', 'gaussian', 'phasemask_fourquadrant' ]

pointing = numpy.dtype([('time', numpy.float64), ('ra', numpy.float64), ('dec', numpy.float64), ('pa', numpy.float64), ('flag', numpy.int64)])

def create_fitsheader(array, extname=None, crval=(0.,0.), crpix=None, ctype=('RA---TAN','DEC--TAN'), cunit='deg', cd=None, cdelt=None, naxis=None):
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
            raise ValueError('An array argument or naxis keyword should be specified.')
        typename = 'float64'
    else:
        if not isinstance(array, numpy.ndarray):
            raise TypeError('The input is not an ndarray.')
        naxis = tuple(reversed(array.shape))
        typename = array.dtype.name
    
    numaxis = len(naxis)
    if extname is None:
        card = pyfits.createCard('simple', True)
    else:
        card = pyfits.createCard('xtension', 'IMAGE', 'Image extension')
    header = pyfits.Header([card])
    header.update('bitpix', pyfits.PrimaryHDU.ImgCode[typename], 'array data type')
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
        crpix = (numpy.array(naxis) + 1) / 2.
    else:
        crpix = numpy.asarray(crpix, dtype=numpy.float64)
    if crpix.size != 2:
        raise ValueError('CRPIX does not have two elements.')

    ctype = numpy.asarray(ctype, dtype=numpy.string_)
    if ctype.size != 2:
        raise ValueError('CTYPE does not have two elements.')

    if _my_isscalar(cunit):
        cunit = (cunit, cunit)
    cunit = numpy.asarray(cunit, dtype=numpy.string_)
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


def mean_degrees(array):
    """
    Returns the mean value of an array of values in degrees, by taking into 
    account the discrepancy at 0 degrees
    """
    return tmf.mean_degrees(numpy.asarray(array, dtype=config.get_default_dtype_float()).ravel())


#-------------------------------------------------------------------------------


def minmax_degrees(array):
    """
    Returns the minimum and maximum value of an array of values in degrees, 
    by taking into account the discrepancy at 0 degrees
    """
    return tmf.minmax_degrees(numpy.asarray(array, dtype=config.get_default_dtype_float()).ravel())


#-------------------------------------------------------------------------------


class MaskPolicy(object):
    def __init__(self, flags, values, description=None):
        self.description = description
        if _my_isscalar(flags):
            if isinstance(flags, str):
                flags = flags.split(',')
            else:
                flags = (flags,)
        self.flags = tuple(flags)
        if _my_isscalar(values):
            values = (values,)
        if len(flags) != len(values):
            raise ValueError('The number of policy flags is different from the number of policies.')

        conversion_policy = { 'keep':0, 'mask':1, 'remove':2 }
        self._array = numpy.ndarray(len(flags), dtype='int')
    
        for i, (flag, value) in enumerate(zip(flags, values)):
            try:
                self._array[i] = conversion_policy[str(value).lower()]
            except KeyError:
                raise KeyError("Invalid policy "+flag+"='" + value + "'. Valid policies are 'keep', 'mask' or 'remove'.")

    def __array__(self):
        return numpy.array(self._array, dtype='int32')

    def __str__(self):
        str = self.description + ': ' if self.description is not None else ''
        conversion = { 0:'keep', 1:'mask', 2:'remove' }
        str_ = []
        for i, flag in enumerate(self.flags):
            str_.append(conversion[self._array[i]] + " '" + flag + "'")
        str += ', '.join(str_)
        return str


#-------------------------------------------------------------------------------
        

def plot_scan(ra, dec, title=None, num=None, figsize=None, dpi=None):
    crval = [mean_degrees(ra), numpy.mean(dec)]
    ra_min,  ra_max  = minmax_degrees(ra)
    dec_min, dec_max = numpy.min(dec), numpy.max(dec)
    cdelt = numpy.max((ra_max-ra_min)/1000., (dec_max-dec_min)/1000.)
    header = create_fitsheader(None, naxis=[1,1], cdelt=cdelt, crval=crval, crpix=[1,1])
    proj = kapteyn.wcs.Projection(header)
    coords = numpy.array([ra, dec]).T
    xy = proj.topixel(coords)
    xmin = int(numpy.round(numpy.min(xy[:,0])))
    xmax = int(numpy.round(numpy.max(xy[:,0])))
    ymin = int(numpy.round(numpy.min(xy[:,1])))
    ymax = int(numpy.round(numpy.max(xy[:,1])))
    xmargin = int(numpy.ceil((xmax - xmin + 1) / 10.))
    ymargin = int(numpy.ceil((ymax - ymin + 1) / 10.))
    xmin -= xmargin
    xmax += xmargin
    ymin -= ymargin
    ymax += ymargin
    naxis = (xmax - xmin + 1, ymax - ymin + 1)
    crpix = (-xmin+2, -ymin+2)
    header = create_fitsheader(None, naxis=naxis, cdelt=cdelt, crval=crval, crpix=crpix)
    fitsobj = kapteyn.maputils.FITSimage(externalheader=header)
    frame = pyplot.gcf().add_axes([0.1,0.1,0.8,0.8])
    if title is not None:
        frame.set_title(title)
    image = fitsobj.Annotatedimage(frame, blankcolor='w')
    grat = image.Graticule()
    image.plot()
    image.interact_toolbarinfo()
    image.interact_writepos()
    pyplot.show()
    x, y = image.topixel(ra, dec)
    p = pyplot.plot(x, y, linewidth=2)
    pyplot.plot(x[0], y[0], 'o', color = p[0]._color)


#-------------------------------------------------------------------------------


def airy_disk(shape, fwhm, origin=None, resolution=1., dtype=None):
    d  = distance(shape, origin=origin, resolution=resolution, dtype=dtype)
    d *= 1.61633 / (fwhm/2.)
    d  = (2 * jn(1,d)/d)**2
    d /= numpy.sum(d)
    return d


#-------------------------------------------------------------------------------


def aperture_circular(shape, diameter, origin=None, resolution=1., dtype=None):
    array = distance(shape, origin=origin, resolution=resolution, dtype=dtype)
    m = array > diameter / 2.
    array[ m] = 0
    array[~m] = 1
    return array


#-------------------------------------------------------------------------------


def distance(shape, origin=None, resolution=1., dtype=None):
    """
    Returns an array whose values are the distances to a given origin
    
    Parameters
    ----------
    shape : tuple of integer
        dimensions of the output array. For a 2d array, the first integer
        is for the Y-axis and the second one for the X-axis.
    origin : array-like
        coordinates of the origin, for which the output array value is
        zero. Default value is the array center
    resolution : inter-pixel distance
    dtype : type of the output array

    Example
    -------
    nx, ny = 3, 3
    print distance((ny,nx))
    [[ 1.41421356  1.          1.41421356]
    [ 1.          0.          1.        ]
    [ 1.41421356  1.          1.41421356]]
    """
    if type(shape) is not tuple and type(shape) is not list:
        shape = (shape,)
    if origin is None:
        origin = (numpy.array(shape) - 1) / 2.
    dim = []
    for length, c in zip(reversed(shape), reversed(origin)):
        dim.append((numpy.arange(length) - c) * resolution)
    return Map(numpy.sqrt(numpy.sum(numpy.square(numpy.meshgrid(*dim)), axis=0)), dtype=dtype)
    

#-------------------------------------------------------------------------------


def gaussian(shape, fwhm, origin=None, resolution=1., dtype=None):
    sigma = fwhm / numpy.sqrt(2*numpy.log(2))
    d = distance(shape, origin=origin, resolution=resolution, dtype=dtype)
    d = numpy.exp(-d**2/(2*sigma**2)) / (2*numpy.pi*sigma**2)
    return d


#-------------------------------------------------------------------------------


def phasemask_fourquadrant(shape, phase=-1):
    array = Map.ones(shape, dtype=get_default_dtype_complex())
    array[0:shape[0]//2,shape[1]//2:] = phase
    array[shape[0]//2:,0:shape[1]//2] = phase
    return array


