import config
import kapteyn
import numpy
import os
import pyfits
import re
import tamasisfortran as tmf
import time
from datatypes import Map, create_fitsheader
from matplotlib import pyplot
from numpyutils import _my_isscalar
from unit import Quantity

__all__ = [ 'ds9', 'hs', 'mean_degrees', 'minmax_degrees', 'plot_scan', 'airy_disk', 'aperture_circular', 'distance', 'gaussian', 'phasemask_fourquadrant' ]


#-------------------------------------------------------------------------------


def ds9(array, origin=None, **keywords):
    array = numpy.asanyarray(array)
    array = Map(array, dtype=array.dtype, copy=False)
    array.ds9(origin=origin, **keywords)
    

#-------------------------------------------------------------------------------


def hs(arg):
    """
    Display the attributes of an object (except methods and those starting
    with an underscore) or an ndarray with composite dtype alongside the
    names of the records. The display is truncated so that each name fits
    in one line.
    """
    import inspect
    if isinstance(arg, numpy.ndarray):
        names = arg.dtype.names
        if names is None:
            print arg
            return
        print str(arg.size) + ' element' + ('s' if arg.size > 1 else '')
    else:
        members = inspect.getmembers(arg, lambda x: not inspect.ismethod(x) and not inspect.isbuiltin(x))
        members = filter(lambda x: x[0][0] != '_', members)
        names = map(lambda x:x[0], members)

    length = numpy.max(map(len, names))
    lnames = numpy.array([names[i].ljust(length)+': ' for i in range(len(names))])
    for name, lname in zip(names, lnames):
        value = str(getattr(arg, name))[0:72-length-2]
        if len(value) == 72-length-2:
            value = value[0:-3] + '...'
        print lname+value

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


def plot_scan(input, map=None, mask=None, title=None, new_figure=True, color='magenta', linewidth=2):
    if hasattr(input, 'pointing'):
        ra   = input.pointing.ra
        dec  = input.pointing.dec
        if mask is None:
            mask = input.pointing.removed
    elif hasattr(input, 'ra') and hasattr(input, 'dec'):
        ra  = input.ra
        dec = input.dec
        if mask is None:
            if hasattr(input, 'mask'):
                mask = input.mask
            elif hasattr(input, 'removed'):
                mask = input.removed
    else:
        ra  = input[0]
        dec = input[1]

    if map is not None and map.header is not None and map.header.has_key('CRVAL1'):
        image = map.imshow(title=title)
        x, y = image.topixel(ra, dec)
        x[mask] = numpy.nan
        p = pyplot.plot(x, y, color=color, linewidth=linewidth)
        pyplot.plot(x[0], y[0], 'o', color=p[0]._color)
        return image

    crval = [mean_degrees(ra), numpy.mean(dec)]
    ra_min,  ra_max  = minmax_degrees(ra)
    dec_min, dec_max = numpy.min(dec), numpy.max(dec)
    cdelt = numpy.max((ra_max-ra_min)/1000., (dec_max-dec_min)/1000.)
    header = create_fitsheader(None, naxis=[1,1], cdelt=cdelt, crval=crval, crpix=[1,1])
    proj = kapteyn.wcs.Projection(header)
    coords = numpy.array([ra, dec]).T
    xy = proj.topixel(coords)
    xy[mask,:] = numpy.nan
    xmin = int(numpy.round(numpy.nanmin(xy[:,0])))
    xmax = int(numpy.round(numpy.nanmax(xy[:,0])))
    ymin = int(numpy.round(numpy.nanmin(xy[:,1])))
    ymax = int(numpy.round(numpy.nanmax(xy[:,1])))
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
    if new_figure:
        fig = pyplot.figure()
        frame = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    else:
        frame = pyplot.gca()
    if title is not None:
        frame.set_title(title)
    image = fitsobj.Annotatedimage(frame, blankcolor='w')
    grat = image.Graticule()
    image.plot()
    image.interact_toolbarinfo()
    image.interact_writepos()
    pyplot.show()
    x, y = image.topixel(ra, dec)
    p = pyplot.plot(x, y, color=color, linewidth=linewidth)
    pyplot.plot(x[0], y[0], 'o', color = p[0]._color)
    return image


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
    origin : array-like of two elements (x0, y0)
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
    for length, c in zip(reversed(shape), origin):
        dim.append((numpy.arange(length) - c) * resolution)
    return Map(numpy.sqrt(numpy.sum(numpy.square(numpy.meshgrid(*dim)), axis=0)), dtype=dtype)
    

#-------------------------------------------------------------------------------


def gaussian(shape, fwhm, origin=None, resolution=1., dtype=None):
    sigma = fwhm / numpy.sqrt(8*numpy.log(2))
    d = distance(shape, origin=origin, resolution=resolution, dtype=dtype)
    d = numpy.exp(-d**2/(2*sigma**2)) / (2*numpy.pi*sigma**2)
    return d


#-------------------------------------------------------------------------------


def phasemask_fourquadrant(shape, phase=-1):
    array = Map.ones(shape, dtype=config.get_default_dtype_complex())
    array[0:shape[0]//2,shape[1]//2:] = phase
    array[shape[0]//2:,0:shape[1]//2] = phase
    return array


