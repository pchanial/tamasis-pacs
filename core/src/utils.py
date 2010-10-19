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

__all__ = [ 'FlatField', 'MaskPolicy', 'Pointing', 'hs', 'mean_degrees', 'minmax_degrees', 'plot_scan', 'airy_disk', 'aperture_circular', 'distance', 'gaussian', 'phasemask_fourquadrant' ]


#-------------------------------------------------------------------------------


def hs(array):
    """
    Display the content of an ndarray with composite dtype alongside the
    names of the records. The display is truncated so that each name fits
    in one line.
    """
    array = numpy.asarray(array)
    names = array.dtype.names
    if names is None:
        print array
        return
    length = numpy.max(map(len, names))
    lnames = numpy.array([names[i].ljust(length)+': ' for i in range(len(names))])
    print str(array.size) + ' element' + ('s' if array.size > 1 else '')

    for name, lname in zip(names, lnames):
        value = str(array[name])[0:72-length-2]
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


class FlatField(object):
    def __init__(self, optical, detector):
        self.optical = Map(optical, origin='upper')
        self.detector = Map(detector, origin='upper')

    @property
    def total(self):
        return Map(self.optical * self.detector, origin='upper')


#-------------------------------------------------------------------------------


class MaskPolicy(object):
    KEEP   = 0
    MASK   = 1
    REMOVE = 2
    def __init__(self, flags, values, description=None):
        self._description = description
        if _my_isscalar(flags):
            if isinstance(flags, str):
                flags = flags.split(',')
            else:
                flags = (flags,)
        if _my_isscalar(values):
            if isinstance(values, str):
                values = values.split(',')
            else:
                values = (values,)
        if len(flags) != len(values):
            raise ValueError('The number of policy flags is different from the number of policies.')

        self._policy = []
        for flag, value in zip(flags, values):
            if flag[0] == '_':
                raise ValueError('A policy flag should not start with an underscore.')
            value = value.strip().lower()
            if value not in ('keep', 'mask', 'remove'):
                raise KeyError("Invalid policy "+flag+"='" + value + "'. Valid policies are 'keep', 'mask' or 'remove'.")
            self._policy.append({flag:value})
            setattr(self, flag, value)
        self._policy = tuple(self._policy)

    def __array__(self, dtype=int):
        conversion = { 'keep':self.KEEP, 'mask':self.MASK, 'remove':self.REMOVE }
        return numpy.array([conversion[policy.values()[0]] for policy in self._policy], dtype=dtype)

    def __str__(self):
        str = self._description + ': ' if self._description is not None else ''
        str_ = []
        for policy in self._policy:
            str_.append(policy.values()[0] + " '" + policy.keys()[0] + "'")
        str += ', '.join(str_)
        return str


#-------------------------------------------------------------------------------


class Pointing(numpy.ndarray):
    def __new__(cls, time, ra, dec, pa, flag, nsamples=None):
        if nsamples is None:
            nsamples = (time.size,)
        else:
            if _my_isscalar(nsamples):
                nsamples = (nsamples,)
            else:
                nsamples = tuple(nsamples)
        nsamples_tot = numpy.sum(nsamples)
        if numpy.any(numpy.array([ra.size,dec.size,pa.size,flag.size]) != nsamples_tot):
            raise ValueError('The pointing inputs do not have the same size.')

        ftype = config.get_default_dtype_float()
        dtype = numpy.dtype([('time', ftype), ('ra', ftype), ('dec', ftype), ('pa', ftype), ('flag', numpy.int64)])
        result = numpy.recarray(nsamples_tot, dtype=dtype)
        result.time = time
        result.ra   = ra
        result.dec  = dec
        result.pa   = pa
        result.flag = flag
        result.nsamples = nsamples
        return result

    @property
    def velocity(self):
        dra  = numpy.diff(self.ra)
        ddec = numpy.diff(self.dec)
        dtime = numpy.diff(self.time)
        vel = numpy.sqrt((dra*numpy.cos(self.dec[0:-1].inunit('rad')))**2 + ddec**2) / dtime
        vel.unit = 'arcsec/s'
        u = vel._unit
        vel = numpy.append(vel, vel[-1])
        # BUG: append eats the unit...
        vel._unit = u
        return vel


#-------------------------------------------------------------------------------


def plot_scan(ra, dec=None, title=None, new_figure=True):
    if hasattr(ra, 'pointing'):
        ra = ra.pointing
    if hasattr(ra, 'ra') and hasattr(ra, 'dec'):
        dec = ra.dec
        ra = ra.ra

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
    p = pyplot.plot(x, y, linewidth=2)
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


