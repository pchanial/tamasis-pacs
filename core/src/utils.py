import kapteyn
import numpy as np
import os
import pyfits
import re
import scipy.special
import tamasisfortran as tmf
import time
import scipy.signal

from matplotlib import pyplot
from . import var
from .datatypes import Map, create_fitsheader
from .numpyutils import _my_isscalar
from .quantity import Quantity

__all__ = [ 
    'airy_disk',
    'aperture_circular',
    'distance',
    'ds9',
    'gaussian',
    'hs',
    'mean_degrees',
    'minmax_degrees',
    'phasemask_fourquadrant',
    'plot_scan',
    'profile',
    'psd2',
]


class Ds9(object):
    """
    Helper around the ds9 package.
    
    Examples
    --------
    >>> ds9.open_in_new_window = False
    >>> d = ds9.current
    >>> d.set('scale linear')
    >>> ds9(my new map)
    """
    def __call__(self, array, origin=None, **keywords):
        if not isinstance(array, str):
            array = np.asanyarray(array)
            dtype = array.dtype
        else:
            dtype = None
        array = Map(array, dtype=dtype, copy=False, origin=origin)
        array.ds9(new=self.open_in_new_window, **keywords)

    @property
    def targets(self):
        """
        Returns a list of the ids of the running ds9 instances.
        """
        import ds9
        return ds9.ds9_targets()

    @property
    def current(self):
        """
        Return current ds9 instance.
        """
        targets = self.targets()
        if targets is None:
            return None
        return self.get(targets[-1])

    def get(self, id):
        """
        Return ds9 instance matching a specified id.
        """
        import ds9
        return ds9.ds9(id)

    open_in_new_window = True

ds9 = Ds9()


#-------------------------------------------------------------------------------


def hs(arg):
    """
    Display the attributes of an object (except methods and those starting
    with an underscore) or an ndarray with composite dtype alongside the
    names of the records. The display is truncated so that each name fits
    in one line.
    """
    import inspect
    if isinstance(arg, np.ndarray):
        names = arg.dtype.names
        if names is None:
            print(arg)
            return
        print(str(arg.size) + ' element' + ('s' if arg.size > 1 else ''))
    else:
        members = inspect.getmembers(arg, lambda x: not inspect.ismethod(x) \
                                     and not inspect.isbuiltin(x))
        members = [x for x in members if x[0][0] != '_']
        names = [x[0] for x in members]

    length = np.max(list(map(len, names)))
    lnames = np.array([names[i].ljust(length)+': ' for i in range(len(names))])
    for name, lname in zip(names, lnames):
        value = str(getattr(arg, name))[0:72-length-2]
        if len(value) == 72-length-2:
            value = value[0:-3] + '...'
        print(lname+value)

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


#-------------------------------------------------------------------------------


def plot_scan(input, map=None, title=None, new_figure=True, color='magenta',
              linewidth=2):
    if hasattr(input, 'pointing'):
        input = input.pointing
    if hasattr(input, 'ra') and hasattr(input, 'dec'):
        ra  = input.ra
        dec = input.dec
        if hasattr(input, 'removed'):
            ra = ra.copy()
            dec = dec.copy()
            ra [input.removed] = np.nan
            dec[input.removed] = np.nan
    elif type(input) in (list,tuple):
        ra  = input[0]
        dec = input[1]
    else:
        ra = None

    if not isinstance(ra, np.ndarray) or not isinstance(dec, np.ndarray):
        raise TypeError("Invalid input type '" + type(input) + "'.")

    nvalids = np.sum(np.isfinite(ra*dec))
    if nvalids == 0:
        raise ValueError('There is no valid pointing.')

    if isinstance(map, Map) and map.has_wcs():
        image = map.imshow(title=title)
        x, y = image.topixel(ra, dec)
        p = pyplot.plot(x, y, color=color, linewidth=linewidth)
        pyplot.plot(x[0], y[0], 'o', color=p[0]._color)
        return image

    crval = [mean_degrees(ra), np.nansum(dec)/nvalids]
    ra_min,  ra_max  = minmax_degrees(ra)
    dec_min, dec_max = np.nanmin(dec), np.nanmax(dec)
    cdelt = np.max((ra_max-ra_min)/1000., (dec_max-dec_min)/1000.)
    header = create_fitsheader(None, naxis=[1,1], cdelt=cdelt, crval=crval,
                               crpix=[1,1])

    proj = kapteyn.wcs.Projection(header)
    coords = np.array([ra, dec]).T
    xy = proj.topixel(coords)
    xmin = int(np.round(np.nanmin(xy[:,0])))
    xmax = int(np.round(np.nanmax(xy[:,0])))
    ymin = int(np.round(np.nanmin(xy[:,1])))
    ymax = int(np.round(np.nanmax(xy[:,1])))
    xmargin = int(np.ceil((xmax - xmin + 1) / 10.))
    ymargin = int(np.ceil((ymax - ymin + 1) / 10.))
    xmin -= xmargin
    xmax += xmargin
    ymin -= ymargin
    ymax += ymargin
    naxis = (xmax - xmin + 1, ymax - ymin + 1)
    crpix = (-xmin+2, -ymin+2)
    header = create_fitsheader(None, naxis=naxis, cdelt=cdelt, crval=crval,
                               crpix=crpix)
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


def profile(input, origin=None, bin=1., nbins=None, histogram=False):
    """
    Returns axisymmetric profile of a 2d image.
    x, y[, n] = profile(image, [origin, bin, nbins, histogram])

    Parameters
    ----------
    input: array
        2d input array
    origin: (x0,y0)
        center of the profile. (Fits convention). Default is the image center
    bin: number
        width of the profile bins (in unit of pixels)
    nbins: integer
        number of profile bins
    histogram: boolean
        if set to True, return the histogram
    """
    input = np.ascontiguousarray(input, var.FLOAT_DTYPE)
    if origin is None:
        origin = (np.array(input.shape[::-1], var.FLOAT_DTYPE) + 1) / 2
    else:
        origin = np.ascontiguousarray(origin, var.FLOAT_DTYPE)
    
    print type(input.shape[0]-origin[1])
    print type(origin[1])
    if nbins is None:
        nbins = int(max(input.shape[0]-origin[1], origin[1],
                        input.shape[1]-origin[0], origin[0]) / bin)

    x, y, n = tmf.profile_axisymmetric_2d(input.T, origin, bin, nbins)
    if histogram:
        return x, y, n
    else:
        return x, y


#-------------------------------------------------------------------------------


def psd2(input):
    input = np.asanyarray(input)
    s = np.abs(scipy.signal.fft2(input))
    for axis, n in zip((-2,-1), s.shape[-2:]):
        s = np.roll(s, n // 2, axis=axis)
    return Map(s)


#-------------------------------------------------------------------------------


def airy_disk(shape, fwhm, origin=None, resolution=1.):
    d  = distance(shape, origin=origin, resolution=resolution)
    index = np.where(d == 0)
    d[index] = 1.e-30
    d *= 1.61633 / (fwhm/2.)
    d  = (2 * scipy.special.jn(1,d)/d)**2
    d /= np.sum(d)
    return d


#-------------------------------------------------------------------------------


def aperture_circular(shape, diameter, origin=None, resolution=1.):
    array = distance(shape, origin=origin, resolution=resolution)
    m = array > diameter / 2.
    array[ m] = 0
    array[~m] = 1
    return array


#-------------------------------------------------------------------------------


def diff(array, axis=0):
    """
    Inplace discrete difference
    """
    array = np.asanyarray(array)
    rank = array.ndim
    
    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.diff(array.ravel(), rank-axis, np.asarray(array.T.shape))
    return array

    
#-------------------------------------------------------------------------------


def diffT(array, axis=0):
    """
    Inplace discrete difference transpose
    """
    array = np.asanyarray(array)
    rank = array.ndim

    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.difft(array.ravel(), rank-axis, np.asarray(array.T.shape))
    return array

    
#-------------------------------------------------------------------------------


def diffTdiff(array, axis=0):
    """
    Inplace discrete difference transpose times discrete difference
    """
    array = np.asanyarray(array)
    rank = array.ndim

    size = int(np.product(array.shape))
    
    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.difftdiff(array.ravel(), rank-axis, np.asarray(array.T.shape))
    return array

    
#-------------------------------------------------------------------------------


def distance(shape, origin=None, resolution=1.):
    """
    Returns an array whose values are the distances to a given origin.
    
    Parameters
    ----------
    shape : tuple of integer
        dimensions of the output array. For a 2d array, the first integer
        is for the Y-axis and the second one for the X-axis.
    origin : array-like in the form of (x0, y0, ...), optional
        The coordinates, according to the FITS convention (i.e. first
        column and row is one), from which the distance is calculated.
        Default value is the array center.
    resolution : array-like in the form of (dx, dy, ...), optional
        Inter-pixel distance. Default is one. If resolution is a Quantity, its
        unit will be carried over to the returned distance array

    Example
    -------
    nx, ny = 3, 3
    print distance((ny,nx))
    [[ 1.41421356  1.          1.41421356]
    [ 1.          0.          1.        ]
    [ 1.41421356  1.          1.41421356]]
    """
    if _my_isscalar(shape):
        shape = (shape,)
    else:
        shape = tuple(shape)
    shape = tuple(int(np.round(s)) for s in shape)
    rank = len(shape)

    if origin is None:
        origin = (np.array(shape[::-1], dtype=var.FLOAT_DTYPE) + 1) / 2
    else:
        origin = np.ascontiguousarray(origin, dtype=var.FLOAT_DTYPE)

    unit = getattr(resolution, '_unit', None)

    if _my_isscalar(resolution):
        resolution = np.resize(resolution, rank)
    resolution = np.asanyarray(resolution, dtype=var.FLOAT_DTYPE)

    if rank == 1:
        d = tmf.distance_1d(shape[0], origin[0], resolution[0])
    elif rank == 2:
        d = tmf.distance_2d(shape[1], shape[0], origin, resolution).T
    elif rank == 3:
        d = tmf.distance_3d(shape[2], shape[1], shape[0], origin, resolution).T
    else:
        d = _distance_slow(shape, origin, resolution, var.FLOAT_DTYPE)

    return Map(d, copy=False, unit=unit)


#-------------------------------------------------------------------------------


def _distance_slow(shape, origin, resolution, dtype):
    """
    Returns an array whose values are the distances to a given origin.

    This routine is written using np.meshgrid routine. It is slower
    than the Fortran-based `distance` routine, but can handle any number
    of dimensions.

    Refer to `tamasis.distance` for full documentation.
    """

    dim = []
    index = []
    for n, o, r in zip(reversed(shape), origin, resolution):
        index.append(slice(0,n))
    d = np.asarray(np.mgrid[index], dtype=var.FLOAT_DTYPE).T
    d -= np.asanyarray(origin) - 1.
    d *= resolution
    np.square(d, d)
    d = Map(np.sqrt(np.sum(d, axis=d.shape[-1])), dtype=dtype, copy=False)
    return d
    

#-------------------------------------------------------------------------------


def gaussian(shape, fwhm, origin=None, resolution=1., unit=None):
    if len(shape) == 2:
        sigma = fwhm / np.sqrt(8*np.log(2))
    else:
        raise NotImplementedError()
    d = distance(shape, origin=origin, resolution=resolution)
    d.unit = ''
    d = np.exp(-d**2/(2*sigma**2))
    d /= np.sum(d)
    if unit:
        d.unit = unit
    return d


#-------------------------------------------------------------------------------


def phasemask_fourquadrant(shape, phase=-1):
    array = Map.ones(shape, dtype=var.COMPLEX_DTYPE)
    array[0:shape[0]//2,shape[1]//2:] = phase
    array[shape[0]//2:,0:shape[1]//2] = phase
    return array
