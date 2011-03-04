import kapteyn
import numpy as np
import os
import re
import scipy.special
import tamasisfortran as tmf
import scipy.signal

from matplotlib import pyplot
from . import var
from .numpyutils import _my_isscalar
from .wcsutils import angle_lonlat, barycenter_lonlat
from .quantity import Quantity
from .datatypes import Map, create_fitsheader

__all__ = [ 
    'airy_disk',
    'aperture_circular',
    'distance',
    'ds9',
    'gaussian',
    'hs',
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
    >>> ds9(mynewmap, 'zoom to fit')
    """
    def __call__(self, array, xpamsg=None, origin=None, **keywords):
        if not isinstance(array, str):
            array = np.asanyarray(array)
            dtype = array.dtype
        else:
            dtype = None
        array = Map(array, dtype=dtype, copy=False, origin=origin)
        array.ds9(xpamsg=xpamsg, new=self.open_in_new_window, **keywords)

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
        targets = self.targets
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
        value = str(getattr(arg, name))[0:72-length-2].replace('\n', ' ')
        if len(value) == 72-length-2:
            value = value[0:-3] + '...'
        print(lname+value)


#-------------------------------------------------------------------------------


def plot_scan(input, map=None, title=None, new_figure=True, linewidth=2, **kw):

    def plot_scan_(ra, dec, nsamples, image, linewidth=linewidth, **kw):
        x, y = image.topixel(ra, dec)
        dest = 0
        for n in nsamples:
            p = pyplot.plot(x[dest:dest+n], y[dest:dest+n],
                            linewidth=linewidth, **kw)
            for i in xrange(dest, dest+n):
                if np.isfinite(x[i]) and np.isfinite(y[i]):
                    pyplot.plot(x[i], y[i], 'o', color=p[0]._color)
                    break
            dest += n

    if type(input) not in (list, tuple):
        input = [input]
    else:
        input = list(input)

    # format the input as a list of (ra, dec, nsamples)
    for islice, slice in enumerate(input):
        if hasattr(slice, 'pointing'):
            slice = slice.pointing
        if hasattr(slice, 'ra') and hasattr(slice, 'dec'):
            ra = slice.ra
            dec = slice.dec
            if hasattr(slice, 'removed'):
                ra = ra.copy()
                dec = dec.copy()
                ra [slice.removed] = np.nan
                dec[slice.removed] = np.nan
            nsamples = getattr(slice, 'nsamples', (ra.size,))
            input[islice] = (ra, dec, nsamples)
        elif type(slice) not in (list,tuple) or len(slice) != 2:
            raise TypeError('Invalid input.')
        else:
            ra = np.array(slice[0], ndmin=1, copy=False)
            dec = np.array(slice[1], ndmin=1, copy=False)
            if ra.shape != dec.shape:
                raise ValueError('Input R.A. and Dec do not have the same len' \
                                 'gth.')
            input[islice] = (ra, dec, (ra.size,))

    npointings = sum([s[0].size for s in input])
    ra = np.empty(npointings)
    dec = np.empty(npointings)

    dest = 0
    nsamples = []
    for slice in input:
        n = slice[0].size
        ra [dest:dest+n] = slice[0]
        dec[dest:dest+n] = slice[1]
        nsamples += slice[2]
        dest += n

    nvalids = np.sum(np.isfinite(ra*dec))
    if nvalids == 0:
        raise ValueError('There is no valid coordinates.')

    if isinstance(map, Map) and map.has_wcs():
        image = map.imshow(title=title)
        plot_scan_(ra, dec, nsamples, image, linewidth=linewidth, **kw)
        return image

    crval = barycenter_lonlat(ra, dec)
    angles = angle_lonlat((ra,dec), crval)
    angle_max = np.nanmax(angles)
    if angle_max >= 90.:
        print('Warning: some coordinates have an angular distance to the proj' \
              'ection point greater than 90 degrees.')
        mask = angles >= 90
        ra[mask] = np.nan
        dec[mask] = np.nan
        angles[mask] = np.nan
        angle_max = np.nanmax(angles)
    cdelt = np.rad2deg(np.tan(np.deg2rad(angle_max))) / 100
    header = create_fitsheader(None, naxis=[201,201], cdelt=cdelt, crval=crval,
                               crpix=[101,101])

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
    plot_scan_(ra, dec, nsamples, image, linewidth=linewidth, **kw)
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

    if not isinstance(array, np.ndarray):
        raise TypeError('Input array is not an ndarray.')

    if array.dtype != var.FLOAT_DTYPE:
        raise TypeError('The data type of the input array is not ' + \
                        str(var.FLOAT_DTYPE.type) + '.')

    rank = array.ndim
    
    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.diff(array.ravel(), rank-axis, np.asarray(array.T.shape))

    
#-------------------------------------------------------------------------------


def diffT(array, axis=0):
    """
    Inplace discrete difference transpose
    """

    if not isinstance(array, np.ndarray):
        raise TypeError('Input array is not an ndarray.')

    if array.dtype != var.FLOAT_DTYPE:
        raise TypeError('The data type of the input array is not ' + \
                        str(var.FLOAT_DTYPE.type) + '.')

    rank = array.ndim

    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.difft(array.ravel(), rank-axis, np.asarray(array.T.shape))

    
#-------------------------------------------------------------------------------


def diffTdiff(array, axis=0, scalar=1.):
    """
    Inplace discrete difference transpose times discrete difference
    """

    if not isinstance(array, np.ndarray):
        raise TypeError('Input array is not an ndarray.')

    if array.dtype != var.FLOAT_DTYPE:
        raise TypeError('The data type of the input array is not ' + \
                        str(var.FLOAT_DTYPE.type) + '.')

    rank = array.ndim
    
    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.difftdiff(array.ravel(), rank-axis, np.asarray(array.T.shape),
                      scalar)

    
#-------------------------------------------------------------------------------


def shift(array, n, axis=0):
    """
    Shift array elements inplace along a given axis.

    Elements that are shifted beyond the last position are not re-introduced
    at the first.

    Parameters
    ----------
    array : float array
        Input array to be modified

    n : integer number or array
        The number of places by which elements are shifted. If it is an array,
        specific offsets are applied along the first dimensions.

    axis : int, optional
        The axis along which elements are shifted. By default, it is the first
        axis.

    Examples
    --------
    >>> a = ones(8)
    >>> shift(a, 3); a
    array([0., 0., 0., 1., 1., 1., 1., 1., 1.])

    >>> a = array([[1.,1.,1.,1.],[2.,2.,2.,2.]])
    >>> shift(a, [1,-1], axis=1); a
    array([[0., 1., 1., 1.],
           [2., 2., 2., 0.]])
    """
    if not isinstance(array, np.ndarray):
        raise TypeError('Input array is not an ndarray.')

    if array.dtype != var.FLOAT_DTYPE:
        raise TypeError('The data type of the input array is not ' + \
                        str(var.FLOAT_DTYPE.name) + '.')

    rank = array.ndim
    n = np.array(n, ndmin=1, dtype='int32').ravel()
    
    if axis < 0:
        axis = rank + axis

    if axis == 0 and n.size > 1 or n.size != 1 and n.size not in \
       np.cumproduct(array.shape[0:axis]):
        raise ValueError('The offset size is incompatible with the first dime' \
                         'nsions of the array')
    if rank == 0:
        array.shape = (1,)
        array[:] = 0
        array.shape = ()
    else:
        tmf.shift(array.ravel(), rank-axis, np.asarray(array.T.shape), n)

    
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
    print(distance((ny,nx)))
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
