import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pyfits
try:
    import ds9
    _imported_ds9 = True
except:
    _imported_ds9 = False

from unit import Quantity
from utils import create_fitsheader, _my_isscalar

__all__ = [ 'FitsArray', 'Map', 'Tod', 'distance', 'gaussian' ]


class FitsArray(Quantity):

    __slots__ = ('_header', '__dict__')
    def __new__(cls, data, header=None, unit=None, dtype=None, copy=True, order='C', subok=False, ndmin=0):

        if type(data) is str:
            ihdu = 0
            while True:
                try:
                    hdu = pyfits.open(data)[ihdu]
                except IndexError:
                    raise IOError('The FITS file has no data.')
                if hdu.data is not None:
                    break
                ihdu += 1
            data = hdu.data
            header = hdu.header
            copy = False
            if unit is None:
                if header.has_key('bunit'):
                    unit = header['bunit']
                elif header.has_key('QTTY____'):
                    unit = header['QTTY____'] # HCSS crap

        # get a new FitsArray instance (or a subclass if subok is True)
        result = Quantity(data, unit, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        # copy header attribute
        if header is not None:
            result.header = header
        elif hasattr(data, '_header') and data._header.__class__ is pyfits.Header:
            if copy:
                result._header = data._header.copy()
            else:
                result._header = data._header
        elif str(result.dtype) not in ('complex64', 'complex128', 'complex256'):
            result._header = create_fitsheader(result)
        else:
            result._header = None

        return result

    def __array_finalize__(self, obj):
        Quantity.__array_finalize__(self, obj)
        # obj might be None without __new__ being called... (ex: append)
        #if obj is None: return
        self._header = getattr(obj, '_header', None)

    def __array_wrap__(self, obj, context=None):
        result = Quantity.__array_wrap__(self, obj, context).view(type(self))
        result._header = self._header
        return result

    @staticmethod
    def empty(shape, header=None, unit=None, dtype=None, order=None):
        return FitsArray(numpy.empty(shape, dtype, order), header, unit, copy=False)

    @staticmethod
    def ones(shape, header=None, unit=None, dtype=None, order=None):
        return FitsArray(numpy.ones(shape, dtype, order), header, unit, copy=False)

    @staticmethod
    def zeros(shape, header=None, unit=None, dtype=None, order=None):
        return FitsArray(numpy.zeros(shape, dtype, order), header, unit, copy=False)

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header ('+str(type(header))+').')
        self._header = header

    def tofile(self, fid, sep='', format='%s'):
        super(FitsArray,self).tofile(fid, sep, format)

    def writefits(self, filename):
        """Save a FitsArray instance to a fits file given a filename
       
        If the same file already exist it overwrites it.
        """

        header = self.header.copy()
       
        unit = self.unit
        if unit != '':
            header.update('BUNIT', unit)

        if hasattr(self, 'nsamples'):
            header.update('NSAMPLES', str(self.nsamples))

        if numpy.rank(self) == 0:
            value = self.reshape((1,))
        else:
            value = self.T if numpy.isfortran(self) else self
        hdu = pyfits.PrimaryHDU(value, header)
        hdu.writeto(filename, clobber=True)

    def imsave(self, filename, colorbar=True, **kw):
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        is_interactive = matplotlib.is_interactive()
        matplotlib.interactive(False)
        dpi = 80.
        figsize = numpy.clip(numpy.max(numpy.array(self.shape[::-1])/dpi), 8, 50)
        figsize = (figsize + (2 if colorbar else 0), figsize)
        self.imshow(colorbar=colorbar, figsize=figsize, dpi=dpi, **kw)
        fig = pyplot.gcf()
        fig.savefig(filename)
        matplotlib.interactive(is_interactive)

    def imshow(self, mask=None, num=None, colorbar=True, title=None, aspect=None, interpolation='nearest', origin='lower', figsize=None, dpi=None, xlabel='', ylabel='', **kw):
        """
        A simple graphical display function for the Tod class

        mask: array-like
            True means masked.
        """

        if mask is None:
            mask = ~numpy.isfinite(self)
        else:
            mask = numpy.logical_or(mask, ~numpy.isfinite(self))

        data = numpy.ma.MaskedArray(self, mask=mask, copy=False)
        if str(data.dtype) in ('complex64', 'complex128', 'complex256'):
            data = abs(data)
        mean   = numpy.mean(data)
        stddev = numpy.std(data)
        # casting to float because of a bug numpy1.4 + matplotlib
        minval = float(max(mean - 2*stddev, numpy.min(data)))
        maxval = float(min(mean + 5*stddev, numpy.max(data)))

        fig = pyplot.figure(num=num, figsize=figsize, dpi=dpi)
        fontsize = 12. * fig.get_figheight() / 6.125

        image = pyplot.imshow(data, aspect=aspect, interpolation=interpolation, origin=origin, **kw)
        image.set_clim(minval, maxval)

        ax = pyplot.gca()
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)

        if title is not None:
            pyplot.title(title, fontsize=fontsize)
        if colorbar:
            colorbar = pyplot.colorbar()
            for tick in colorbar.ax.get_yticklabels():
                tick.set_fontsize(fontsize)

        pyplot.draw()
        return image

    _ds9_id = 1
    def ds9(self, id=None, cmap='heat', scale=('scope local', 'mode 99.5'), wait=10, zoom='to fit'):
        """
        Display the array using ds9

        Parameters
        ----------
        id : string
            A string which identifies the ds9 process. If not provided,
            the id will be 'ds9_#', where '#' is incremented each time
            a ds9 process is launched.
        cmap : string, tuple of string
            Specify the cmap access point (ex: 'heat')
        scale : string , tuple of string
            Specify the scale access point (ex: 'mode 99')
        wait : integer
            Seconds to wait for ds9 to start
        zoom : string
            Specify the zoom access point (ex: 'to fit')

        Examples
        --------
        >>> map = Map('myfits.fits')
        >>> map.ds9(id='myfits.fits', scale='histequ', cmap='invert yes')
        """
        if not _imported_ds9:
            raise RuntimeError('The library pyds9 has not been installed.')
        import os, time, xpa
        if isinstance(cmap, str):
            cmap = (cmap,)
        if isinstance(scale, str):
            scale = (scale,)

        if id is None:
            id = 'ds9_' + str(self._ds9_id)
            self._ds9_id += 1

        command = 'ds9 -title ' + id
        for option in cmap:
            command += ' -cmap ' + option
        for option in scale:
            command += ' -scale ' + option
            
        command += '&'
        os.system(command)
        for i in range(wait):
            list = xpa.xpaaccess(id, None, 1024)
            if list: break
            time.sleep(1)
	if not list:
	    raise ValueError, 'no active ds9 running for target: %s' % list
        d = ds9.ds9(id)
        d.set_np2arr(self.view(numpy.ndarray).T)
        d.set('zoom ' + zoom)
        if self.header is not None:
            d.set('wcs append', str(self.header))
        return d

#-------------------------------------------------------------------------------


class Map(FitsArray):

    __slots__ = ['coverage', 'error']

    """
    Represent a map, complemented with unit and FITS header.
    """
    def __new__(cls, data, coverage=None, error=None, header=None, unit=None, dtype=None, copy=True, order='C', subok=False, ndmin=0):

        # get a new Map instance (or a subclass if subok is True)
        result = FitsArray(data, header, unit, dtype, copy, order, True, ndmin)

        if type(data) is str:
            try:
                error = pyfits.open(data)['Error'].data
            except:
                pass
            try:
                coverage = pyfits.open(data)['Coverage'].data
            except:
                pass

        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)
        
        result.coverage = coverage
        result.error = error

        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None:
            self.coverage = None
            self.error = None
            return

        # coverage
        if hasattr(obj, 'coverage'):
            self.coverage = obj.coverage
        else:
            self.coverage = None

        # error
        if hasattr(obj, 'error'):
            self.error = obj.error
        else:
            self.error = None

    @staticmethod
    def empty(shape, coverage=None, error=None, header=None, unit=None, dtype=None, order=None):
        return Map(numpy.empty(shape, dtype, order), coverage, error, header, unit, copy=False)

    @staticmethod
    def ones(shape, coverage=None, error=None, header=None, unit=None, dtype=None, order=None):
        return Map(numpy.ones(shape, dtype, order), coverage, error, header, unit, copy=False)

    @staticmethod
    def zeros(shape, coverage=None, error=None, header=None, unit=None, dtype=None, order=None):
        return Map(numpy.zeros(shape, dtype, order), coverage, error, header, unit, copy=False)

    def copy(self, order='C'):
        return Map(self, copy=True, order=order)

    def imshow(self, mask=None, num=None, title=None, figsize=None, dpi=None, aspect='equal', **kw):
        """A simple graphical display function for the Map class"""

        if mask is None and self.coverage is not None:
            mask = self.coverage <= 0

        if self.header is not None and self.header.has_key('CRPIX1'):
            xlabel = 'Right Ascension [pixels]'
            ylabel = 'Declination [pixels]'
        else:
            xlabel = 'X'
            ylabel = 'Y'

        image = super(Map, self).imshow(mask=mask, num=num, title=title, figsize=figsize, dpi=dpi, xlabel=xlabel, ylabel=ylabel, aspect=aspect, **kw)
        return image

    def writefits(self, filename):
        super(Map, self).writefits(filename)
        if self.error is not None:
            header = create_fitsheader(self.error, extname='Error')
            pyfits.append(filename, self.error, header)
        if self.coverage is not None:
            header = create_fitsheader(self.coverage, extname='Coverage')
            pyfits.append(filename, self.coverage, header)


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    __slots__ = ('mask', 'nsamples')

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None, dtype=None, copy=True, order='C', subok=False, ndmin=0):

        # get a new Tod instance (or a subclass if subok is True)
        result = FitsArray(data, header, unit, dtype, copy, order, True, ndmin)
        
        if type(data) is str:
            try:
                mask = pyfits.open(data)['Mask'].data.view('int8')
            except:
                pass

        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        # mask attribute
        if mask is not None:
            result.mask = numpy.asarray(mask, dtype='int8')
        elif copy and hasattr(data, 'mask') and data.mask is not None and data.mask is not numpy.ma.nomask:
            result.mask = data.mask.copy()

        # nsamples attribute
        if type(data) is str and nsamples is None:
            if result.header.has_key('nsamples'):
                nsamples = result.header['nsamples'][1:-1].replace(' ', '')
                if len(nsamples) > 0:
                    nsamples = map(lambda x: int(float(x)), nsamples.split(','))
        if nsamples is None:
            return result
        shape = validate_sliced_shape(result.shape, nsamples)
        result.nsamples = shape[-1]
        if type(result.nsamples) is not tuple:
            result.nsamples = (result.nsamples,)

        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        # obj might be None without __new__ being called (ex: append, concatenate)
        if obj is None:
            self.nsamples = () if numpy.rank(self) == 0 else (self.shape[-1],)
            self.mask = None
            return
        if hasattr(obj, 'mask') and obj.mask is not numpy.ma.nomask:
            self.mask = obj.mask
        else:
            self.mask = None
        self.nsamples = getattr(obj, 'nsamples', () if numpy.rank(self) == 0 else (self.shape[-1],))

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.mask = self.mask
        result.nsamples = self.nsamples
        return result

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if not isinstance(item, Tod):
            return item
        if item.mask is not None:
            item.mask = item.mask[key]
        if not isinstance(key, tuple):
            return item
        if len(key) > 1:
            if not isinstance(key[-1], slice):
                return item
            else:
                if key[-1].start is not None or key[-1].stop is not None or key[-1].step is not None:
                    item.nsamples = (item.shape[-1],)
        return item

    def reshape(self, newdims, order='C'):
        result = super(Tod, self).reshape(newdims, order=order)
        if self.mask is not None:
            result.mask = self.mask.reshape(newdims, order=order)
        return result

    @staticmethod
    def empty(shape, mask=None, nsamples=None, header=None, unit=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(numpy.empty(shape_flat, dtype, order), mask, shape[-1], header, unit, copy=False)

    @staticmethod
    def ones(shape, mask=None, nsamples=None, header=None, unit=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(numpy.ones(shape_flat, dtype, order), mask, shape[-1], header, unit, copy=False)

    @staticmethod
    def zeros(shape, mask=None, nsamples=None, header=None, unit=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(numpy.zeros(shape_flat, dtype, order), mask, shape[-1], header, unit, copy=False)
   
    def copy(self, order='C'):
        return Tod(self, copy=True, order=order)

    def flatten(self, order='C'):
        """
        Return a copy of the array collapsed into one dimension.
        """
        result = super(self.__class__, self).flatten(order)
        result.nsamples = None
        return result

    def ravel(self, order='C'):
        """
        Return a flattened view of the array
        """
        result = super(self.__class__, self).ravel(order)
        result.nsamples = None
        return result
        
    def imshow(self, num=None, title=None, figsize=None, dpi=None, aspect='auto', **kw):
        """
        A simple graphical display function for the Map class
        """

        mask = self.mask
        xlabel = 'Sample'
        ylabel = 'Detector number'
        image = super(Tod, self).imshow(mask=mask, num=num, title=title, origin='upper', figsize=figsize, dpi=dpi, xlabel=xlabel, ylabel=ylabel, aspect=aspect, **kw)
        return image

    def __str__(self):
        if numpy.rank(self) == 0:
            if self.dtype.type in (numpy.complex64, numpy.complex128, numpy.complex256):
                return 'Tod ' + str(complex(self))
            else:
                return 'Tod ' + str(float(self))
        output = 'Tod ['
        if numpy.rank(self) > 1:
            output += str(self.shape[-2])+' detector'
            if self.shape[-2] > 1:
                output += 's'
            output += ', '
        output += str(self.shape[-1]) + ' sample'
        if self.shape[-1] > 1:
            output += 's'
        nslices = len(self.nsamples)
        if nslices > 1:
            strsamples = ','.join((str(self.nsamples[i]) for i in range(nslices)))
            output += ' in ' + str(nslices) + ' slices ('+strsamples+')'
        return output + ']'

    def writefits(self, filename):
        super(Tod, self).writefits(filename)
        if self.mask is None:
            return
        mask = numpy.abs(self.mask).view('uint8')
        header = create_fitsheader(mask, extname='Mask')
        pyfits.append(filename, mask, header)
   

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


def gaussian(shape, sigma, origin=None, resolution=1., dtype=None):
    d = distance(shape, origin=origin, resolution=resolution, dtype=dtype)
    d = numpy.exp(-d**2/(2*sigma**2)) / (2*numpy.pi*sigma**2)
    return d


#-------------------------------------------------------------------------------


def flatten_sliced_shape(shape):
    if shape is None: return shape
    if _my_isscalar(shape):
        return (int(shape),)
    return tuple(map(numpy.sum, shape))

   
#-------------------------------------------------------------------------------


def combine_sliced_shape(shape, nsamples):
    if _my_isscalar(shape):
        shape = [ shape ]
    else:
        shape = list(shape) # list makes a shallow copy
    if _my_isscalar(nsamples):
        nsamples = (int(nsamples),)
    else:
        nsamples = tuple(nsamples)
    shape.append(nsamples)
    return tuple(shape)

   
#-------------------------------------------------------------------------------

    
def validate_sliced_shape(shape, nsamples=None):
    # convert shape and nsamples to tuple
    if shape is None:
        if nsamples is None:
            return None
        shape = ()
    elif _my_isscalar(shape):
        shape = (int(shape),)
    else:
        shape = tuple(shape)
    if nsamples is not None:
        if _my_isscalar(nsamples):
            nsamples = (int(nsamples),)
        else:
            nsamples = tuple(nsamples)
    
    if len(shape) == 0:
        if nsamples is not None and len(nsamples) != 0:
            raise ValueError("The input is scalar, but nsamples is equal to '"+str(nsamples)+"'.")
        return (shape,)
    
    if nsamples is None:
        if _my_isscalar(shape[-1]):
            nsamples = int(shape[-1])
        else:
            nsamples = tuple(map(int, shape[-1]))
    else:
        if len(nsamples) == 0:
            raise ValueError("The input is not scalar and is incompatible with nsamples.")
        if numpy.any(numpy.array(nsamples) < 0):
            raise ValueError('The input nsamples has negative values.')
        elif numpy.sum(nsamples) != numpy.sum(shape[-1]):
            raise ValueError("The sliced input has an incompatible number of samples '" + str(nsamples) + "' instead of '" + str(shape[-1]) + "'.")
        if len(nsamples) == 1:
            nsamples = nsamples[0]
    
    l = list(shape[0:-1])
    l.append(nsamples)
    return tuple(l)
