import matplotlib.pyplot as pyplot
import numpy
import pyfits

from unit import Quantity
from utils import create_fitsheader, _my_isscalar

__all__ = [ 'FitsArray', 'Map', 'Tod' ]


class FitsArray(Quantity):

    __slots__ = ('_header',)
    def __new__(cls, data, header=None, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

        # get a new FitsArray instance (or a subclass if subok is True)
        result = Quantity(data, unit, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        # copy header attribute
        if header is not None:
            result.header = header
        elif copy and hasattr(data, '_header') and data._header.__class__ is pyfits.Header:
            result.header = data._header.copy()

        return result

    def __array_finalize__(self, obj):
        Quantity.__array_finalize__(self, obj)
        if obj is None: return
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
        from os.path import isfile
        from os import remove
  
        if self.header is None:
            header = create_fitsheader(reversed(self.shape))
        else:
            header = self.header
       
        unit = self.unit
        if unit != '':
            header.update('BUNIT', unit)

        hdu = pyfits.PrimaryHDU(self.T if self.flags.f_contiguous else self, header)
        try:
            remove(filename)
        except OSError:
            pass
        try:
            hdu.writeto(filename)
        except IOError:
            pass


#-------------------------------------------------------------------------------


class Map(FitsArray):
    
    """
    Represent a map, complemented with unit and FITS header.
    """
    def __new__(cls, data, header=None, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

        # get a new Map instance (or a subclass if subok is True)
        result = FitsArray(data, header, unit, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None: return

    @staticmethod
    def empty(shape, header=None, unit=None, dtype=None, order=None):
        return Map(numpy.empty(shape, dtype, order), header, unit, copy=False)

    @staticmethod
    def ones(shape, header=None, unit=None, dtype=None, order=None):
        return Map(numpy.ones(shape, dtype, order), header, unit, copy=False)

    @staticmethod
    def zeros(shape, header=None, unit=None, dtype=None, order=None):
        return Map(numpy.zeros(shape, dtype, order), header, unit, copy=False)

    def copy(self, order='C'):
        return Map(self, copy=True, order=order)

    def imshow(self, num=None, axis=True, title=None):
        """A simple graphical display function for the Map class"""

        finite = self[numpy.isfinite(self)]
        mean   = numpy.mean(finite)
        stddev = numpy.std(finite)
        minval = max(mean - 2*stddev, numpy.min(finite))
        maxval = min(mean + 5*stddev, numpy.max(finite))

        pyplot.figure(num=num)
        # HACK around bug in imshow !!!
        pyplot.imshow(self.T if self.flags.f_contiguous else self, interpolation='nearest')
        pyplot.clim(minval, maxval)
        if self.header is not None and self.header.has_key('CRPIX1'):
            pyplot.xlabel('Right Ascension')
            pyplot.ylabel("Declination")
        if title is not None:
            pyplot.title(title)
        pyplot.colorbar()
        pyplot.draw()


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    __slots__ = ('mask', 'nsamples')

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

        # get a new Tod instance (or a subclass if subok is True)
        result = FitsArray(data, header, unit, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        # mask attribute
        if mask is not None:
            result.mask = numpy.asarray(mask, dtype='int8')
        elif copy and hasattr(data, 'mask') and data.mask is not None and data.mask is not numpy.ma.nomask:
            result.mask = data.mask.copy()

        # nsamples attribute
        if nsamples is None:
            return result
        junk = validate_sliced_shape(result.shape, nsamples)
        if _my_isscalar(nsamples):
            result.nsamples = (int(nsamples),)
        else:
            result.nsamples = tuple(nsamples)

        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None: return
        if hasattr(obj, 'mask') and obj.mask is not numpy.ma.nomask:
            self.mask = obj.mask
        else:
            self.mask = None
        self.nsamples = getattr(obj, 'nsamples', (obj.shape[-1],) if numpy.rank(obj) > 0 else None)

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.mask = self.mask
        result.nsamples = self.nsamples
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

    def imshow(self, num=None, axis=True, title=None):
        """A simple graphical display function for the Tod class"""

        mean   = numpy.mean(self[numpy.isfinite(self)])
        stddev = numpy.std(self[numpy.isfinite(self)])
        minval = mean - 2*stddev
        maxval = mean + 5*stddev

        data = numpy.ma.MaskedArray(self, mask=self.mask, copy=False)
        pyplot.figure(num=num)
        pyplot.imshow(data, aspect='auto', interpolation='nearest')
        pyplot.clim(minval, maxval)
        pyplot.xlabel("Signal")
        pyplot.ylabel('Detector number')
        if title is not None:
            pyplot.title(title)
        pyplot.colorbar()
        pyplot.draw()

    def __str__(self):
        if numpy.rank(self) == 0:
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
        nsamples = int(nsamples)
    else:
        nsamples = tuple(nsamples)
        if len(nsamples) == 1:
            nsamples = nsamples[0]
    shape.append(nsamples)
    return tuple(shape)

   
#-------------------------------------------------------------------------------

    
def validate_sliced_shape(shape, nsamples=None):
    if isinstance(shape, numpy.ndarray) and len(shape.shape) == 0 or len(shape) == 0:
        if nsamples is not None and numpy.sum(nsamples) != 1:
            raise ValueError("The input is scalar, but nsamples is equal to '"+str(nsamples)+"'.")
        return ()
        
    if nsamples is not None:
        if numpy.sum(nsamples) != numpy.sum(shape[-1]):
            raise ValueError("The sliced input has an incompatible number of samples '" + str(nsamples) + "' instead of '" + str(shape[-1]) + "'.")
    else:
        nsamples = shape[-1]
    return combine_sliced_shape(shape[0:-1], nsamples)
