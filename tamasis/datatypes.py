import matplotlib.pyplot as pyplot
import numpy
import pyfits

from unit import Quantity
from utils import create_fitsheader, _my_isscalar

__all__ = [ 'FitsArray', 'Map', 'Tod' ]


class FitsArray(Quantity):

    __slots__ = ('_header',)
    def __new__(cls, data, header=None, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

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
                    unit = header('QTTY____') # HCSS crap

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

        if self.header is None:
            header = create_fitsheader(self)
        else:
            header = self.header
       
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


#-------------------------------------------------------------------------------


class Map(FitsArray):

    __slots__ = ['coverage', 'error']

    """
    Represent a map, complemented with unit and FITS header.
    """
    def __new__(cls, data, coverage=None, error=None, header=None, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

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

    def imshow(self, num=None, axis=True, title=None, coverage=None):
        """A simple graphical display function for the Map class"""

        mask = ~numpy.isfinite(self)
        if coverage is None:
            coverage = self.coverage
        if coverage is not None:
            mask = numpy.logical_or(mask, coverage <= 0)
        finite = self[~mask]
        mean   = numpy.mean(finite)
        stddev = numpy.std(finite)
        minval = max(mean - 2*stddev, numpy.min(finite))
        maxval = min(mean + 5*stddev, numpy.max(finite))

        data = numpy.ma.MaskedArray(self, mask=mask, copy=False)
        fig = pyplot.figure(num=num)
        image=pyplot.imshow(data, interpolation='nearest', origin='lower')
        image.set_clim(minval, maxval)
        if self.header is not None and self.header.has_key('CRPIX1'):
            pyplot.xlabel('Right Ascension [pixels]')
            pyplot.ylabel("Declination [pixels]")
        if title is not None:
            pyplot.title(title)
        pyplot.colorbar()
        pyplot.draw()
        return fig

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

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

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
            try:
                nsamples = result.header['nsamples'][1:-1].replace(' ', '')
                if len(nsamples) > 0:
                    nsamples = map(lambda x: int(x), nsamples.split(','))
                else:
                    nsamples = ()
            except:
                pass
        if nsamples is None:
            return result
        shape = validate_sliced_shape(result.shape, nsamples)
        result.nsamples = shape[-1]
        if _my_isscalar(nsamples):
            result.nsamples = (int(nsamples),)
        else:
            result.nsamples = tuple(nsamples)

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
        fig = pyplot.figure(num=num)
        image = pyplot.imshow(data, aspect='auto', interpolation='nearest')
        image.set_clim(minval, maxval)
        pyplot.xlabel("Sample")
        pyplot.ylabel('Detector number')
        if title is not None:
            pyplot.title(title)
        pyplot.colorbar()
        pyplot.draw()
        return fig

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

    def writefits(self, filename):
        super(Tod, self).writefits(filename)
        if self.mask is None:
            return
        mask = numpy.abs(self.mask).view('uint8')
        header = create_fitsheader(mask, extname='Mask')
        pyfits.append(filename, mask, header)
   
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
    # convert shape and nsamples to tuple
    if shape is None:
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
            nsamples = (int(shape[-1]),)
        else:
            nsamples = tuple(shape[-1])
    else:
        if len(nsamples) == 0:
            raise ValueError("The input is not scalar and is incompatible with nsamples.")
        if numpy.any(numpy.array(nsamples) < 0):
            raise ValueError('The input nsamples has negative values.')
        elif numpy.sum(nsamples) != numpy.sum(shape[-1]):
            raise ValueError("The sliced input has an incompatible number of samples '" + str(nsamples) + "' instead of '" + str(shape[-1]) + "'.")
    
    l = list(shape[0:-1])
    l.append(nsamples)
    return tuple(l)
