try:
    import fftw3
except:
    print 'Warning: Library PyFFTW3 is not installed.'

import matplotlib.pyplot as pyplot
import numpy
import os
import pyfits
import re
import scipy.sparse.linalg   
import tamasisfortran as tmf
from   unit import Quantity, UnitError

__version_info__ = (1, 1, 0)
__version__ = '.'.join((str(i) for i in __version_info__))
__verbose__ = False
tamasis_dir = os.path.dirname(__file__)+'/' if os.path.dirname(__file__) != '' else './'



#-------------------------------------------------------------------------------
#
# Acquisition models
#
#-------------------------------------------------------------------------------


class ValidationError(Exception): pass


#-------------------------------------------------------------------------------


class AcquisitionModel(object):
    """
    Class representing the Acquisition model M as in the equation
         y = M.x,
    where x is the map and y is the instrumental response. An acquisition model can be the
    combination of several submodels describing various parts of the instrumental chain:
         M = M3 * M2 * M1 ...
    Two methods are provided to the user:
         y = M.direct(x)s
         x = M.transpose(y)
    they return M.x (submodels are applied in right-to-left order) and transpose(M).y.
    Author: P. Chanial
    """
    def __init__(self, description):
        self.description = description

    def __call__(self, data, reusein=False, reuseout=False):
        return self.direct(data, reusein=reusein, reuseout=reuseout)

    def direct(self, data, reusein=False, reuseout=False):
        import inspect
        for i, model in enumerate(self):
            reusein_  = reusein  or i != 0
            reuseout_ = reuseout or i != len(self)-1
            if 'reuseout' not in inspect.getargspec(model.direct)[0]:
                reusein_ = reusein_  and reuseout_
            data = model.direct(data, reusein=reusein_, reuseout=reuseout_)
        return data

    def transpose(self, data, reusein=False, reuseout=False):
        import inspect
        for i, model in enumerate(reversed(self)):
            reusein_  = reusein  or i != 0
            reuseout_ = reuseout or i != len(self)-1
            if 'reuseout' not in inspect.getargspec(model.transpose)[0]:
                reusein_ = reusein_  and reuseout_
            data = model.transpose(data, reusein=reusein_, reuseout=reuseout_)
        return data

    @property
    def T(self):
        return AcquisitionModelTranspose(self)

    @property
    def shapein(self):
        if len(self) == 0:
            return self._shapein
        shapein = None
        for model in reversed(self):
            shapein = model._validate_shape_transpose(shapein)
        return shapein

    @shapein.setter
    def shapein(self, shape):
        self._shapein = shape

    @property
    def shapeout(self):
        if len(self) == 0:
            return self._shapeout
        shapeout = None
        for model in self:
            shapeout = model._validate_shape_direct(shapeout)
        return shapeout

    @shapeout.setter
    def shapeout(self, shape):
        self._shapeout = shape

    def aslinearoperator(self, shape=None, unpacking=None):
        if numpy.isscalar(shape):
            shape = (shape, shape)
        shapein  = self._flatten_sliced_shape(self.shapein )
        shapeout = self._flatten_sliced_shape(self.shapeout)
        if shape is None:
            if (shapein is None or shapeout is None):
                raise ValidationError('You must specify the shape of the linear operator: the Acquisition model is not constrained.')
            shape = (shapein, shapeout)
        elif shapein is not None and shapein  != shape[1] or shapeout is not None and shapeout != shape[0]:
            raise ValidationError('The input shape is incompatible with the AcquisitionModel.')

        linopshape = (numpy.product(shape[0]), numpy.product(shape[1]))
        model = Reshaping(shape[0], linopshape[0]) * self * Reshaping(linopshape[1], shape[1])
        operator = scipy.sparse.linalg.interface.LinearOperator(linopshape, matvec=model.direct, rmatvec=model.transpose, dtype='double')
        operator.unpacking = Reshaping(linopshape[1], shape[1])
        return operator
        
    def _validate_shape(self, shapein):
        if shapein is None or self.shapein is None:
            return
        if self._flatten_sliced_shape(shapein) != self._flatten_sliced_shape(self.shapein):
            raise ValidationError('The input of '+self.__class__.__name__+' has incompatible shape '+str(shapein)+' instead of '+str(self.shapein)+'.')

    def _validate_shape_direct(self, shapein):
        if shapein is None or self.shapein is None:
            if self.shapeout is not None:
               return self.shapeout
            else:
               return shapein
        if self._flatten_sliced_shape(shapein) != self._flatten_sliced_shape(self.shapein):
            raise ValidationError('The input of '+self.__class__.__name__+' has incompatible shape '+str(shapein)+' instead of '+str(self.shapein)+'.')
        return self.shapeout

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None or self.shapeout is None:
            if self.shapein is not None:
                return self.shapein
            else:
                return shapeout
        if self._flatten_sliced_shape(shapeout) != self._flatten_sliced_shape(self.shapeout):
            raise ValidationError("The input of '"+self.__class__.__name__+".T' has incompatible shape "+str(shapeout)+' instead of '+str(self.shapeout)+'.')
        return self.shapein

    def _validate_input(self, cls, data):
        if not isinstance(data, cls):
            if issubclass(cls, data.__class__):
                data = cls(data, copy=False)
            else:
                raise TypeError("The input of '"+self.__class__.__name__+"' has an invalid type '"+data.__class__.__name__+"' instead of '"+cls.__name__+"'.")
        self._validate_shape(self._get_sliced_shape(data))
        return data

    def _validate_input_direct(self, cls, data, reusein):
        if not isinstance(data, cls):
            if issubclass(cls, data.__class__):
                data = cls(data, copy=False)
            else:
                raise TypeError("The input of '"+self.__class__.__name__+"' has an invalid type '"+data.__class__.__name__+"' instead of '"+cls.__name__+"'.")
        shapeout = self._validate_shape_direct(self._get_sliced_shape(data))
        if reusein: self._output_transpose = data
        return data, shapeout

    def _validate_input_transpose(self, cls, data, reusein):
        if not isinstance(data, cls):
            if issubclass(cls, data.__class__):
                data = cls(data, copy=False)
            else:
                raise TypeError("The input of '"+self.__class__.__name__+".T' has an invalid type '"+data.__class__.__name__+"' instead of '"+cls.__name__+"'.")
        shapein = self._validate_shape_transpose(self._get_sliced_shape(data))
        if reusein: self._output_direct = data
        return data, shapein

    # allocate memory for the output of the acquisition models for which input and output have the same shape
    # These are the models for which the direct and transpose routines only have the reusein keyword
    def _validate_output(self, array, reusein):
        if reusein:
            return array
        if __verbose__: print 'Info: Allocating '+str(numpy.product(array.shape)/2.**17)+' MiB in '+self.__class__.__name__+'.'
        return array.copy('a')

    # allocate memory for the output of the direct operation.
    # **options should not depend on the input array, since it might only be used for the first call and it would not be
    # propagated during the next calls
    def _validate_output_direct(self, cls, shapeout, reuseout, **options):
        if shapeout is None:
            raise ValueError('The shape of the output of '+self.__class__.__name__+' is not known.')
        shapeout_flat = self._flatten_sliced_shape(shapeout)
        if self._output_direct is None or shapeout_flat != self._output_direct.shape:
            if __verbose__: print 'Info: Allocating '+str(numpy.product(shapeout_flat)/2.**17)+' MiB for the output of '+self.__class__.__name__+'.'
            if cls == numpy.ndarray:
                self._output_direct = numpy.empty(shapeout, dtype=numpy.float64, **options)
            else:
                self._output_direct = cls.empty(shapeout_flat, dtype=numpy.float64, **options)
        output = self._output_direct
        if not reuseout:
            self._output_direct = None
        return output

    # allocate memory for the output of the transpose operation.
    # **options should not depend on the input array
    def _validate_output_transpose(self, cls, shapein, reuseout, **options):
        if shapein is None:
            raise ValueError('The shape of the input of '+self.__class__.__name__+' is not known.')
        shapein_flat = self._flatten_sliced_shape(shapein)
        if self._output_transpose is None or shapein_flat != self._output_transpose.shape:
            if __verbose__: print 'Info: Allocating '+str(numpy.product(shapein_flat)/2.**17)+' MiB for the output of '+self.__class__.__name__+'.T.'
            if cls == numpy.ndarray:
                self._output_transpose = numpy.empty(shapein, dtype=numpy.float64, **options)
            else:
                self._output_transpose = cls.empty(shapein_flat, dtype=numpy.float64, **options)
        output = self._output_transpose
        if not reuseout:
            self._output_transpose = None
        return output

    _shapein          = None   # input is unconstrained
    _shapeout         = None   # output is unconstrained
    _output_direct    = None   # stores the input of the transpose model. Its memory allocation is re-used as the output of the direct model
    _output_transpose = None   # stores the input of the direct model. Its memory allocation is re-used as the output of the transpose model
    blocks            = []     # components are ordered as in the direct order

    @staticmethod
    def _flatten_sliced_shape(shape):
        if shape is None: return shape
        if not isinstance(shape, tuple) and not isinstance(shape, list) and not isinstance(shape, numpy.ndarray):
            return (int(shape),)
        return tuple(map(numpy.sum, shape))

    @staticmethod
    def _combine_sliced_shape(shape, nsamples):
        if numpy.isscalar(shape):
            shape = [ shape ]
        else:
            shape = list(shape) # list makes a shallow copy
        if numpy.rank(nsamples) == 0:
            nsamples = int(nsamples)
        elif numpy.rank(nsamples) != 1:
            raise ValueError('The input nsamples is neither a scalar nor a vector.')
        else:
            if isinstance(nsamples, numpy.ndarray) and nsamples.size == 1 or len(nsamples) == 1:
                nsamples = nsamples[0]
            else:
                nsamples = tuple(nsamples)
        shape.append(nsamples)
        return tuple(shape)

    @staticmethod
    def _get_sliced_shape(data):
        if isinstance(data, Tod):
            return AcquisitionModel._combine_sliced_shape(data.shape[0:-1], data.nsamples)
        return data.shape

    @staticmethod
    def _smart_reshape(array, shape):
        curr = array
        shape = AcquisitionModel._flatten_sliced_shape(shape)
        while True:
            if curr.shape == shape:
                return curr
            if curr.base is None or curr.base.dtype != array.dtype or curr.base.__class__ != array.__class__:
                return curr.reshape(shape)
            curr = curr.base

    def __mul__(self, other):
        if numpy.isscalar(other):
            other = Scalar(other)
        elif isinstance(other, numpy.ndarray):
            return self(other)
        elif not isinstance(other, AcquisitionModel):
            raise TypeError("It is not possible to compose an AcquisitionModel with '"+str(type(other))+"'.")

        newmodel = AcquisitionModel(None)
        newmodel.blocks = []
        for block in other:
            newmodel.blocks.append(block)
        for block in self:
            newmodel.blocks.append(block)
        # shape validation
        shapein = self.shapein
        shapeout = self.shapeout
        return newmodel

    def __rmul__(self, other):
        if numpy.isscalar(other):
            return Scalar(other) * self
        raise TypeError("It is not possible to compose '"+str(type(other))+"' with an AcquisitionModel.")

    def __add__(self, other):
        return Addition(self, other)

    def __getitem__(self, index):
        if len(self) == 0:
            raise IndexError('This acquisition model has no component.')
        return self.blocks[index]

    def __setitem__(self, index, value):
        if len(self) == 0:
            raise IndexError('This acquisition model has no component.')

        try:
            oldvalue = self[index]
        except:
            raise IndexError('Only substitutions are allowed in an acquisition model.')

        self.blocks[index] = value
        # shape validation
        try:
            shapein = self.shapein
            shapeout = self.shapeout
        except ValidationError as inst:
            self.blocks[index] = oldvalue
            raise inst

    def __iter__(self):
        if len(self) == 0:
            return iter((self,))
        return self.blocks.__iter__()

    def __reversed__(self):
        if len(self) == 0:
            return iter((self,))
        return self.blocks.__reversed__()

    def __len__(self):
        return self.blocks.__len__()

    def __str__(self):
        result = self.__class__.__name__
        if self.description is not None:
            result = self.description+' ('+self.__class__.__name__+')'
        else:
            result = self.__class__.__name__
        if self.shapein is not None or self.shapeout is not None:
            result += ' [input:'
            result += 'unconstrained' if self.shapein is None else str(self.shapein).replace(' ','')
            result += ', output:'
            result += 'unconstrained' if self.shapeout is None else str(self.shapeout).replace(' ','')
            result += ']'
        if len(self) == 0:
            return result
        result += '\n  '+'\n  '.join((str(block) for block in self))
        return result


#-------------------------------------------------------------------------------


class AcquisitionModelTranspose(AcquisitionModel):

    def __init__(self, model):
        if model.description is None:
            description = model.__class__.__name__
        else:
            description = model.description
        AcquisitionModel.__init__(self, description+'.T')
        self.model = model

    def direct(self, data, reusein=False, reuseout=False):
        return self.model.transpose(data, reusein=reusein, reuseout=reuseout)

    def transpose(self, data, reusein=False, reuseout=False):
        return self.model.direct(data, reusein=reusein, reuseout=reuseout)

    @property
    def shapein(self):
        return self.model.shapeout

    @property
    def shapeout(self):
        return self.model.shapein


#-------------------------------------------------------------------------------


class Addition(AcquisitionModel):
    """
    This class adds two AcquisitionModel together.
    """
    def __init__(self, model1, model2, description=None):
        AcquisitionModel.__init__(self, description)
        self.model1 = model1
        self.model2 = model2
   
    def direct(self, data, reusein=False, reuseout=False):
        return self.model1(data, reusein=False, reuseout=reuseout) + self.model2(data, reusein=reusein, reuseout=reuseout)

    def transpose(self, data, reusein=False, reuseout=False):
        return self.model1.T(data, reusein=False, reuseout=reuseout) + self.model2.T(data, reusein=reusein, reuseout=reuseout)

    @property
    def shapein(self):
        shapein1 = self.model1.shapein
        shapein2 = self.model1.shapein
        if shapein1 is None:
            return shapein2
        if shapein2 is None:
            return shapein1
        if shapein1 != shapein2:
            raise ValidationError("The input shape of acquisitionModels '" + str(self.model1)+"' and '" + str(self.model2) + "' " \
                                  'are not compatible.')
        return shapein1

    @property
    def shapeout(self):
        shapeout1 = self.model1.shapeout
        shapeout2 = self.model1.shapeout
        if shapeout1 is None:
            return shapeout2
        if shapeout2 is None:
            return shapeout1
        if shapeout1 != shapeout2:
            raise ValidationError("The input shape of acquisitionModels '" + str(self.model1)+"' and '" + str(self.model2) + "' " \
                                  'are not compatible.')
        return shapeout1


#-------------------------------------------------------------------------------


class DiscreteDifference(AcquisitionModel):
    """Calculate the nth order discrete difference along given axis."""
    def __init__(self, n=1, axis=0, description=None):
        AcquisitionModel.__init__(self, description)
        if n != 1:
            raise NotImplemented('DiscreteDifference is not implemented for n > 1.')
        self.n = n
        self.axis = axis

    def direct(self, data, reusein=False, reuseout=False):
        data = self._validate_input(numpy.ndarray, data)
        return numpy.diff(data, self.n, self.axis)

    def transpose(self, data, reusein=False, reuseout=False):
        data = self._validate_input(numpy.ndarray, data)
        shape = list(data.shape)
        shape[self.axis] = 1
        border = numpy.zeros(shape)
        out = numpy.concatenate((border, data, border), axis=self.axis)
        return - numpy.diff(out, axis=self.axis)

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        shapeout = list(shapein)
        shapeout[self.axis] -= self.n
        return tuple(shapeout)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return None
        shapein = list(shapeout)
        shapein[self.axis] += self.n
        return tuple(shapein)


#-------------------------------------------------------------------------------


class Projection(AcquisitionModel):
    """
    This class handles the direct and transpose operations by the pointing matrix
    The input observation has the following required attributes/methods:
        - nfinesamples
        - nsamples
        - ndetectors
        - get_pointing_matrix()
        - unit
    The instance has the following specific attributes:
        - header: the FITS header of the map
        - pmatrix: opaque representation of the pointing matrix
        - npixels_per_sample: maximum number of sky map pixels that can be intercepted by a detector
    Author: P. Chanial
    """
    def __init__(self, observation, header=None, resolution=None, npixels_per_sample=None, finer_sampling=True, description=None):
        AcquisitionModel.__init__(self, description)
        self._pmatrix, self.header, ndetectors, nsamples, self.npixels_per_sample = \
            observation.get_pointing_matrix(header, resolution, npixels_per_sample, finer_sampling=finer_sampling)
        self.shapein = tuple([self.header['naxis'+str(i+1)] for i in reversed(range(self.header['naxis']))])
        self.pmatrix = self._pmatrix.view(dtype=[('weight', 'f4'), ('pixel', 'i4')])
        self.pmatrix.resize((observation.ndetectors, nsamples, self.npixels_per_sample))
        self.shapeout = self._combine_sliced_shape(observation.ndetectors, nsamples)

    def direct(self, map2d, reusein=False, reuseout=False):
        map2d, shapeout = self._validate_input_direct(Map, map2d, reusein)
        output = self._validate_output_direct(Tod, shapeout, reuseout)
        output._unit = map2d._unit
        tmf.pointing_matrix_direct(self._pmatrix, map2d.T, output.T, self.npixels_per_sample)
        return output

    def transpose(self, signal, reusein=False, reuseout=False):
        signal, shapein = self._validate_input_transpose(Tod, signal, reusein)
        output = self._validate_output_transpose(Map, shapein, reuseout)
        output.header = self.header
        output._unit = signal._unit
        tmf.pointing_matrix_transpose(self._pmatrix, signal.T, output.T,  self.npixels_per_sample)
        return output

    def get_ptp(self):
        ndetectors = self.shapeout[0]
        nsamples = numpy.sum(self.shapeout[1])
        npixels = numpy.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample, nsamples, ndetectors, npixels).T

    def __str__(self):
        return super(Projection, self).__str__()


#-------------------------------------------------------------------------------


class PacsMultiplexing(AcquisitionModel):
    """
    Performs the multiplexing of the PACS subarrays. The subarray columns are read one after the
    other, in a 0.025s cycle (40Hz).
    Author: P. Chanial
    """
    def __init__(self, pacs, description=None):
        AcquisitionModel.__init__(self, description)
        self.fine_sampling_factor = pacs.fine_sampling_factor
        self.ij = pacs.ij

    def direct(self, signal, reusein=False, reuseout=False):
        signal, shapeout = self._validate_input_direct(Tod, signal, reusein)
        output = self._validate_output_direct(Tod, shapeout, reuseout)
        output.nsamples = tuple(numpy.divide(signal.nsamples, self.fine_sampling_factor))
        tmf.pacs_multiplexing_direct(signal.T, output.T, self.fine_sampling_factor, self.ij)
        return output

    def transpose(self, signal, reusein=False, reuseout=False):
        signal, shapein = self._validate_input_transpose(Tod, signal, reusein)
        output = self._validate_output_transpose(Tod, shapein, reuseout)
        output.nsamples = tuple(numpy.multiply(signal.nsamples, self.fine_sampling_factor))
        tmf.pacs_multiplexing_transpose(signal.T, output.T, self.fine_sampling_factor, self.ij)
        return output

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        #super(PacsMultiplexing, self)._validate_shape_direct(shapein)
        if shapein[1] % self.fine_sampling_factor != 0:
            raise ValidationError('The input timeline size ('+str(shapein[1])+') is not an integer times the fine sampling factor ('+str(self.fine_sampling_factor)+').')
        shapeout = list(shapein)
        shapeout[1] = shapeout[1] / self.fine_sampling_factor
        return tuple(shapeout)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return
        super(PacsMultiplexing, self)._validate_shape_transpose(shapeout)
        shapein = list(shapeout)
        shapein[1] = shapein[1] * self.fine_sampling_factor
        return tuple(shapein)


#-------------------------------------------------------------------------------


class Compression(AcquisitionModel):
    """
    Superclass for compressing the input signal.
    Author: P. Chanial
    """
    def __init__(self, compression_direct, compression_transpose, compression_factor, description):
        AcquisitionModel.__init__(self, description)
        self._compression_direct = compression_direct
        self._compression_transpose = compression_transpose
        self.compression_factor = compression_factor

    def direct(self, signal, reusein=False, reuseout=False):
        if self.compression_factor == 1:
            if not reusein: return signal.copy('a')
            return signal
        signal, shapeout = self._validate_input_direct(Tod, signal, reusein)
        output = self._validate_output_direct(Tod, shapeout, reuseout)
        output.nsamples = tuple(numpy.divide(signal.nsamples, self.compression_factor))
        self._compression_direct(signal, output, self.compression_factor)
        return output

    def transpose(self, compressed, reusein=False, reuseout=False):
        if self.compression_factor == 1:
            if not reusein: return compressed.copy('a')
            return compressed
        compressed, shapein = self._validate_input_transpose(Tod, compressed, reusein)
        output = self._validate_output_transpose(Tod, shapein, reuseout)
        output.nsamples = tuple(numpy.multiply(compressed.nsamples, self.compression_factor))
        self._compression_transpose(compressed, output, self.compression_factor)
        return output

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        if numpy.any(numpy.array(shapein[-1]) % self.compression_factor != 0):
            raise ValidationError('The input timeline size ('+str(shapein[-1])+') is not an integer times the compression factor ('+str(self.compression_factor)+').')
        return self._combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) / self.compression_factor)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return
        return self._combine_sliced_shape(shapeout[0:-1], numpy.array(shapeout[-1]) * self.compression_factor)

    def __str__(self):
        return super(Compression, self).__str__()+' (x'+str(self.compression_factor)+')'
       

#-------------------------------------------------------------------------------


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, self.compression_average_direct, self.compression_average_transpose, compression_factor, description)

    def compression_average_direct(self, signal, output, compression_factor):
        tmf.compression_average_direct(signal.T, output.T, compression_factor)
        
    def compression_average_transpose(self, signal, output, compression_factor):
        tmf.compression_average_transpose(signal.T, output.T, compression_factor)
        

#-------------------------------------------------------------------------------


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number specified by the compression factor
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, tmf.downsampling_direct, tmf.downsampling_transpose, compression_factor, description)
      

#-------------------------------------------------------------------------------


class Identity(AcquisitionModel):
    """
    Do nothing class.
    Author: P. Chanial
    """
    def __init__(self, description=None):
        AcquisitionModel.__init__(self, description)

    def direct(self, data, reusein=False, **kw):
        data = self._validate_input(numpy.ndarray, data)
        return self._validate_output(data, reusein)

    def transpose(self, data, reusein=False, **kw):
        data = self._validate_input(numpy.ndarray, data)
        return self._validate_output(data, reusein)
       

#-------------------------------------------------------------------------------


class Scalar(AcquisitionModel):
    """
    Do nothing class.
    Author: P. Chanial
    """
    def __init__(self, value, description=None):
        AcquisitionModel.__init__(self, description)
        self.value = float(value)

    def direct(self, data, reusein=False, **kw):
        data = self._validate_input(numpy.ndarray, data)
        output = self._validate_output(data, reusein)
        output *= self.value
        return output

    def transpose(self, data, reusein=False, **kw):
        return self.direct(data, reusein=reusein)
       

#-------------------------------------------------------------------------------


class Masking(AcquisitionModel):
    """
    Apply a mask: where the mask is non-null or True, the input values are set to zero
    Author: P. Chanial
    """
    def __init__(self, mask, description=None):
        AcquisitionModel.__init__(self, description)
        if mask is None:
            self.mask = None
            self.shapein = None
            self.shapeout = None
        else:
            self.mask = mask
            self.shapein = mask.shape
            self.shapeout = mask.shape

    def direct(self, data, reusein=False, **kw):
        data = self._validate_input(numpy.ndarray, data)
        output = self._validate_output(data, reusein)
        status = tmf.masking(output.T, self.mask.T)
        if status != 0: raise RuntimeError()
        return output

    def transpose(self, data, reusein=False, **kw):
        return self.direct(data, reusein=reusein)
  
    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = None
        elif not isinstance(mask, numpy.ndarray):
            raise TypeError('Incorrect type for the input mask ('+str(type(mask))+').')
        elif numpy.rank(mask) == 0:
            raise ValueError('The input mask should not be scalar.')
        elif self.shapein is not None and self.shapein != mask.shapein:
            raise ValueError('The input mask has an incompatible shape '+str(mask.shape)+' instead of '+str(self.shapein)+'.')
        else:
            self._mask = numpy.array(mask, dtype='int8', copy=False)


#-------------------------------------------------------------------------------


class Unpacking(AcquisitionModel):
    """
    Convert 1d map into 2d map, under the control of a mask (true means non-observed)
    Author: P. Chanial
    """
    def __init__(self, mask, field=0., description=None):
        AcquisitionModel.__init__(self, description)
        self.mask = numpy.array(mask, dtype='int8', copy=False)
        self.field = field
        self.shapein  = (numpy.sum(self.mask == 0),)
        self.shapeout = tuple(self.mask.shape)

    def direct(self, packed, reusein=False, reuseout=False):
        packed, shapeout = self._validate_input_direct(numpy.ndarray, packed, reusein)
        output = self._validate_output_direct(Map, shapeout, reuseout)
        tmf.unpack_direct(packed.T, self.mask.T, output.T, self.field)
        return output

    def transpose(self, unpacked, reusein=False, reuseout=False):
        unpacked, shapein = self._validate_input_transpose(Map, unpacked, reusein)
        output = self._validate_output_transpose(numpy.ndarray, shapein, reuseout)
        tmf.unpack_transpose(unpacked.T, self.mask.T, output.T)
        return output


#-------------------------------------------------------------------------------


class Reshaping(AcquisitionModel):
    """
    Reshape arrays
    Author: P. Chanial
    """
    def __init__(self, shapein, shapeout, description=None):
        AcquisitionModel.__init__(self, description)
        self.shapein  = self._flatten_sliced_shape(shapein)
        self.shapeout = self._flatten_sliced_shape(shapeout)
        if numpy.product(self.shapein) != numpy.product(self.shapeout):
            raise ValueError('The number of elements of the input and output of the Reshaping operator are incompatible.')

    def direct(self, array, reusein=False, **kw):
        self._validate_input_direct(numpy.ndarray, array, False)
        output = self._validate_output(array, reusein)
        output = self._smart_reshape(output, self.shapeout)
        return output

    def transpose(self, array, reusein=False, **kw):
        self._validate_input_transpose(numpy.ndarray, array, False)
        output = self._validate_output(array, reusein)
        return self._smart_reshape(output, self.shapein)


#-------------------------------------------------------------------------------


class Padding(AcquisitionModel):
    "Pads before and after a Tod."
    def __init__(self, left=0, right=0, value=0., description=None):
        AcquisitionModel.__init__(self, description)
        left = numpy.array(left, ndmin=1)
        right = numpy.array(right, ndmin=1)
        if numpy.any(left < 0):
            raise ValueError('Left padding is not positive.')
        if numpy.any(right < 0):
            raise ValueError('Right padding is not positive.')
        if numpy.rank(left) != 1 or numpy.rank(right) != 1:
            raise ValueError('Padding must be scalar or a vector.')
        self.left  = tuple(left)
        self.right = tuple(right)
        self.value = value
   
    def direct(self, array, reusein=False, reuseout=False):
        array, shapeout = self._validate_input_direct(Tod, array, reusein)
        output = self._validate_output_direct(Tod, shapeout, reuseout)
        output.nsamples = tuple(numpy.array(array.nsamples) + self.left + self.right)
        dest = 0
        dest_padded = 0
        for islice in range(len(array.nsamples)):
            nsamples = array.nsamples[islice]
            left = self.left[islice if len(self.left) > 1 else 0]
            output[...,dest_padded:dest_padded+left] = self.value
            output[...,dest_padded+left:dest_padded+left+nsamples] = array[...,dest:dest+nsamples]
            output[...,dest_padded+left+nsamples:] = self.value
            dest += nsamples
            dest_padded += output.nsamples[islice]
        return output
   
    def transpose(self, array, reusein=False, reuseout=False):
        array, shapeout = self._validate_input_transpose(Tod, array, reusein)
        output = self._validate_output_transpose(Tod, shapeout, reuseout)
        output.nsamples = tuple(numpy.array(array.nsamples) - self.left - self.right)
        dest = 0
        dest_padded = 0
        for islice in range(len(array.nsamples)):
            nsamples = output.nsamples[islice]
            left  = self.left [islice if len(self.left)  > 1 else 0]
            output[...,dest:dest+nsamples] = array[...,dest_padded+left:dest_padded+left+nsamples]
            dest += nsamples
            dest_padded += array.nsamples[islice]
        return output

    def _validate_input_direct(self, cls, data, reusein):
        data, shapeout = super(Padding, self)._validate_input_direct(cls, data, reusein)
        if len(self.left) != 1 and len(self.left) != len(data.nsamples):
            raise ValueError("The input Tod has a number of slices '" + str(len(data.nsamples)) +
                             "' incompatible with the specified padding.")
        return data, shapeout
       
    def _validate_input_transpose(self, cls, data, reusein):
        data, shapein = super(Padding, self)._validate_input_transpose(cls, data, reusein)
        if len(self.left) != 1 and len(self.left) != len(data.nsamples):
            raise ValueError("The input Tod has a number of slices '" + str(len(data.nsamples)) +
                             "' incompatible with the specified padding.")
        return data, shapein       

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        return self._combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) + self.left + self.right)
       
    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return None
        return self._combine_sliced_shape(shapeout[0:-1], numpy.array(shapeout[-1]) - self.left - self.right)
       

#-------------------------------------------------------------------------------


class Fft(AcquisitionModel):
    """
    Performs real fft
    """
    def __init__(self, nsamples, description=None):
        AcquisitionModel.__init__(self, description)
        self.nsamples = tuple(numpy.array(nsamples, ndmin=1))
        self.tarray = []
        self.farray = []
        self.forward_plan = []
        self.backward_plan = []
        for n in self.nsamples:
            tarray = numpy.empty(n, dtype='double')
            farray = numpy.empty(n, dtype='double')
            self.tarray.append(tarray)
            self.farray.append(farray)
            self.forward_plan.append(fftw3.Plan(tarray, farray, direction='forward', flags=['measure'], realtypes=['halfcomplex r2c']))
            self.backward_plan.append(fftw3.Plan(farray, tarray, direction='backward', flags=['measure'], realtypes=['halfcomplex c2r']))

    def direct(self, array, reusein=False, **kw):
        array = self._validate_input(numpy.ndarray, array)
        output = self._validate_output(array, reusein)
        output = self._smart_reshape(output, (numpy.product(array.shape[:-1]), array.shape[-1]))
        dest = 0
        for iobs, n in enumerate(self.nsamples):
            for i in range(output.shape[0]):
                self.tarray[iobs][:] = output[i,dest:dest+n]
                fftw3.execute(self.forward_plan[iobs])
                output[i,dest:dest+n] = self.farray[iobs]
            dest += n
        output = self._smart_reshape(output, array.shape)
        return output

    def transpose(self, array, reusein=False, **kw):
        array = self._validate_input(numpy.ndarray, array)
        output = self._validate_output(array, reusein)
        output = self._smart_reshape(output, (numpy.product(array.shape[:-1]), array.shape[-1]))
        dest = 0
        for iobs, n in enumerate(self.nsamples):
            for i in range(output.shape[0]):
                self.farray[iobs][:] = output[i,dest:dest+n]
                fftw3.execute(self.backward_plan[iobs])
                output[i,dest:dest+n] = self.tarray[iobs]
            output[:,dest:dest+n] /= n
            dest += n
        output = self._smart_reshape(output, array.shape)
        return output

    def _validate_shape(self, shapein):
        if shapein is None:
            return None
        nsamples = shapein[-1]
        if not isinstance(nsamples, tuple) and not isinstance(nsamples, list) and not isinstance(nsamples, numpy.ndarray):
            nsamples = (int(nsamples), )
        else:
            nsamples = tuple(nsamples)
        if nsamples != self.nsamples:
            raise ValidationError("Invalid FFT size '"+str(nsamples)+"' instead of '"+str(self.nsamples)+"'.")
        return shapein


#-------------------------------------------------------------------------------


class InvNtt(AcquisitionModel):

    def __init__(self, shape, filename, convert='native', description=None):
        AcquisitionModel.__init__(self, description)
        ndetectors = shape[0]
        nsamples = numpy.array(shape[1], ndmin=1, dtype='int32')
        nsamples_tot = numpy.sum(nsamples)
        tod_filter, ncorrelations, status = tmf.invntt_madmap1(filename, convert, nsamples, nsamples_tot, ndetectors)
        if status != 0: raise RuntimeError()
        self.data = Tod(tod_filter.T, copy=False)
        self.ncorrelations = ncorrelations
        self.shapein  = shape
        self.shapeout = shape
       
    def direct(self, data, reusein=False, **kw):
        return data

    def transpose(self, data, reusein=False, **kw):
        data = self._validate_input(Tod, data)
        output = self._validate_output(data, reusein)
        output *= self.data
        return output



#-------------------------------------------------------------------------------


class SqrtInvNtt(AcquisitionModel):

    def __init__(self, observation, filename, convert='native', description=None):
        AcquisitionModel.__init__(self, description)
        nslices = len(observation.nsamples)
        nsamples = numpy.array(observation.nsamples)
        nsamples_tot = numpy.sum(nsamples)
        ndetectors = observation.ndetectors
        tod_filter, ncorrelations, status = tmf.sqrt_invntt_madmap1(filename, convert, nsamples, nsamples_tot, ndetectors)
        if status != 0: raise RuntimeError()
        self.data = Tod(tod_filter.T, copy=False)
        self.ncorrelations = ncorrelations
        self.shapein  = (ndetectors, nsamples_tot)
        self.shapeout = self.shapein
       
    def direct(self, data, reusein=False, **kw):
        data = self._validate_input(numpy.ndarray, data)
        output = self._validate_output(data, reusein)
        output *= self.data
        return output

    def transpose(self, data, reusein=False, **kw):
        return self.direct(data, reusein=reusein)



#-------------------------------------------------------------------------------
#
# PACS observations & simulations
#
#-------------------------------------------------------------------------------



class _Pacs():

    """
    Base class which encapsulates handy information about the PACS instrument and the processed
    observation. See subclasses PacsObservation and PacsSimulation for details.
    - observing_mode      : 'prime', 'parallel' or 'transparent'
    - fine_sampling_factor : number of frames to be processed during the 0.025s period.
    - compression_factor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    Author: P. Chanial
    """

    sampling = 0.025 # 40 Hz sampling xXX not exact value! rather 0.024996

    def __init__(self, array, pointing_time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_detector_mask, keep_bad_detectors):

        from os import getenv
        from os.path import join

        if array.lower() not in ('blue', 'green', 'red'):
            raise ValueError("The input array is not 'blue', 'green' nor 'red'.")
        self.shapein  = tuple(numpy.array(self._flatten_sliced_shape(shapein), ndmin=1))
        self.shapeout = tuple(numpy.array(self._flatten_sliced_shape(shapeout), ndmin=1))

        if observing_mode is not None:
            observing_mode = observing_mode.lower()

        if observing_mode not in (None, 'prime', 'parallel', 'transparent'):
            raise ValueError("Observing mode is not 'prime', 'parallel' nor 'transparent'.")

        if compression_factor is None:
            if observing_mode is None:
                raise ValueError('The compression factor is not specified.')
            compression_factor = {'prime':4, 'parallel':8, 'transparent':1}[observing_mode]

        self.transparent_mode = observing_mode == 'transparent'

        if (fine_sampling_factor not in (1,2,4,8,16)):
            raise ValueError('Invalid fine sampling factor. It may be 1, 2, 4, 8 or 16')

        # 'blue', 'green' or 'red' array
        self.array = array.lower()

        # astrometry
        if pointing_time is None or ra is None or dec is None or pa is None or chop is None:
            raise ValueError('The simulated scan astrometry is not defined.')
        shape = pointing_time.shape
        if ra.shape != shape or dec.shape != shape or pa.shape != shape or chop.shape != shape:
            raise ValueError('The input time, ra, dec, pa, chope do not have the same shape.')
        self.pointing_time = pointing_time
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.chop = chop
        self.header = None

        # sampling information
        self.fine_time = fine_time
        self.nfinesamples = fine_time.size
        self.npixels_per_sample = npixels_per_sample
        self.fine_sampling_factor = fine_sampling_factor
        self.compression_factor = compression_factor
       
        if bad_detector_mask is None:
            bad_detector_maskFile = join(getenv('TAMASIS_DIR'),'data','PCalPhotometer_BadPixelMask_FM_v3.fits')
            bad_detector_mask = numpy.array(pyfits.fitsopen(bad_detector_maskFile)[self.array].data, dtype='int8')
        else:
            array_shapes = {'blue':(32,64), 'green':(32,64), 'red':(16,32)}
            if bad_detector_mask.shape != array_shapes[self.array]:
                raise ValueError('Input bad pixel mask has incorrect shape '+str(bad_detector_mask.shape)+' instead of '+str(array_shapes[self.array]))
            bad_detector_mask = numpy.array(bad_detector_mask, dtype='int8', copy=False)
        self.bad_detector_mask = bad_detector_mask
        if self.transparent_mode:
            self.bad_detector_mask[:,0:16] = True
            self.bad_detector_mask[16:32,16:32] = True
            self.bad_detector_mask[:,32:] = True
        self.keep_bad_detectors = keep_bad_detectors

        self.ndetectors = self.bad_detector_mask.size
        if not self.keep_bad_detectors:
            self.ndetectors -= int(numpy.sum(self.bad_detector_mask))
       
        #XXX NOT IMPLEMENTED!
        self.ij = None


#-------------------------------------------------------------------------------


class PacsObservation(_Pacs):
    """
    Class which encapsulates handy information about the PACS instrument and the processed
    observation.
    It is assumed that the detectors have sharp edges, so that the integration is simply equal
    to the sum of the sky pixels values weighted by their intersections with the detector surface.
    It contains the following attributes:
    - filename           : name of the file name, including the array colour, but excluding
                           the type. Example: '1342184520_blue'.
    - npixels_per_sample : number of sky pixels which intersect a PACS detector
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - bad_detector_mask  : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, filename, fine_sampling_factor=1, bad_detector_mask=None, keep_bad_detectors=False, mask_bad_line=False):

        filename_, nfilenames = self._files2tmf(filename)

        channel, status = tmf.pacs_info_channel(filename_, nfilenames)
        if status != 0: raise RuntimeError()

        nrows, ncolumns = (16,32) if channel == 'r' else (32,64)

        if bad_detector_mask is not None:
            if bad_detector_mask.shape != (nrows, ncolumns):
                raise ValueError('Invalid shape of the input '+('red' if channel == 'r' else 'blue')+' bad detector mask: '+str(bad_detector_mask.shape)+'.')
            bad_detector_mask = numpy.array(bad_detector_mask, dtype='int8', copy=False)

        else:
            bad_detector_mask = numpy.ones((nrows,ncolumns), dtype='int8')

        # retrieve information from the observations
        ndetectors, bad_detector_mask, transparent_mode, compression_factor, nsamples, unit, detector_area, dflat, oflat, status = tmf.pacs_info(tamasis_dir, filename_, nfilenames, fine_sampling_factor, keep_bad_detectors, numpy.asfortranarray(bad_detector_mask), mask_bad_line)
        if status != 0: raise RuntimeError()

        self.filename = filename
        self.channel = channel
        self.nrows = nrows
        self.ncolumns = ncolumns
        self.default_npixels_per_sample = 11 if channel == 'r' else 6
        self.default_resolution = 6. if channel == 'r' else 3.
        self.nobservations = nfilenames
        self.nsamples = tuple(nsamples)
        self.nfinesamples = tuple(nsamples * compression_factor * fine_sampling_factor)
        self.ndetectors = ndetectors
        self.bad_detector_mask = numpy.ascontiguousarray(bad_detector_mask)
        self.keep_bad_detectors = keep_bad_detectors
        self.mask_bad_line = mask_bad_line
        self.fine_sampling_factor = fine_sampling_factor
        self.transparent_mode = transparent_mode
        self.compression_factor = compression_factor
        self.unit = unit.strip()+' / detector'
        self.detector_area = Map(detector_area, unit='arcsec^2/detector')
        self.flatfield = {
            'total'   : Map(dflat*oflat),
            'detector': Map(numpy.ascontiguousarray(dflat)),
            'optical' : Map(numpy.ascontiguousarray(oflat))
            }
        
        #_Pacs.__init__(self, channel, npixels, nsamples, ndetectors_per_sample, fine_sampling_factor, compression_factor, bad_detector_mask, keep_bad_detectors)

    def get_map_header(self, resolution=None, finer_sampling=True):
        if resolution is None:
            resolution = self.default_resolution
        filename_, nfilenames = self._files2tmf(self.filename)
        header, status = tmf.pacs_map_header(tamasis_dir, filename_, nfilenames, finer_sampling, self.fine_sampling_factor, self.keep_bad_detectors, numpy.asfortranarray(self.bad_detector_mask), self.mask_bad_line, resolution)
        if status != 0: raise RuntimeError()
        header = str2fitsheader(header)
        return header
   
    def get_tod(self, unit=None, do_flatfielding=True, do_subtraction_mean=True):
        """
        Returns the signal and mask timelines.
        """
        filename_, nfilenames = self._files2tmf(self.filename)
        signal, mask, status = tmf.pacs_timeline(tamasis_dir, filename_, self.nobservations, numpy.sum(self.nsamples), self.ndetectors, self.keep_bad_detectors, numpy.asfortranarray(self.bad_detector_mask), self.mask_bad_line, do_flatfielding, do_subtraction_mean)
        if status != 0: raise RuntimeError()
       
        tod = Tod(signal.T, mask.T, nsamples=self.nsamples, unit=self.unit)

        # the flux calibration has been done by using HCSS photproject and assuming that the central detectors had a size of 3.2x3.2
        # squared arcseconds. To be consistent with detector sharp edge model, we need to adjust the Tod.
        if tod.unit == 'Jy / detector':
            tod *= numpy.mean(self.detector_area[15:18,31:33])/3.2**2

        if unit is None:
            return tod

        newunit = Quantity(1., unit)
        newunit_si = newunit.unit_si._unit
        if 'sr' in newunit_si and newunit_si['sr'] == -1:
            area = self.detector_area[self.bad_detector_mask != True].reshape((self.ndetectors,1))
            tod /= area
           
        tod.unit = newunit._unit
        return tod

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, finer_sampling=True):
        nsamples = numpy.sum(self.nfinesamples if finer_sampling else self.nsamples)
        if npixels_per_sample is None:
            npixels_per_sample = self.default_npixels_per_sample
        if header is None:
            header = self.get_map_header(resolution, finer_sampling)
        elif isinstance(header, str):
            header = str2fitsheader(header)

        filename_, nfilenames = self._files2tmf(self.filename)
        sizeofpmatrix = npixels_per_sample * nsamples * self.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
           
        status = tmf.pacs_pointing_matrix_filename(tamasis_dir, filename_, self.nobservations, finer_sampling, self.fine_sampling_factor, npixels_per_sample, nsamples, self.ndetectors, self.keep_bad_detectors, numpy.asfortranarray(self.bad_detector_mask), self.mask_bad_line, str(header).replace('\n', ''), pmatrix)
        if status != 0: raise RuntimeError()

        return pmatrix, header, self.ndetectors, nsamples, npixels_per_sample

    @staticmethod
    def _files2tmf(filename):
        if isinstance(filename, str):
            return filename, 1

        nfilenames = len(filename)
        length = max(len(f) for f in filename)
        filename_ = ''
        for f in filename:
            filename_ += f + (length-len(f))*' '
        return filename_, nfilenames


#-------------------------------------------------------------------------------


class PacsSimulation(_Pacs):
    """
    Class which encapsulates handy information describing the simulation and the setting of the PACS instrument.
    It contains the following attributes:
    - time               : time vector. time[0] is the time of the first sky projection
    - ra                 : boresight right ascension vector
    - dec                : boresight declination vector
    - pa                 : position Angle vector
    - chop               : chop angle vector
    - header             : pyfits header of the input sky map
    - ndetectors         : number of detectors involved in the processing
    - npixels_per_sample   : number of sky pixels which intersect a PACS detector
    - observing_mode      : 'prime', 'parallel' or 'transparent'
    - fine_sampling_factor : number of frames to be processed during the 1/40Hz period.
    - compression_factor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - bad_detector_mask  : (nx,ny) mask of int8 values (0 or 1). 1 means dead pixel.
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, inputmap, time=None, ra=None, pa=None, chop=None, array='blue', npixels_per_sample=9, observing_mode='prime', fine_sampling_factor=1, compression_factor=None, bad_detector_mask=None, keep_bad_detectors=False):

        if not isinstance(inputmap, Map):
            raise TypeError('The input is not a Map.')

        if inputmap.header is None:
            raise TypeError('The input map header is not known.')

        fine_time = numpy.arange(time[0], floor((time[-1] - time[0]) / pacs.sampling) * fine_sampling_factor, _Pacs.sampling/fine_sampling_factor)

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_detector_mask, keep_bad_detectors)

        self.inputmap
        self.header = inputmap.header

    def get_tod(self, model, noise=None):
        """
        Returns simulated signal and mask (=None) timelines.
        """
        signal = model.direct(self.inputmap)
        if self.keep_bad_detectors:
            mask = numpy.zeros(signal.shape, dtype=bool)
            for idetector, badpixel in enumerate(self.bad_detector_mask.flat):
                if badpixel != 0:
                    mask[idetector,:] = True
        else:
            mask = None

        return Tod(signal, mask=mask, nsamples=self.nsamples)

   

#-------------------------------------------------------------------------------
#
# MADmap1 observations
#
#-------------------------------------------------------------------------------


class MadMap1Observation(object):
    """Class for the handling of an observation in the MADMAP1 format"""
    def __init__(self, todfile, invnttfile, mapmaskfile, convert, ndetectors, missing_value=None):
        import re
        nslices, status = tmf.read_madmap1_nslices(invnttfile, ndetectors)
        if (status != 0): raise RuntimeError()
        self.npixels_per_sample, nsamples, status = tmf.read_madmap1_info(todfile, invnttfile, convert, ndetectors, nslices)
        if (status != 0): raise RuntimeError()

        self.todfile = todfile
        self.invnttfile = invnttfile
        m=re.search(r'(?P<filename>.*)\[(?P<extname>\w+)\]$', mapmaskfile)
        if m is None:
            mask = pyfits.fitsopen(mapmaskfile)[0].data
        else:
            filename = m.group('filename')
            extname  = m.group('extname')
            mask = pyfits.fitsopen(filename)[extname].data
        if mask is None:
            raise IOError('HDU '+mapmaskfile+' has no data.')
        self.mapmask = numpy.zeros(mask.shape, dtype='int8')
        if missing_value is None:
            self.mapmask[mask != 0] = 1
        elif numpy.isnan(missing_value):
            self.mapmask[numpy.isnan(mask)] = 1
        elif numpy.isinf(missing_value):
            self.mapmask[numpy.isinf(mask)] = 1
        else:
            self.mapmask[mask == missing_value] = 1
        self.convert = convert
        self.ndetectors = ndetectors
        self.nsamples = tuple(nsamples)
        self.nfinesamples = tuple(nsamples)

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, finer_sampling=False):
        if npixels_per_sample is not None and npixels_per_sample != self.npixels_per_sample:
            raise ValueError('The npixels_per_sample value is incompatible with the MADMAP1 file.')
        if header is not None:
            raise ValueError('The map header cannot be specified for MADmap1 observations.')
        if resolution is not None:
            raise ValueError('The map resolution cannot be specified for MADmap1 observations.')
        header = pyfits.Header()
        header.update('simple', True)
        header.update('bitpix', -64)
        header.update('extend', True)
        header.update('naxis', 1)
        header.update('naxis1', numpy.sum(self.mapmask == 0))

        tod = Tod.zeros((self.ndetectors,self.nsamples))
        sizeofpmatrix = self.npixels_per_sample * numpy.sum(self.nsamples) * self.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.read_madmap1(self.todfile, self.invnttfile, self.convert, self.npixels_per_sample, tod.T, pmatrix)
        if (status != 0): raise RuntimeError()
        return pmatrix, header, self.ndetectors, numpy.sum(self.nsamples), self.npixels_per_sample

    def get_tod(self):
        tod = Tod.zeros((self.ndetectors,self.nsamples))
        sizeofpmatrix = self.npixels_per_sample * numpy.sum(self.nsamples) * self.ndetectors
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.read_madmap1(self.todfile, self.invnttfile, self.convert, self.npixels_per_sample, tod.T, pmatrix)
        if (status != 0): raise RuntimeError()
        return tod



#-------------------------------------------------------------------------------
#
# Map & Tod class
#
#-------------------------------------------------------------------------------



class FitsArray(Quantity):

    def __new__(cls, data, header=None, unit=None, dtype=None, order='C', copy=True):
        from copy import copy as cp
        result = Quantity(data, unit, dtype=dtype, order=order, copy=copy)
        if not isinstance(result, FitsArray):
            result = result.view(cls)
        if header is not None:
            result.header = header
        elif isinstance(data, FitsArray):
            result.header = cp(data.header) if copy else data.header
        else:
            result.header = None
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
            header = FitsArray.create_header(self.shape)
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

    @staticmethod
    def create_header(shape):
        if numpy.isscalar(shape):
            shape = (shape,)
        else:
            shape = tuple(shape)
        naxis = len(shape)
        header = pyfits.Header()
        header.update('simple', True)
        header.update('bitpix', -64)
        header.update('extend', True)
        header.update('naxis', naxis)
        for dim in range(naxis):
            header.update('naxis'+str(dim+1), shape[naxis-dim-1])
        return header

#-------------------------------------------------------------------------------


class Map(FitsArray):

    def __new__(cls, data, header=None, unit=None, dtype=None, order='C', copy=True):
        result = FitsArray(data, header, unit, dtype, order, copy)
        if not isinstance(data,Map): result = result.view(cls)
        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None: return

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.header = self.header
        return result

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
        return Map(self, order=order, copy=True)

    @staticmethod
    def create_header(naxis, crval=(0.,0.), crpix=None, ctype=('RA---TAN','DEC--TAN'), cunit='deg', cdelt=None, cd=None):

        if numpy.isscalar(naxis):
            naxis = (naxis, naxis)
        else:
            naxis = tuple(naxis)

        if len(naxis) != 2:
            raise ValueError("Invalid dimension '"+str(len(naxis))+"' instead of 2.")

        header = super(Map, self).create_header(naxis)

        if cd is not None:
            if cd.shape != (2,2):
                raise ValueError('The CD matrix is not a 2x2 matrix.')
        else:
            if cdelt is None:
                return header
            if numpy.isscalar(cdelt):
                cdelt = (-cdelt, cdelt)
            cd = numpy.array(((cdelt[0], 0), (0,cdelt[1])))

        if len(crval) != 2:
            raise ValueError('CRVAL does not have two elements.')
        if crpix is None:
            crpix = numpy.array(naxis) / 2 + 1
        if len(crpix) != 2:
            raise ValueError('CRPIX does not have two elements.')
        if len(ctype) != 2:
            raise ValueError('CTYPE does not have two elements.')
        if numpy.isscalar(cunit):
            cunit = (cunit, cunit)

        header.update('crval1', crval[0])
        header.update('crval2', crval[1])
        header.update('crpix1', crpix[0])
        header.update('crpix2', crpix[1])
        header.update('cd1_1', cd[0,0])
        header.update('cd2_1', cd[1,0])
        header.update('cd1_2', cd[0,1])
        header.update('cd2_2', cd[1,1])
        header.update('ctype1', ctype[0])
        header.update('ctype2', ctype[1])
        header.update('cunit1', cunit[0])
        header.update('cunit2', cunit[1])
        return header

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

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None, dtype=None, order='C', copy=False):
        result = FitsArray(data, header, unit, dtype, order, copy)
        if not isinstance(data, Tod):
            result = result.view(cls)
        if mask is None and isinstance(data, Tod):
            mask = data.mask
        if mask is not None:
            mask = numpy.array(mask, dtype='int8', copy=copy, order=order)
        result.mask = mask
        if nsamples is None:
            if isinstance(data, Tod):
                result.nsamples = data.nsamples
            else:
                result.nsamples = (data.shape[-1],)
        else:
            junk, result.nsamples = Tod._validate_shape(data.shape, nsamples)
       
        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None: return
        self.mask = getattr(obj, 'mask', None)
        self.nsamples = getattr(obj, 'nsamples', None)

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.mask = self.mask
        result.nsamples = self.nsamples
        return result

    @staticmethod
    def empty(shape, mask=None, nsamples=None, header=None, unit=None, dtype=None, order=None):
        shape, nsamples = Tod._validate_shape(shape, nsamples)
        return Tod(numpy.empty(shape, dtype, order), mask, nsamples, header, unit, copy=False)

    @staticmethod
    def ones(shape, mask=None, nsamples=None, header=None, unit=None, dtype=None, order=None):
        shape, nsamples = Tod._validate_shape(shape, nsamples)
        return Tod(numpy.ones(shape, dtype, order), mask, nsamples, header, unit, copy=False)

    @staticmethod
    def zeros(shape, mask=None, nsamples=None, header=None, unit=None, dtype=None, order=None):
        shape, nsamples = Tod._validate_shape(shape, nsamples)
        return Tod(numpy.zeros(shape, dtype, order), mask, nsamples, header, unit, copy=False)
   
    def copy(self, order='C'):
        return Tod(self, order=order, copy=True)

    def imshow(self, num=None, axis=True, title=None):
        """A simple graphical display function for the Tod class"""

        mean   = numpy.mean(self[numpy.isfinite(self)])
        stddev = numpy.std(self[numpy.isfinite(self)])
        minval = mean - 2*stddev
        maxval = mean + 5*stddev

        pyplot.figure(num=num)
        pyplot.imshow(self.T if self.flags.f_contiguous else self, aspect='auto', interpolation='nearest')
        pyplot.clim(minval, maxval)
        pyplot.xlabel("Signal")
        pyplot.ylabel('Detector number')
        if title is not None:
            pyplot.title(title)
        pyplot.colorbar()
        pyplot.draw()

    def __str__(self):
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

    @staticmethod
    def _validate_shape(shape, nsamples):
        if not isinstance(shape[-1], tuple) and not isinstance(shape[-1], list) and not isinstance(shape[-1], numpy.ndarray):
            if nsamples is not None and shape[-1] != numpy.sum(nsamples):
                raise ValueError("The total number of samples '"+str(shape[-1])+"' is not equal to the number of samples specified by the keyword nsamples '"+str(nsamples)+"="+str(numpy.sum(nsamples))+"'.")
            if nsamples is None:
                nsamples = (shape[-1],)
            return tuple(shape), tuple(nsamples)
        if nsamples is not None and tuple(shape[-1]) != tuple(nsamples):
            raise ValueError("The slices specified by the shape '"+str(shape)+"' are not compatible with these specified by the keyword nsamples '"+str(nsamples)+"'.")
        shape_ = list(shape[0:-1])
        shape_.append(numpy.sum(shape[-1]))
        return tuple(shape_), tuple(shape[-1])



#-------------------------------------------------------------------------------
#
# Linear operator
#
#-------------------------------------------------------------------------------


class LeastSquareMatvec(object):
    def __init__(self, model, unpacking=None):
        #cgs solves for vectors, so we need to reshape them into something (a map) that the model can handle
        if unpacking is None:
            unpacking = Reshaping(numpy.product(model.shapein), model.shapein)
        self.unpacking = unpacking
        self.criterion = unpacking.T * model.T * model * unpacking
    def __call__(self, xvec):
        xout = self.criterion(xvec, reusein=False)
        return xout

class PcgCallback():
    def __init__(self):
        self.count = 1
    def __call__(self, x):
        print 'PCG Iteration '+str(self.count)
        #if numpy.any(numpy.isnan(x)): print '... NaN in there.'
        self.count += 1

class RegularizedLeastSquareMatvec(LeastSquareMatvec):
    def __init__(self, model, unpacking=None, hyper=1.):
        if unpacking is None:
            unpacking = Reshaping(numpy.product(model.shapein), model.shapein)
        self.unpacking = unpacking
        # define smooth prior
        dX = DiscreteDifference(axis=1)
        dY = DiscreteDifference(axis=0)
        self.criterion = self.unpacking.T * (model.T * model + hyper * (dX.T * dX + dY.T * dY)) * self.unpacking

    def __call__(self, xvec):
        return self.criterion(xvec)


#-------------------------------------------------------------------------------
#
# Miscellaneous routines
#
#-------------------------------------------------------------------------------


def mapper_naive(tod, model, unit=None):
    """
    Returns a naive map, i.e.: model.transpose(tod) / model.transpose(1)
    """
    mymap = model.transpose(tod)   
    unity = tod.copy()
    unity[:] = 1.
    unity.unit = ''
    map_weights = model.transpose(unity, reusein=True)
    mymap /= map_weights
    mymap.coverage = map_weights
   
    if unit is None:
        return mymap

    newunit = Quantity(1., unit)
    newunit_si = newunit.unit_si._unit
    if 'pixel' in newunit_si and newunit_si['pixel'] == -1:
        h = mymap.header
        try:
            cd = numpy.array([ [h['cd1_1'],h['cd1_2']], [h['cd2_1'],h['cd2_2']] ])
        except:
            raise KeyError('The pixel size cannot be determined from the map header: The CD matrix is missing.')
        pixel_area = Quantity(abs(numpy.linalg.det(cd)), 'deg^2/pixel')
        mymap *= pixel_area
       
    mymap.unit = newunit._unit
       
    return mymap


#-------------------------------------------------------------------------------


def mapper_ls(tod, model, noise=Identity(), tol=1.e-6, maxiter=300):
    from scipy.sparse import dia_matrix
    from scipy.sparse.linalg import LinearOperator, cgs

    matvec = LeastSquareMatvec(noise * model)
    shape = 2*(numpy.product(model.shapein),)
    operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)

    rhs = matvec.unpacking.transpose(model.transpose(noise.transpose(noise.direct(tod))))
    x0 = matvec.unpacking.transpose(mapper_naive(tod, model))
    x0[numpy.isnan(x0)] = 0.
    solution, nit = cgs(operator, rhs, x0=x0, tol=tol, maxiter=maxiter, callback=PcgCallback())
    return Map(matvec.unpacking.direct(solution))


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, noise=Identity(), hyper=1e1, tol=1.e-6, maxiter=300):
    #from scipy.sparse import dia_matrix
    from scipy.sparse.linalg import LinearOperator, cgs

    if numpy.any(numpy.isnan(tod)) or numpy.any(numpy.isinf(tod)): raise ValueError('Input Tod contains not finite values.')

    matvec = RegularizedLeastSquareMatvec(noise * model, hyper=hyper)
    shape = 2*(numpy.product(model.shapein),)
    operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)

    rhs = matvec.unpacking.transpose(model.transpose(noise.transpose(noise.direct(tod))))
    if numpy.any(numpy.isnan(rhs)) or numpy.any(numpy.isinf(rhs)): raise ValueError('RHS contains not finite values.')
    x0 = matvec.unpacking.transpose(mapper_naive(tod, model))
    x0[numpy.isnan(x0)] = 0.
    solution, nit = cgs(operator, rhs, x0=x0, tol=tol, maxiter=maxiter, callback=PcgCallback())
    return Map(matvec.unpacking.direct(solution))


#-------------------------------------------------------------------------------


def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the standard deviation to the mean.
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    mask = tod.mask.copy()
    tmf.deglitch_l2b_std(projection._pmatrix, nx, ny, tod.T, mask.T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def deglitch_l2mad(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the MAD (median absolute deviation to the median)
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    mask = tod.mask.copy()
    tmf.deglitch_l2b_mad(projection._pmatrix, nx, ny, tod.T, mask.T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def filter_median(tod, length=10):
    """
    Median filtering, O(1) in window length
    """
    filtered = tod.copy()
    status = tmf.filter_median(filtered.T, length, numpy.array(tod.nsamples))
    if status != 0:
        raise RuntimeError()
    return filtered


#-------------------------------------------------------------------------------


def filter_polynomial(tod, degree):
    """
    Filter by subtracting a fitted polynomial of arbitrary degree.
    """
    filtered = tod.copy()
    dest = 0
    for islice in range(len(tod.nsamples)):
        nsamples = tod.nsamples[islice]
        x = numpy.arange(nsamples)
        slope = scipy.polyfit(x, tod[:,dest:dest+nsamples].T, deg=degree)

        for idetector in range(tod.shape[0]):
            filtered[idetector,dest:dest+nsamples] -= scipy.polyval(slope[:,idetector], x)
       
        dest = dest + nsamples

    return filtered


#-------------------------------------------------------------------------------


def str2fitsheader(string):
    """
    Convert a string into a pyfits.Header object
    All cards are extracted from the input string until the END keyword is reached.
    """
    header = pyfits.Header()
    cards = header.ascardlist()
    iline = 0
    while (iline*80 < len(string)):
        line = string[iline*80:(iline+1)*80]
        if line[0:3] == 'END': break
        cards.append(pyfits.Card().fromstring(line))
        iline += 1
    return header
   

#-------------------------------------------------------------------------------


def any_neq(a,b, precision):
    mask = numpy.isnan(a)
    if numpy.any(mask != numpy.isnan(b)):
        return True
    return numpy.any((abs(a-b) > 10.**(-precision) * abs(a))[numpy.isnan(a).__neg__()])
