class ValidationError(Exception): pass

try:
    import fftw3
except:
    print('Warning: Library PyFFTW3 is not installed.')

import multiprocessing
import numpy
import scipy.sparse.linalg
import tamasisfortran as tmf

from config import __verbose__, get_default_dtype, get_default_dtype_float, get_default_dtype_complex
from datatypes import Map, Tod, combine_sliced_shape, flatten_sliced_shape, validate_sliced_shape
from mpi4py import MPI
from numpyutils import _my_isscalar
from processing import interpolate_linear

__all__ = ['AcquisitionModel', 'AcquisitionModelTranspose', 'Composition', 
           'Addition', 'Square', 'Symmetric', 'Diagonal', 'Scalar', 'Identity',
           'DiscreteDifference', 'Projection', 'CompressionAverage', 'DownSampling', 
           'Masking', 'Unpacking', 'Unpacking2', 'Reshaping', 'Padding', 'Fft', 
           'FftHalfComplex', 'InvNtt', 'InterpolationLinear', 'AllReduce', 'CircularShift', 
           'ValidationError', 'asacquisitionmodel']


class ValidationError(Exception): pass

class AcquisitionModel(object):
    """Abstract class representing an instrument acquisition model.

    The acquisition model M is such as
         y = M.x,
    where x is the map and y is the instrumental response. An acquisition model can be the
    combination of several submodels describing various parts of the instrumental chain:
         M = M3 * M2 * M1 ...
    Two methods are provided to the user:
         y = M.direct(x)
         x = M.transpose(y)
    they return M.x (submodels are applied in right-to-left order) and transpose(M).y.
    Author: P. Chanial
    """

    def __init__(self, description=None):
        if description is None:
            description = self.__class__.__name__
        self.description = description
        # input and output shape are unconstrained
        self._shapein          = None
        self._shapeout         = None
        # stores the input of the transpose model. Its memory allocation is re-used as the output of the direct model
        self._output_direct    = None
        # stores the input of the direct model. Its memory allocation is re-used as the output of the transpose model
        self._output_transpose = None
        # dtype
        self._dtype = None
        # Acquisition model components
        self.blocks = []
    
    def __call__(self, data, reusein=False, reuseout=False):
        return self.direct(data, reusein=reusein, reuseout=reuseout)

    def direct(self, data, reusein=False, reuseout=False):
        raise NotImplementedError()

    def transpose(self, data, reusein=False, reuseout=False):
        raise NotImplementedError()

    @property
    def T(self):
        return AcquisitionModelTranspose(self)

    @property
    def shape(self):
        return (self.shapeout, self.shapein)

    @property
    def shapein(self):
        return self._shapein

    @shapein.setter
    def shapein(self, shape):
        self._shapein = validate_sliced_shape(shape)

    @property
    def shapeout(self):
        return self._shapeout

    @shapeout.setter
    def shapeout(self, shape):
        self._shapeout = validate_sliced_shape(shape)

    @property
    def dtype(self):
        if self._dtype is not None:
            return self._dtype
        for block in self.blocks:
            if block.dtype.type in (numpy.complex64, numpy.complex128, numpy.complex256):
                return block.dtype
        return get_default_dtype_float()

    @dtype.setter
    def dtype(self, dtype):
        if len(self.blocks) != 0:
            raise ValueError('It is not possible to set the dtype of a composite AcquisitionModel.')
        self._dtype = numpy.dtype(dtype)

    def aslinearoperator(self, shape=None, packing=None, unpacking=None):
        """Returns the AcquisitionModel as a LinearOperator"""

        if unpacking is not None:
            if packing is None:
                packing = unpacking.T
                model = packing * self * unpacking
        elif packing is not None:
            unpacking = packing.T
            model = packing * self * unpacking
        else:
            packing = Identity()
            unpacking = Identity()
            model = self
            
        shapein  = numpy.product(flatten_sliced_shape(model.shapein))
        shapeout = numpy.product(flatten_sliced_shape(model.shapeout))
        if shape is None:
            if (model.shapein is None or model.shapeout is None):
                raise ValidationError('You must specify the shape of the linear operator: the AcquisitionModel is not constrained.')
            shape = (shapeout, shapein)
        else:
            if shapein is not None and shapein != shape[1] or \
               shapeout is not None and shapeout != shape[0]:
                raise ValidationError("The input shape '" + str(shape) + 
                      "' is incompatible with the AcquisitionModel '" + str((shapein,shapeout)) + "'.")

        # ensure that input and output of model are vectors
        if model.shapein is not None and model.shapein != (shape[1],):
            reshape = Reshaping(shape[1], model.shapein)
            unpacking = unpacking * reshape
            model = model * reshape
        if model.shapeout is not None and model.shapeout != (shape[0],):
            reshape = Reshaping(model.shapeout, shape[0])
            packing = reshape * packing
            model = reshape * model

        operator = scipy.sparse.linalg.interface.LinearOperator(shape, matvec=model.direct, rmatvec=model.transpose, dtype=self.dtype)
        operator.packing   = packing
        operator.unpacking = unpacking
        operator.model = model
        return operator
        
    def validate_shape(self, shapein):
        if shapein is None:
            return self.shapein
        if self.shapein is None:
            return shapein
        if flatten_sliced_shape(shapein) != flatten_sliced_shape(self.shapein):
            raise ValidationError('The input of '+self.__class__.__name__+' has incompatible shape '+str(shapein)+' instead of '+str(self.shapein)+'.')
        return self.shapein

    def validate_shapein(self, shapein):
        if shapein is None or self.shapein is None:
            if self.shapeout is not None:
               return self.shapeout
            else:
               return shapein
        if flatten_sliced_shape(shapein) != flatten_sliced_shape(self.shapein):
            raise ValidationError('The input of '+self.__class__.__name__+' has an incompatible shape '+str(shapein)+' instead of '+str(self.shapein)+'.')
        return self.shapeout

    def validate_shapeout(self, shapeout):
        if shapeout is None or self.shapeout is None:
            if self.shapein is not None:
                return self.shapein
            else:
                return shapeout
        if flatten_sliced_shape(shapeout) != flatten_sliced_shape(self.shapeout):
            raise ValidationError("The input of '" + type(self).__name__ + ".T' has an incompatible shape " + str(shapeout) + ' instead of ' + str(self.shapeout) + '.')
        return self.shapein

    def validate_input(self, cls, data):
        data = numpy.asanyarray(data)
        if not _is_scientific_dtype(data.dtype):
            raise TypeError("The input of '" + type(self).__name__ + "' has a non-numeric type '" + str(data.dtype) + "'.")
        data = numpy.asanyarray(data, _get_dtype(self.dtype, data.dtype))
        if not isinstance(data, cls):
            data = data.view(cls)
        shapeout = self.validate_shape(data.shape)
        return data

    def validate_input_direct(self, cls, data, reusein):
        data = numpy.asanyarray(data)
        if not _is_scientific_dtype(data.dtype):
            raise TypeError("The input of '" + type(self).__name__ + "' has a non-numeric type '" + str(data.dtype) + "'.")
        data = numpy.asanyarray(data, _get_dtype(self.dtype, data.dtype))
        if not isinstance(data, cls):
            data = data.view(cls)
        shapeout = self.validate_shapein(validate_sliced_shape(data.shape, getattr(data, 'nsamples', None)))
        if reusein: self._output_transpose = data
        return data, shapeout

    def validate_input_transpose(self, cls, data, reusein):
        data = numpy.asanyarray(data)
        if not _is_scientific_dtype(data.dtype):
            raise TypeError("The input of '" + type(self).__name__ + ".T' has a non-numeric type '" + str(data.dtype) + "'.")
        data = numpy.asanyarray(data, _get_dtype(self.dtype, data.dtype))
        if not isinstance(data, cls):
            data = data.view(cls)
        shapein = self.validate_shapeout(validate_sliced_shape(data.shape, getattr(data, 'nsamples', None)))
        if reusein: self._output_direct = data
        return data, shapein

    def validate_output(self, array, reusein):
        """
        Allocate memory for the output of the acquisition models for which
        input and output have the same shape. These are the models for which
        the direct and transpose routines only have the reusein keyword.
        """
        if reusein:
            return array
        if __verbose__: print('Info: Allocating ' + str(array.dtype.itemsize * numpy.product(array.shape)/2.**20) + ' MiB in ' + type(self).__name__ + '.')
        return array.copy('a')

    def validate_output_direct(self, cls, shapeout, dtype, reuseout, **options):
        """
        Allocate memory for the output of the direct operation.
        **options should not depend on the input array, since it might
        only be used for the first call and it would not be propagated
        during the next calls.
        """
        if shapeout is None:
            raise ValueError('The shape of the output of ' + type(self).__name__+' is not known.')
        shapeout_flat = flatten_sliced_shape(shapeout)
        if self._output_direct is None or shapeout_flat != self._output_direct.shape or self._output_direct.dtype != dtype:
            if __verbose__: print('Info: Allocating '+str(dtype.itemsize*numpy.product(shapeout_flat)/2.**20)+' MiB for the output of ' + type(self).__name__+'.')
            if cls == numpy.ndarray:
                self._output_direct = numpy.empty(shapeout, dtype, **options)
            else:
                self._output_direct = cls.empty(shapeout_flat, dtype=dtype, **options)
        output = self._output_direct
        if not reuseout:
            self._output_direct = None
        return output

    def validate_output_transpose(self, cls, shapein, dtype, reuseout, **options):
        """
        Allocate memory for the output of the transpose operation.
        **options should not depend on the input array
        """
        if shapein is None:
            raise ValueError('The shape of the input of ' + type(self).__name__+' is not known.')
        shapein_flat = flatten_sliced_shape(shapein)
        if self._output_transpose is None or shapein_flat != self._output_transpose.shape or self._output_transpose.dtype != dtype:
            if __verbose__: print('Info: Allocating '+str(dtype.itemsize*numpy.product(shapein_flat)/2.**20)+' MiB for the output of ' + type(self).__name__+'.T.')
            if cls == numpy.ndarray:
                self._output_transpose = numpy.empty(shapein, dtype, **options)
            else:
                self._output_transpose = cls.empty(shapein_flat, dtype=dtype, **options)
        output = self._output_transpose
        if not reuseout:
            self._output_transpose = None
        return output

    def __mul__(self, other):
        """
        Returns the composition of two Acquisition models.
        If the operands already are the result of composition, by associativity, a flattened list is created
        to make use of the AcquisitionModel's caching system.
        """
        if isinstance(other, numpy.ndarray):
            return self(other)
        return Composition([self, other])

    def __rmul__(self, other):
        if not _my_isscalar(other):
            raise NotImplementedError("It is not possible to compose '"+str(type(other))+"' with an AcquisitionModel.")
        return Composition([other, self])

    def __imul__(self, other):
        _toacquisitionmodel(self, Composition)
        self *= other
        return self

    def __add__(self, other):
        return Addition([self, other])

    def __radd__(self, other):
        return Addition([other, self])

    def __iadd__(self, other):
        _toacquisitionmodel(self, Addition)
        self += other
        return self

    def __sub__(self, other):
        return Addition([self, -other])

    def __rsub__(self, other):
        return Addition([other, -self])

    def __isub__(self, other):
        _toacquisitionmodel(self, Addition)
        self -= other
        return self

    def __neg__(self):
        return Scalar(-1.) * self

    def __str__(self):
        result = self.description
        if type(self) == Identity:
            result += ' (Identity)'
        if self.shapein is not None or self.shapeout is not None:
            result += ' [input:'
            result += 'unconstrained' if self.shapein is None else str(self.shapein).replace(' ','')
            result += ', output:'
            result += 'unconstrained' if self.shapeout is None else str(self.shapeout).replace(' ','')
            result += ']'
        return result


#-------------------------------------------------------------------------------


class AcquisitionModelTranspose(AcquisitionModel):

    def __init__(self, model):
        self.model = model
        AcquisitionModel.__init__(self, model.description)

    def direct(self, data, reusein=False, reuseout=False):
        return self.model.transpose(data, reusein=reusein, reuseout=reuseout)

    def transpose(self, data, reusein=False, reuseout=False):
        return self.model.direct(data, reusein=reusein, reuseout=reuseout)

    @property
    def T(self):
        return self.model

    @property
    def description(self):
        return self.model.description+'.T'

    @description.setter
    def description(self, value):
        self.model.description = value

    @property
    def shapein(self):
        return self.model.shapeout

    @shapein.setter
    def shapein(self, value):
        self.model.shapeout = value

    @property
    def shapeout(self):
        return self.model.shapein

    @shapeout.setter
    def shapeout(self, value):
        self.model.shapein = value

    @property
    def _output_direct(self):
        return self.model._output_transpose

    @_output_direct.setter
    def _output_direct(self, value):
        self.model._output_transpose = value

    @property
    def _output_transpose(self):
        return self.model._output_direct

    @_output_transpose.setter
    def _output_transpose(self, value):
        self.model._output_direct = value

    @property
    def blocks(self):
        if len(self.model.blocks) != 0:
            raise NotImplementedError('This case should not happen.')
        return self.model.blocks

    @blocks.setter
    def blocks(self, value):
        pass

    def validate_shapein(self, shapein):
        return self.model.validate_shapeout(shapein)

    def validate_shapeout(self, shapeout):
        return self.model.validate_shapein(shapeout)

    def validate_input_direct(self, cls, data, reusein):
        raise NotImplementedError()

    def validate_input_transpose(self, cls, data, reusein):
        raise NotImplementedError()

    def validate_output_direct(self, cls, shapeout, dtype, reuseout, **options):
        raise NotImplementedError()

    def validate_output_transpose(self, cls, shapeout, dtype, reuseout, **options):
        raise NotImplementedError()


#-------------------------------------------------------------------------------


class Addition(AcquisitionModel):
    """
    This class adds two AcquisitionModel together.
    """
    def __init__(self, models, description=None):
        AcquisitionModel.__init__(self, description)
        if len(models) == 0:
            models = [models]
        for model in models:
            model = asacquisitionmodel(model)
            if isinstance(model, Addition):
                for block in model:
                    self.blocks.append(block)
            else:
                self.blocks.append(model)
        shapein = self.shapein
        shapeout = self.shapeout
   
    def direct(self, data, reusein=False, reuseout=False):
        data, shapeout = self.validate_input_direct(numpy.ndarray, data, reusein)
        self._output_direct = self.blocks[0].direct(data, reusein=False, reuseout=False)
        output = self._output_direct
        if not reuseout:
            self._output_direct = None
        for i, model in enumerate(self.blocks[1:]):
            reusein_ = reusein and i == len(self.blocks)-1
            output += model(data, reusein=reusein_, reuseout=True)
        return output

    def transpose(self, data, reusein=False, reuseout=False):
        data, shapein = self.validate_input_transpose(numpy.ndarray, data, reusein)
        self._output_transpose = self.blocks[0].transpose(data, reusein=False, reuseout=False)
        output = self._output_transpose
        if not reuseout:
            self._output_transpose = None
        for i, model in enumerate(self.blocks[1:]):
            reusein_ = reusein and i == len(self.blocks)-1
            output += self.blocks[i].transpose(data, reusein=reusein_, reuseout=True)
        return output

    @property
    def T(self):
        return Addition([model.T for model in self.blocks])

    @property
    def shapein(self):
        shapein = None
        for model in reversed(self.blocks):
            shapein_ = model.validate_shapeout(None)
            if shapein_ is None:
                continue
            if shapein is None:
                shapein = shapein_
                continue
            if flatten_sliced_shape(shapein) != flatten_sliced_shape(shapein_):
                raise ValidationError("Incompatible shape in operands: '" + str(shapein) +"' and '" + str(shapein_) + "'.")
        return shapein

    @property
    def shapeout(self):
        shapeout = None
        for model in reversed(self.blocks):
            shapeout_ = model.validate_shapein(None)
            if shapeout_ is None:
                continue
            if shapeout is None:
                shapeout = shapeout_
                continue
            if flatten_sliced_shape(shapeout) != flatten_sliced_shape(shapeout_):
                raise ValidationError("Incompatible shape in operands: '" + str(shapeout) +"' and '" + str(shapeout_) + "'.")
        return shapeout

    def __iadd__(self, other):
        oldblocks = self.blocks
        self.blocks.append(asacquisitionmodel(other))
        try:
            shapein = self.shapein
            shapeout = self.shapeout
        except ValidationError, errmsg:
            self.blocks = oldblocks
            raise ValidationError(errmsg)
        return self

    def __isub__(self, other):
        return self.__iadd__(-other)

    def __str__(self):
        result = super(self.__class__, self).__str__() + ':'
        components = []
        for block in self:
            components.extend(str(block).split('\n'))
        result += '\n    '+'\n    '.join(components)
        return result


#-------------------------------------------------------------------------------


class Composition(AcquisitionModel):
    """ This class represents the composition of AcquisitionModel instances."""
    def __init__(self, models, description=None):
        AcquisitionModel.__init__(self, description)
        if len(models) == 0:
            models = [models]
        for model in reversed(models):
            model = asacquisitionmodel(model)
            if isinstance(model, Composition):
                for block in model.blocks:
                    self.blocks.append(block)
            else:
                self.blocks.append(model)
        shapein = self.shapein
        shapeout = self.shapeout

    def direct(self, data, reusein=False, reuseout=False):
        for i, model in enumerate(self.blocks):
            reusein_  = reusein  or i != 0
            reuseout_ = reuseout or i != len(self.blocks)-1
            data = model.direct(data, reusein=reusein_, reuseout=reuseout_)
        return data

    def transpose(self, data, reusein=False, reuseout=False):
        for i, model in enumerate(reversed(self.blocks)):
            reusein_  = reusein  or i != 0
            reuseout_ = reuseout or i != len(self.blocks)-1
            data = model.transpose(data, reusein=reusein_, reuseout=reuseout_)
        return data

    @property
    def T(self):
        return Composition([model.T for model in self.blocks])

    @property
    def shapein(self):
        shapein = None
        for model in reversed(self.blocks):
            shapein = model.validate_shapeout(shapein)
        return shapein

    @property
    def shapeout(self):
        shapeout = None
        for model in self.blocks:
            shapeout = model.validate_shapein(shapeout)
        return shapeout

    def __imul__(self, other):
        oldblocks = self.blocks
        self.blocks.append(asacquisitionmodel(other))
        try:
            shapein = self.shapein
            shapeout = self.shapeout
        except ValidationError, errmsg:
            self.blocks = oldblocks
            raise ValidationError(errmsg)
        return self
        
    def __str__(self):
        result = super(self.__class__, self).__str__() + ':'
        components = []
        for block in self.blocks:
            components.extend(str(block).split('\n'))
        result += '\n    '+'\n    '.join(components)
        return result


#-------------------------------------------------------------------------------


class Square(AcquisitionModel):
    """Square operator"""

    @property
    def shapein(self):
        return self._shapein

    @shapein.setter
    def shapein(self, shape):
        self._shapein  = shape
        self._shapeout = shape

    @property
    def shapeout(self):
        return self._shapeout

    @shapeout.setter
    def shapeout(self, shape):
        self._shapein  = shape
        self._shapeout = shape


#-------------------------------------------------------------------------------


class Symmetric(Square):
    """Symmetric operator"""

    def transpose(self, data, reusein=False, reuseout=False):
        return self.direct(data, reusein=reusein, reuseout=reuseout)

    @property
    def T(self):
        return self


#-------------------------------------------------------------------------------


class Diagonal(Symmetric):
    """Diagonal operator"""
    def __init__(self, diagonal, nsamples=None, description=None):
        AcquisitionModel.__init__(self, description)
        self.diagonal = numpy.asarray(diagonal, dtype=get_default_dtype(diagonal))
        if not _my_isscalar(self.diagonal):
            self.shapein = validate_sliced_shape(self.diagonal.shape, nsamples)

    def direct(self, data, reusein=False, reuseout=False):
        data = self.validate_input(numpy.ndarray, data)
        output = self.validate_output(data, reusein and reuseout)
        output *= self.diagonal
        return output

    @property
    def dtype(self):
        return self.diagonal.dtype


#-------------------------------------------------------------------------------


class DiscreteDifference(AcquisitionModel):
    """Calculate the nth order discrete difference along given axis."""
    def __init__(self, n=1, axis=0, description=None):
        AcquisitionModel.__init__(self, description)
        if n != 1:
            raise NotImplementedError('DiscreteDifference is not implemented for n > 1.')
        self.n = n
        self.axis = axis

    def direct(self, data, reusein=False, reuseout=False):
        data = self.validate_input(numpy.ndarray, data)
        return numpy.diff(data, self.n, self.axis)

    def transpose(self, data, reusein=False, reuseout=False):
        data = self.validate_input(numpy.ndarray, data)
        shape = list(data.shape)
        shape[self.axis] = 1
        border = numpy.zeros(shape)
        out = numpy.concatenate((border, data, border), axis=self.axis)
        return - numpy.diff(out, axis=self.axis)

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        shapeout = list(shapein)
        shapeout[self.axis] -= self.n
        return tuple(shapeout)

    def validate_shapeout(self, shapeout):
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
        - pmatrix: transparent view of the pointing matrix
        - _pmatrix: opaque representation of the pointing matrix
        - npixels_per_sample: maximum number of sky map pixels that can be intercepted by a detector
    Author: P. Chanial
    """
    def __init__(self, observation, method=None, header=None, resolution=None, npixels_per_sample=None, oversampling=True, description=None):
        AcquisitionModel.__init__(self, description)

        self._pmatrix, self.header, ndetectors, nsamples, self.npixels_per_sample = \
            observation.get_pointing_matrix(header, resolution, npixels_per_sample, method=method, oversampling=oversampling)
        self.pmatrix = self._pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]).view(numpy.recarray)
        self.pmatrix.resize((ndetectors, numpy.sum(nsamples), self.npixels_per_sample))
        self.shapein = tuple([self.header['naxis'+str(i+1)] for i in reversed(range(self.header['naxis']))])
        self.shapeout = combine_sliced_shape(ndetectors, nsamples)

    def direct(self, map2d, reusein=False, reuseout=False):
        map2d, shapeout = self.validate_input_direct(Map, map2d, reusein)
        output = self.validate_output_direct(Tod, shapeout, map2d.dtype, reuseout)
        output.nsamples = self.shapeout[-1]
        output._unit = map2d._unit
        tmf.pointing_matrix_direct(self._pmatrix, map2d.T, output.T, self.npixels_per_sample)
        return output

    def transpose(self, signal, reusein=False, reuseout=False):
        signal, shapein = self.validate_input_transpose(Tod, signal, reusein)
        output = self.validate_output_transpose(Map, shapein, signal.dtype, reuseout)
        output.header = self.header
        output._unit = signal._unit
        tmf.pointing_matrix_transpose(self._pmatrix, signal.T, output.T,  self.npixels_per_sample)
        return output

    def get_ptp(self):
        ndetectors = self.shapeout[0]
        nsamples = numpy.sum(self.shapeout[1])
        npixels = numpy.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample, nsamples, ndetectors, npixels).T


#-------------------------------------------------------------------------------


class Compression(AcquisitionModel):
    """
    Abstract class for compressing the input signal.
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description):
        AcquisitionModel.__init__(self, description)
        if hasattr(compression_factor, 'compression_factor'):
            compression_factor = compression_factor.compression_factor
        elif hasattr(compression_factor, 'slice'):
            compression_factor = compression_factor.slice.compression_factor
        if _my_isscalar(compression_factor):
            self.factor = int(compression_factor)
        else:
            self.factor = tuple(compression_factor)

    def compression_direct(signal, output, compression_factor):
        raise NotImplementedError()

    def compression_transpose(signal, output, compression_factor):
        raise NotImplementedError()

    def direct(self, signal, reusein=False, reuseout=False):
        if numpy.all(numpy.array(self.factor) == 1):
            if not reusein: return signal.copy('a')
            return signal
        signal, shapeout = self.validate_input_direct(Tod, signal, reusein)
        output = self.validate_output_direct(Tod, shapeout, signal.dtype, reuseout)
        output.nsamples = tuple(numpy.divide(signal.nsamples, self.factor))
        self.compression_direct(signal, output, self.factor)
        return output

    def transpose(self, compressed, reusein=False, reuseout=False):
        if numpy.all(numpy.array(self.factor) == 1):
            if not reusein: return compressed.copy('a')
            return compressed
        compressed, shapein = self.validate_input_transpose(Tod, compressed, reusein)
        output = self.validate_output_transpose(Tod, shapein, compressed.dtype, reuseout)
        output.nsamples = tuple(numpy.multiply(compressed.nsamples, self.factor))
        self.compression_transpose(compressed, output, self.factor)
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        if numpy.any(numpy.array(shapein[-1]) % self.factor != 0):
            raise ValidationError('The input timeline size ('+str(shapein[-1])+') is not an integer times the compression factor ('+str(self.factor)+').')
        return combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) / self.factor)

    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return
        return combine_sliced_shape(shapeout[0:-1], numpy.array(shapeout[-1]) * self.factor)

    def __str__(self):
        return super(Compression, self).__str__()+' (x'+str(self.factor)+')'
       

#-------------------------------------------------------------------------------


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, compression_factor, description)

    def compression_direct(self, signal, output, compression_factor):
        tmf.compression_average_direct(signal.T, output.T, numpy.resize(compression_factor, [len(signal.nsamples)]))
        
    def compression_transpose(self, signal, output, compression_factor):
        tmf.compression_average_transpose(signal.T, output.T, numpy.resize(compression_factor, [len(signal.nsamples)]))
        

#-------------------------------------------------------------------------------


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number specified by the compression factor
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, compression_factor, description)

    def compression_direct(self, signal, output, compression_factor):
        tmf.downsampling_direct(signal.T, output.T, numpy.resize(compression_factor, [len(signal.nsamples)]))
        
    def compression_transpose(self, signal, output, compression_factor):
        tmf.downsampling_transpose(signal.T, output.T, numpy.resize(compression_factor, [len(signal.nsamples)]))
      

#-------------------------------------------------------------------------------


class Identity(Symmetric):
    """
    Identity class.
    """
    def __init__(self, description=None):
        AcquisitionModel.__init__(self, description)

    def direct(self, data, reusein=False, reuseout=False):
        data = self.validate_input(numpy.ndarray, data)
        return self.validate_output(data, reusein and reuseout)
       

#-------------------------------------------------------------------------------


class Scalar(Diagonal):
    """
    This class represents a scalar operator.
    """
    def __init__(self, value, description=None):
        if not numpy.iscomplex(value):
            value = numpy.real(value)
        Diagonal.__init__(self, value, description)
       
    def __str__(self):
        return super(self.__class__, self).__str__()+' (' + str(self.diagonal) + ')'
    

#-------------------------------------------------------------------------------


class Masking(Symmetric):
    """
    Apply a mask.
    Applying a boolean (or int8) mask sets to zero values whose mask is
    True (non-null). Otherwise, the input is multiplied by the mask.
    """
    def __init__(self, mask, description=None):
        AcquisitionModel.__init__(self, description)
        if mask is None:
            self.mask = None
            print('Warning: input mask is None.')
        else:
            self.mask = mask
        # shapein is not set, since it would fail for Tod with more than one slice

    def direct(self, data, reusein=False, reuseout=False):
        data = self.validate_input(numpy.ndarray, data)
        output = self.validate_output(data, reusein and reuseout)
        if self.mask is None:
            return output
        if self.mask.dtype.itemsize == 1:
            status = tmf.masking(output.T, self.mask.view(numpy.int8).T)
            if status != 0: raise RuntimeError()
        else:
            output *= self.mask
        return output

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = None
            return
        mask = numpy.asanyarray(mask)
        if not _is_scientific_dtype(mask.dtype):
            raise TypeError("Invalid type for the mask: '" + str(mask.dtype) + "'.")
        if numpy.rank(mask) == 0:
            raise TypeError('The input mask should not be scalar.')
        self._mask = mask

    @property
    def dtype(self):
        if self.mask is not None and numpy.iscomplexobj(self.mask):
            return self.mask.dtype
        return get_default_dtype_float()


#-------------------------------------------------------------------------------


class Unpacking(AcquisitionModel):
    """
    Convert 1d map into 2d map, under the control of a mask (true means non-observed)
    Author: P. Chanial
    """
    def __init__(self, mask, field=0., description=None):
        AcquisitionModel.__init__(self, description)
        self.mask = numpy.array(mask, numpy.bool8, copy=False)
        self.field = field
        self.shapein  = (numpy.sum(self.mask == 0),)
        self.shapeout = self.mask.shape

    def direct(self, packed, reusein=False, reuseout=False):
        packed, shapeout = self.validate_input_direct(numpy.ndarray, packed, reusein)
        output = self.validate_output_direct(Map, shapeout, packed.dtype, reuseout)
        tmf.unpack_direct(packed.T, self.mask.view(numpy.int8).T, output.T, self.field)
        return output

    def transpose(self, unpacked, reusein=False, reuseout=False):
        unpacked, shapein = self.validate_input_transpose(Map, unpacked, reusein)
        output = self.validate_output_transpose(numpy.ndarray, shapein, unpacked.dtype, reuseout)
        tmf.unpack_transpose(unpacked.T, self.mask.view(numpy.int8).T, output.T)
        return output


#-------------------------------------------------------------------------------


class Unpacking2(AcquisitionModel):
    """
    Convert 1d map into 2d map, under the control of a mask (true means non-observed)
    Author: P. Chanial
    """
    def __init__(self, mask, field=0., description=None):
        AcquisitionModel.__init__(self, description)
        self.mask = numpy.array(mask, numpy.bool8, copy=False)
        self.field = field
        self.shapein  = (numpy.sum(self.mask),)
        self.shapeout = tuple(self.mask.shape)

    def direct(self, packed, reusein=False, reuseout=False):
        output = numpy.zeros(self.shapeout) + self.field
        output[self.mask] = packed
        return output

    def transpose(self, unpacked, reusein=False, reuseout=False):
        return unpacked[self.mask]


#-------------------------------------------------------------------------------


class Reshaping(AcquisitionModel):
    """
    Reshape arrays
    Author: P. Chanial
    """
    def __init__(self, shapein, shapeout, description=None):
        AcquisitionModel.__init__(self, description)
        self.shapein  = shapein
        self.shapeout = shapeout
        if numpy.product(flatten_sliced_shape(self.shapein)) != numpy.product(flatten_sliced_shape(self.shapeout)):
            raise ValueError('The number of elements of the input and output of the Reshaping operator are incompatible.')

    def direct(self, array, reusein=False, reuseout=False):
        array, junk = self.validate_input_direct(numpy.ndarray, array, False)
        output = self.validate_output(array, reusein and reuseout)
        output = _smart_reshape(output, self.shapeout)
        if type(self.shapeout[-1]) is tuple:
            output = Tod(output, nsamples=self.shapeout[-1], copy=False)
        return output

    def transpose(self, array, reusein=False, reuseout=False):
        array, junk = self.validate_input_transpose(numpy.ndarray, array, False)
        output = self.validate_output(array, reusein and reuseout)
        output = _smart_reshape(output, self.shapein)
        if type(self.shapein[-1]) is tuple:
            output = Tod(output, nsamples=self.shapein[-1], copy=False)
        return output


#-------------------------------------------------------------------------------


class Padding(AcquisitionModel):
    "Pads before and after a Tod."
    def __init__(self, left=0, right=0, value=0., description=None):
        AcquisitionModel.__init__(self, description)
        left = numpy.array(left, ndmin=1, dtype=int)
        right = numpy.array(right, ndmin=1, dtype=int)
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
        array, shapeout = self.validate_input_direct(Tod, array, reusein)
        output = self.validate_output_direct(Tod, shapeout, array.dtype, reuseout)
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
        array, shapeout = self.validate_input_transpose(Tod, array, reusein)
        output = self.validate_output_transpose(Tod, shapeout, array.dtype, reuseout)
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

    def validate_input_direct(self, cls, data, reusein):
        data, shapeout = super(Padding, self).validate_input_direct(cls, data, reusein)
        if len(self.left) != 1 and len(self.left) != len(data.nsamples):
            raise ValueError("The input Tod has a number of slices '" + str(len(data.nsamples)) +
                             "' incompatible with the specified padding.")
        return data, shapeout
       
    def validate_input_transpose(self, cls, data, reusein):
        data, shapein = super(Padding, self).validate_input_transpose(cls, data, reusein)
        if len(self.left) != 1 and len(self.left) != len(data.nsamples):
            raise ValueError("The input Tod has a number of slices '" + str(len(data.nsamples)) +
                             "' incompatible with the specified padding.")
        return data, shapein       

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        return combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) + self.left + self.right)
       
    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return None
        return combine_sliced_shape(shapeout[0:-1], numpy.array(shapeout[-1]) - self.left - self.right)
       

#-------------------------------------------------------------------------------


class CircularShift(Square):
    
    def __init__(self, n, axes=None, description=None):
        Square.__init__(self, description)
        if _my_isscalar(n):
            n = (n,)
        if axes is None:
            axes = tuple(numpy.arange(-len(n), 0))
        elif _my_isscalar(axes):
            axes = (axes,)
        self.n = tuple(map(int, n))
        self.axes = tuple(map(int, axes))

    def direct(self, array, reusein=False, reuseout=False):
        for axis, n in zip(self.axes, self.n):
            array = numpy.roll(array, -n, axis=axis)
        return array

    def transpose(self, array, reusein=False, reuseout=False):
        for axis, n in zip(self.axes, self.n):
            array = numpy.roll(array, n, axis=axis)
        return array


#-------------------------------------------------------------------------------


class Fft(Square):
    """
    Performs complex fft
    """
    def __init__(self, shape, axes=None, flags=['estimate'], description=None):
        AcquisitionModel.__init__(self, description)
        if fftw3.planning.lib_threads is None:
            nthreads = 1
        else:
            nthreads = tmf.info_nthreads()
        self.shapein = shape
        self.n = numpy.product(shape)
        self.axes = axes
        self._in  = numpy.zeros(shape, dtype=get_default_dtype_complex())
        self._out = numpy.zeros(shape, dtype=get_default_dtype_complex())
        self.forward_plan = fftw3.Plan(self._in, self._out, direction='forward', flags=flags, nthreads=nthreads)
        self.backward_plan = fftw3.Plan(self._in, self._out, direction='backward', flags=flags, nthreads=nthreads)
        self.dtype = get_default_dtype_complex()

    def direct(self, array, reusein=False, reuseout=False):
        self._in[:] = array
        fftw3.execute(self.forward_plan)
        return Map(self._out)

    def transpose(self, array, reusein=False, reuseout=False):
        self._in[:] = array
        fftw3.execute(self.backward_plan)
        return Map(self._out / self.n, copy=False)


#-------------------------------------------------------------------------------


class FftHalfComplex(Square):
    """
    Performs real-to-half-complex fft
    """
    def __init__(self, nsamples, description=None):
        AcquisitionModel.__init__(self, description)
        self.nsamples = numpy.array(nsamples, ndmin=1, dtype='int64')
        self.forward_plan = numpy.empty(self.nsamples.size, dtype='int64')
        self.backward_plan = numpy.empty(self.nsamples.size, dtype='int64')
        nsamples_max = numpy.max(self.nsamples)
        for i, n in enumerate(self.nsamples):
            tarray = numpy.empty(n, dtype=get_default_dtype_float())
            farray = numpy.empty(n, dtype=get_default_dtype_float())
            self.forward_plan[i] = fftw3.Plan(tarray, farray, direction='forward', flags=['measure'], realtypes=['halfcomplex r2c'], nthreads=1)._get_parameter()
            self.backward_plan[i] = fftw3.Plan(farray, tarray, direction='backward', flags=['measure'], realtypes=['halfcomplex c2r'], nthreads=1)._get_parameter()

    def direct(self, array, reusein=False, reuseout=False):
        array = self.validate_input(Tod, array)
        output = self.validate_output(array, reusein and reuseout)
        output_ = _smart_reshape(output, (numpy.product(array.shape[:-1]), array.shape[-1]))
        tmf.fft_plan(output_.T, self.nsamples, self.forward_plan)
        return output

    def transpose(self, array, reusein=False, reuseout=False):
        array = self.validate_input(Tod, array)
        output = self.validate_output(array, reusein and reuseout)
        output_ = _smart_reshape(output, (numpy.product(array.shape[:-1]), array.shape[-1]))
        tmf.fft_plan(output_.T, self.nsamples, self.backward_plan)
        dest = 0
        for n in self.nsamples:
            output_[:,dest:dest+n] /= n
            dest += n
        return output

    def validate_shape(self, shapein):
        if shapein is None:
            return None
        nsamples = shapein[-1]        
        if not _my_isscalar(nsamples) and tuple(nsamples) != self.nsamples:
            raise ValidationError("Invalid FFT size '"+str(nsamples)+"' instead of '"+str(self.nsamples)+"'.")
        return combine_sliced_shape(shapein[0:-1], self.nsamples)


#-------------------------------------------------------------------------------


class InvNtt(Diagonal):

    def __init__(self, nsamples, filter, description=None):
        nsamples = numpy.asarray(nsamples)
        ndetectors = filter.shape[-2]
        ncorrelations = filter.shape[-1] - 1
        nslices = nsamples.size
        if numpy.rank(filter) == 2:
            filter = numpy.resize(filter, (nslices, ndetectors, ncorrelations+1))
        tod_filter, status = tmf.fft_filter_uncorrelated(filter.T, numpy.asarray(nsamples, dtype='int32'), numpy.sum(nsamples))
        if status != 0: raise RuntimeError()
        Diagonal.__init__(self, tod_filter.T, nsamples, description)
        self.ncorrelations = ncorrelations
        self.diagonal /= numpy.max(self.diagonal)


#-------------------------------------------------------------------------------


class InterpolationLinear(Square):
    def __init__(self, mask, description=None):
        Square.__init__(self, description)
        self.mask = mask
    def direct(self, data, reusein=False, reuseout=False):
        data.mask = self.mask
        return interpolate_linear(data)
    def transpose(self, data, reusein=False, reuseout=False):
        raise NotImplementedError()


#-------------------------------------------------------------------------------


class AllReduce(Square):

    def direct(self, array, reusein=False, reuseout=False, op=MPI.SUM):
        array = self.validate_output(array, reusein and reuseout)
        array[:] = MPI.COMM_WORLD.allreduce(array, op=MPI.SUM)
        return array
    def transpose(self, array, reusein=False, reuseout=False):
        array = self.validate_output(array, reusein and reuseout)
        return array


#-------------------------------------------------------------------------------


def asacquisitionmodel(operator, description=None):
    if isinstance(operator, AcquisitionModel):
        return operator
    if _my_isscalar(operator):
        return Scalar(operator)
    if isinstance(operator, scipy.sparse.linalg.interface.LinearOperator):
        model = AcquisitionModel(description)
        model.direct    = lambda data, reusein=False, reuseout=False: operator.matvec(data)
        model.transpose = lambda data, reusein=False, reuseout=False: operator.rmatvec(data)
        model.shapein   = (operator.shape[1],)
        model.shapeout  = (operator.shape[0],)
        model.dtype     = operator.dtype
        return model
    return asacquisitionmodel(scipy.sparse.linalg.aslinearoperator(operator))


#-------------------------------------------------------------------------------


def _toacquisitionmodel(model, cls, description=None):
    import copy
    if model.__class__ == cls:
        return model
    if description is None:
        description = cls.__name__
    model2 = copy.copy(model)
    model.__class__ = cls
    model.__dict__ = {'description': description, 'blocks': [model2], '_shapein': None, '_shapeout': None, '_dtype': None}


#-------------------------------------------------------------------------------


def _get_dtype(type1, type2):
    t1 = type1.type()
    t2 = type2.type()
    t = t1 * t2
    return t.dtype


#-------------------------------------------------------------------------------


def _is_scientific_dtype(dtype):
    """Return true if the data type is """
    return issubclass(dtype.type, numpy.number) or dtype.type == numpy.bool8


#-------------------------------------------------------------------------------


def _smart_reshape(array, shape):
    curr = array
    shape = flatten_sliced_shape(shape)
    while True:
        if curr.shape == shape:
            return curr
        if curr.base is None or curr.base.dtype != array.dtype or curr.base.__class__ != array.__class__:
            return curr.reshape(shape)
        curr = curr.base
