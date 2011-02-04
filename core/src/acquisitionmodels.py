try:
    import fftw3
except:
    print('Warning: Library PyFFTW3 is not installed.')

import multiprocessing
import numpy
import scipy.signal
import scipy.sparse.linalg
import tamasisfortran as tmf

from mpi4py import MPI
from scipy.sparse.linalg.interface import LinearOperator
from . import var
from .datatypes import Map, Tod, combine_sliced_shape, flatten_sliced_shape, validate_sliced_shape
from .numpyutils import _my_isscalar
from .processing import interpolate_linear
from .quantity import Quantity, UnitError, _divide_unit, _multiply_unit
from .utils import diff, diffT, diffTdiff

__all__ = [
    'AcquisitionModel',
    'AcquisitionModelLinear',
    'AllReduce',
    'CircularShift',
    'CompressionAverage',
    'Convolution',
    'DdTdd',
    'Diagonal',
    'DiscreteDifference',
    'DownSampling',
    'Fft',
    'FftHalfComplex',
    'Identity',
    'InterpolationLinear',
    'InvNtt',
    'Masking',
    'Padding',
    'Projection',
    'Reshaping',
    'ResponseTruncatedExponential',
    'Scalar',
    'Unpacking',
    'asacquisitionmodel',
]


class ValidationError(Exception): pass

class AcquisitionModel(object):
    """Abstract class representing an instrument acquisition model.

    The response y from an input signal x by an acquisition model M is given by
         y = M.direct(x) or y = M(x)
    where x and y can be multidimensional numpy.ndarray.
    An acquisition model can be the combination of several submodels
    describing various parts of the instrumental chain:
         M = M3 * M2 * M1 ...
    """

    def __init__(self, direct=None, dtype=None, description=None, 
                 attrin=None, attrout=None, shapein=None, shapeout=None, 
                 typein=None, typeout=None, unitin=None, unitout=None):

        if direct is not None:
            if not hasattr(direct, '__call__'):
                raise TypeError('The input direct method is not callable.')
            self.direct = direct
        self.dtype = dtype
        if description is None:
            description = self.__class__.__name__
        self.description = description
        self.attrin = {} if attrin is None else attrin
        self.attrout = attrout or self.attrin
        self.shapein = validate_sliced_shape(shapein)
        self.shapeout = validate_sliced_shape(shapeout or shapein)
        self.typein = typein
        self.typeout = typeout or typein
        self.unitin = Quantity(1., unitin )._unit
        self.unitout = Quantity(1., unitout or unitin)._unit

        # store the input of the transpose model.
        # Its memory allocation is re-used as the output of the direct model
        self._cacheout = None
        # store the input of the direct model. 
        # Its memory allocation is re-used as the output of the transpose model
        self._cachein = None

    def __call__(self, input, reusein=False, reuseout=False):
        return self.direct(input, reusein, reuseout)

    def direct(self, input, reusein=False, reuseout=False):
        raise NotImplementedError()

    @property
    def shape(self):
        shape = (numpy.product(flatten_sliced_shape(self.shapeout)),
                 numpy.product(flatten_sliced_shape(self.shapein)))
        if shape[0] is None or shape[1] is None:
            return None
        return shape

    @property
    def dtype(self):
        if self._dtype is not None:
            return self._dtype
        return var.FLOAT_DTYPE

    @dtype.setter
    def dtype(self, dtype):
        self._dtype = numpy.dtype(dtype)
    
    def validate_shapein(self, shapein):
        """
        Validate input shape and return the output shape of the direct model
        """
        if shapein is None or self.shapein is None:
            if self.shapeout is not None:
               return self.shapeout
            else:
               return shapein
        if flatten_sliced_shape(shapein) != flatten_sliced_shape(self.shapein):
            raise ValidationError('The input of ' + self.description + \
                ' has an incompatible shape ' + str(shapein) + ' instead of ' + \
                str(self.shapein) + '.')
        return self.shapeout

    def validate_shapeout(self, shapeout):
        """
        Validate input shape and return the output shape of the transpose model
        """
        if shapeout is None or self.shapeout is None:
            if self.shapein is not None:
                return self.shapein
            else:
                return shapeout
        if flatten_sliced_shape(shapeout) != flatten_sliced_shape(self.shapeout):
            raise ValidationError("The input of '" + self.description + \
                ".T' has an incompatible shape " + str(shapeout) + \
                ' instead of ' + str(self.shapeout) + '.')
        return self.shapein

    def validate_input_direct(self, input, reusein, reuseout, **options):

        # validate input
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + \
                "' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        if self.typein is not None and not isinstance(input, self.typein):
            input = input.view(self.typein)
        _validate_input_unit(input, self.unitin)

        # store input for the transpose's output
        if reusein:
            self._cachein = input

        # validate output
        shape = self.validate_shapein(
            validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + \
                type(self).__name__+' is not known.')
        shape_flat = flatten_sliced_shape(shape)
        typeout = self.typeout or input.__class__
        if self._cacheout is None or shape_flat != self._cacheout.shape or \
           self._cacheout.dtype != input.dtype or typeout != type(self._cachein):
            self._cacheout = None
            if var.verbose: print('Info: Allocating '+str(input.dtype.itemsize*numpy.product(shape_flat)/2.**20)+' MiB for the output of ' + type(self).__name__+'.')
            if typeout == numpy.ndarray:
                output = numpy.empty(shape, input.dtype, **options)
            else:
                output = typeout.empty(shape, dtype=input.dtype, **options)
            if reuseout:
                self._cacheout = output
        else:
            output = self._cacheout
            if not reuseout:
                self._cacheout = None
        #XXX the cached value should not need the following lines
        if type(shape[-1]) is tuple:
            output = Tod(output, nsamples=shape[-1], copy=False)
        _propagate_attributes(input, output)
        _validate_output_unit(input, output, self.unitin, self.unitout)
        for k,v in self.attrout.items():
            setattr(output, k, v)

        return input, output

    def validate_input_transpose(self, input, reusein, reuseout, **options):

        # validate input
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + \
                ".T' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        if self.typeout is not None and not isinstance(input, self.typeout):
            input = input.view(self.typeout)
        _validate_input_unit(input, self.unitout)

        # store input for the direct's output
        if reusein:
            self._cacheout = input

        # validate output
        shape = self.validate_shapeout(
            validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + \
                self.description + ' is not known.')
        shape_flat = flatten_sliced_shape(shape)
        typeout = self.typein or input.__class__
        if self._cachein is None or shape_flat != self._cachein.shape or \
           self._cachein.dtype != input.dtype or typeout != type(self._cachein):
            self._cachein = None
            if var.verbose: print('Info: Allocating '+str(input.dtype.itemsize*numpy.product(shape_flat)/2.**20)+' MiB for the output of ' + type(self).__name__+'.T.')
            if typeout == numpy.ndarray:
                output = numpy.empty(shape_flat, input.dtype, **options)
            else:
                output = typeout.empty(shape_flat, dtype=input.dtype, **options)
            if reuseout:
                self._cachein = output
        else:
            output = self._cachein
            if not reuseout:
                self._cachein = None
        #XXX the cached value should not need the following lines
        if type(shape[-1]) is tuple:
            output = Tod(output, nsamples=shape[-1], copy=False)
        _propagate_attributes(input, output)
        _validate_output_unit(input, output, self.unitout, self.unitin)
        for k,v in self.attrin.items():
            setattr(output, k, v)

        return input, output

    def __mul__(self, other):
        if isinstance(other, numpy.ndarray):
            return self.matvec(other)
        return Composition([self, other])

    def __rmul__(self, other):
        if not _my_isscalar(other):
            raise NotImplementedError("It is not possible to multiply '" + \
                str(type(other)) + "' with an AcquisitionModel.")
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
            if self.shapein is None:
                result += 'unconstrained'
            else:
                result += str(self.shapein).replace(' ','')
            result += ', output:'
            if self.shapeout is None:
                result += 'unconstrained'
            else:
                result += str(self.shapeout).replace(' ','')
            result += ']'
        return result


#-------------------------------------------------------------------------------


class AcquisitionModelLinear(AcquisitionModel, LinearOperator):
    """Abstract class representing a linear instrument acquisition model.

    The response y from an input signal x by an acquisition model M is given by
         y = M.direct(x) or y = M(x)
    where x and y can be multidimensional numpy.ndarray.
    The transpose of the acquisition model is
        x = M.transpose(y) or M.T(y)
    This class subclasses the LinearOperator class, so it also provides the
    methods matvec and rmatvec which operate on 1d ndarray.

    An acquisition model can be the combination of several submodels
    describing various parts of the instrumental chain:
         M = M3 * M2 * M1 ...
    """
    def __init__(self, transpose=None, **kw):
        AcquisitionModel.__init__(self, **kw)
        if transpose is not None:
            if not hasattr(transpose, '__call__'):
                raise TypeError('The input transpose method is not callable.')
            self.transpose = transpose

    def transpose(self, input, reusein=False, reuseout=False):
        raise NotImplementedError()

    def matvec(self, v):
        v = v.reshape(flatten_sliced_shape(self.shapein))
        return self.direct(v).ravel()

    def rmatvec(self, v):
        v = v.reshape(flatten_sliced_shape(self.shapeout))
        return self.transpose(v).ravel()

    def dense(self):
        d = numpy.ndarray(self.shape, dtype=self.dtype)
        v = numpy.zeros(self.shape[1], dtype=var.FLOAT_DTYPE)
        for i in range(self.shape[1]):
            v[:] = 0
            v[i] = 1
            d[:,i] = self.matvec(v)
        return d

    @property
    def T(self):
        return AcquisitionModelTranspose(self)


#-------------------------------------------------------------------------------


class AcquisitionModelTranspose(AcquisitionModelLinear):

    def __init__(self, model):
        self.model = model
        AcquisitionModelLinear.__init__(self,
                                        direct=model.transpose,
                                        transpose=model.direct,
                                        dtype=model.dtype,
                                        description=model.description + '.T',
                                        attrin=model.attrout,
                                        attrout=model.attrin,
                                        shapein=model.shapeout,
                                        shapeout=model.shapein,
                                        typein=model.typeout,
                                        typeout=model.typein,
                                        unitin=model.unitout,
                                        unitout=model.unitin)

    @property
    def T(self):
        return self.model

    def validate_shapein(self, shapein):
        return self.model.validate_shapeout(shapein)

    def validate_shapeout(self, shapeout):
        return self.model.validate_shapein(shapeout)


#-------------------------------------------------------------------------------


class Composite(AcquisitionModel):
    """
    Class for grouping acquisition models
    """

    def __init__(self, models, description=None):

        self.blocks = []

        if len(models) == 0:
            models = [models]
        models = [ asacquisitionmodel(m) for m in models ]

        for model in models:
            model = asacquisitionmodel(model)
            if isinstance(model, self.__class__):
                for block in model.blocks:
                    self.blocks.append(block)
            else:
                self.blocks.append(model)

        AcquisitionModel.__init__(self, description=description)

        if all([hasattr(m, 'matvec') for m in self.blocks]):
            self.matvec = lambda v: \
                self.direct(v.reshape(flatten_sliced_shape(self.shapein))) \
                    .ravel()
        if all([hasattr(m, 'rmatvec') for m in self.blocks]):
            self.rmatvec = lambda v: \
                self.transpose(v.reshape(flatten_sliced_shape(self.shapeout))) \
                    .ravel()

    @property
    def dtype(self):
        for block in self.blocks:
            if block.dtype.type in (numpy.complex64, numpy.complex128, numpy.complex256):
                return block.dtype
        return var.FLOAT_DTYPE

    @dtype.setter
    def dtype(self, dtype):
        pass

    @property
    def T(self):
        return AcquisitionModelTranspose(self)

    @property
    def typein(self):
        for model in reversed(self.blocks):
            if model.typein is not None:
                return model.typein
        return numpy.ndarray

    @typein.setter
    def typein(self, value):
        pass

    @property
    def typeout(self):
        for model in reversed(self.blocks):
            if model.typeout is not None:
                return model.typeout
        return numpy.ndarray

    @typeout.setter
    def typeout(self, value):
        pass

    @property
    def unitin(self):
        for model in reversed(self.blocks):
            if len(model.unitin) > 0:
                return model.unitin
        return {}

    @unitin.setter
    def unitin(self, value):
        pass

    @property
    def unitout(self):
        for model in self.blocks:
            if len(model.unitout) > 0:
                return model.unitout
        return {}

    @unitout.setter
    def unitout(self, value):
        pass

    def validate_input(self, input, shape, inplace):
        input = numpy.asanyarray(input)
        if shape is not None and type(shape[-1]) is tuple:
            input = Tod(input, nsamples=shape[-1], copy=False)
        if not inplace:
            if var.verbose: print('Info: Allocating '+str(input.dtype.itemsize * numpy.product(flatten_sliced_shape(shape))/2.**20)+' MiB for the output of Addition.')
            input = input.copy()
        return input


#-------------------------------------------------------------------------------


class Addition(Composite):
    """
    Class for acquisition models addition

    If at least one of the input already is the result of an addition,
    a flattened list of operators is created by associativity, in order to
    benefit from the AcquisitionModel's caching mechanism.
    """

    def direct(self, input, reusein=False, reuseout=False):
        input = self.validate_input(input, self.shapein, False)
        output = self.blocks[0].direct(input, reusein=False, reuseout=False)
        for i, model in enumerate(self.blocks[1:]):
            reusein_ = reusein and i == len(self.blocks)-1
            output += model.direct(input, reusein=reusein_, reuseout=True)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        input = self.validate_input(input, self.shapeout, False)
        output = self.blocks[0].transpose(input, reusein=False, reuseout=False)
        for i, model in enumerate(self.blocks[1:]):
            reusein_ = reusein and i == len(self.blocks)-1
            output += model.transpose(input, reusein=reusein_, reuseout=True)
        return output

    @property
    def shapein(self):
        shapein = None
        for model in self.blocks:
            shapein_ = model.validate_shapeout(None)
            if shapein_ is None:
                continue
            if shapein is None or type(shapein_[-1]) is tuple:
                shapein = shapein_
                continue
            if flatten_sliced_shape(shapein) != flatten_sliced_shape(shapein_):
                raise ValidationError("Incompatible shape in operands: '" + \
                          str(shapein) +"' and '" + str(shapein_) + "'.")
        return shapein

    @shapein.setter
    def shapein(self, value):
        pass

    @property
    def shapeout(self):
        shapeout = None
        for model in self.blocks:
            shapeout_ = model.validate_shapein(None)
            if shapeout_ is None:
                continue
            if shapeout is None or type(shapeout_[-1]) is tuple:
                shapeout = shapeout_
                continue
            if flatten_sliced_shape(shapeout) != flatten_sliced_shape(shapeout_):
                raise ValidationError("Incompatible shape in operands: '" + \
                          str(shapeout) +"' and '" + str(shapeout_) + "'.")
        return shapeout

    @shapeout.setter
    def shapeout(self, value):
        pass

    def __iadd__(self, other):
        oldblocks = self.blocks
        self.blocks.append(asacquisitionmodel(other))
        try:
            shapein = self.shapein
            shapeout = self.shapeout
        except ValidationError as errmsg:
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


class Composition(Composite):
    """
    Class for acquisition models composition

    If at least one of the input already is the result of a composition,
    a flattened list of operators is created by associativity, in order to
    benefit from the AcquisitionModel's caching mechanism.
    """

    def direct(self, input, reusein=False, reuseout=False):
        input = self.validate_input(input, self.shapein, reusein)
        for i, model in enumerate(reversed(self.blocks)):
            input = model.direct(input, True, reuseout or i != len(self.blocks)-1)
        return input

    def transpose(self, input, reusein=False, reuseout=False):
        input = self.validate_input(input, self.shapeout, reusein)
        for i, model in enumerate(self.blocks):
            input = model.transpose(input, True, reuseout or i != len(self.blocks)-1)
        return input

    @property
    def shapein(self):
        shapeout = None
        for model in self.blocks:
            shapeout = model.validate_shapeout(shapeout)
        return shapeout

    @shapein.setter
    def shapein(self, value):
        pass

    @property
    def shapeout(self):
        shapein = None
        for model in reversed(self.blocks):
            shapein = model.validate_shapein(shapein)
        return shapein

    @shapeout.setter
    def shapeout(self, value):
        pass

    def __imul__(self, other):
        oldblocks = self.blocks
        self.blocks.append(asacquisitionmodel(other))
        try:
            shapein = self.shapein
            shapeout = self.shapeout
        except ValidationError as errmsg:
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


class SquareBase(AcquisitionModelLinear):
    """
    Square operator
    
    The input and output must have the same number of elements but not
    necessarily the same shape.
    This operator does not implement the cache mechanism, but operation on
    the input can be done inplace or on a copy.
    """

    def validate_input_direct(self, input, inplace):
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + \
                "' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))

        shape = self.validate_shapein(validate_sliced_shape(input.shape, 
                    getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + \
                      type(self).__name__+' is not known.')

        typein = self.typein
        if type(shape[-1]) is tuple and (typein is None or \
           not issubclass(typein, Tod)):
            typein = Tod
        if typein is not None and not isinstance(input, typein):
            input = input.view(typein)

        if type(shape[-1]) is tuple:
            input.nsamples = shape[-1]

        if not inplace:
            if var.verbose:
                print('Info: Allocating ' + str(input.dtype.itemsize * \
                      numpy.product(input.shape)/2.**20) + ' MiB in ' + \
                      self.description + '.')
            input = input.copy()

        for k,v in self.attrout:
            setattr(input, k, v)

        return input

    def validate_input_transpose(self, input, inplace):
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + ".T' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))

        shape = self.validate_shapeout(validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + type(self).__name__+' is not known.')

        typeout = self.typeout
        if type(shape[-1]) is tuple and (typeout is None or \
           not issubclass(typeout, Tod)):
            typeout = Tod
        if typeout is not None and not isinstance(input, typeout):
            input = input.view(typeout)

        if type(shape[-1]) is tuple:
            input.nsamples = shape[-1]

        if not inplace:
            if var.verbose:
                print('Info: Allocating ' + str(input.dtype.itemsize * \
                      numpy.product(input.shape)/2.**20) + ' MiB in ' + \
                      self.description + '.T.')
            input = input.copy()

        for k,v in self.attrin:
            setattr(input, k, v)

        return input


#-------------------------------------------------------------------------------


class Square(SquareBase):
    """
    Square operator
    
    The input and output must have the same shape
    This operator does not implement the cache mechanism, but operation on
    the input can be done inplace or on a copy.
    """

    def __init__(self, shapein=None, **kw):
        AcquisitionModelLinear.__init__(self, shapein=shapein, shapeout=shapein, **kw)

    def validate_shapeout(self, shapeout):
        return self.validate_shapein(shapeout)

    def validate_input(self, input, inplace):
        return super(Square,self).validate_input_direct(input, inplace)
    

#-------------------------------------------------------------------------------


class Symmetric(Square):
    """Symmetric operator"""

    def transpose(self, input, reusein=False, reuseout=False):
        return self.direct(input, reusein=reusein, reuseout=reuseout)

    @property
    def T(self):
        return self


#-------------------------------------------------------------------------------


class Diagonal(Symmetric):
    """Diagonal operator"""

    def __init__(self, diagonal, shapein=None, description=None):
        diagonal = numpy.asarray(diagonal, dtype=var.get_default_dtype(diagonal))
        Symmetric.__init__(self, dtype=diagonal.dtype, description=description, 
                           shapein=shapein)
        self.diagonal = diagonal

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        if output.ndim == 0:
            return output * self.diagonal
        output.T[:] *= self.diagonal.T
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return self.shapein
        if flatten_sliced_shape(shapein[0:self.diagonal.ndim]) != self.diagonal.shape:
            raise ValueError('The input has an incompatible shape ' + str(shapein) + '.')
        return shapein


#-------------------------------------------------------------------------------


class DiscreteDifference(Square):
    """Calculate the nth order discrete difference along given axis."""

    def __init__(self, n=1, axis=0, shapein=None, description=None):
        Square.__init__(self, shapein=shapein, description=description)
        self.n = n
        self.axis = axis

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        for i in range(self.n):
            diff(output, self.axis)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        for i in range(self.n):
            diffT(output, self.axis)
        return output        


class DdTdd(Symmetric):
    """Calculate operator dX.T dX along a given axis."""

    def __init__(self, shapein=None, axis=0, description=None):
        Symmetric.__init__(self, shapein=shapein, description=description)
        self.axis = axis

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        diffTdiff(output, self.axis)
        return output

   
#-------------------------------------------------------------------------------


class Projection(AcquisitionModelLinear):
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
    """

    def __init__(self, observation, method=None, header=None, resolution=None, npixels_per_sample=0, oversampling=True, description=None):

        self._pmatrix, self.header, ndetectors, nsamples, \
        self.npixels_per_sample, (unitout, unitin), (duout, duin) = \
            observation.get_pointing_matrix(header,
                                            resolution,
                                            npixels_per_sample,
                                            method=method,
                                            oversampling=oversampling)

        attrin = {'header' : self.header}
        if duin is not None:
            attrin['derived_units'] = duin
        attrout = {}
        if duout is not None:
            attrout['derived_units'] = duout
        shapein = tuple([self.header['naxis'+str(i+1)] \
                         for i in reversed(list(range(self.header['naxis'])))])
        shapeout = combine_sliced_shape(ndetectors, nsamples)
        AcquisitionModelLinear.__init__(self,
                                        attrin=attrin,
                                        attrout=attrout,
                                        shapein=shapein,
                                        shapeout=shapeout,
                                        typein=Map,
                                        typeout=Tod,
                                        unitin=unitin,
                                        unitout=unitout,
                                        description=description)

        self.pmatrix = self._pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]) \
                           .view(numpy.recarray)
        self.pmatrix.resize((ndetectors, numpy.sum(nsamples),
                             self.npixels_per_sample))

    def direct(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_direct(input, reusein, reuseout)
        tmf.pointing_matrix_direct(self._pmatrix, input.T, output.T,
                                   self.npixels_per_sample)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        tmf.pointing_matrix_transpose(self._pmatrix, input.T, output.T, 
                                      self.npixels_per_sample)
        return output

    def get_ptp(self):
        ndetectors = self.shapeout[0]
        nsamples = numpy.sum(self.shapeout[1])
        npixels = numpy.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample,
                                       nsamples, ndetectors, npixels).T


#-------------------------------------------------------------------------------


class Compression(AcquisitionModelLinear):
    """
    Abstract class for compressing the input signal.
    """

    def __init__(self, compression_factor, shapein=None, description=None):
        shapeout = self.validate_shapein(shapein)
        AcquisitionModelLinear.__init__(self, description=description,
                                        shapein=shapein, shapeout=shapeout,
                                        typein=Tod)
        if _my_isscalar(compression_factor):
            self.factor = int(compression_factor)
        else:
            self.factor = tuple(compression_factor)

    def compression_direct(input, output, compression_factor):
        raise NotImplementedError()

    def compression_transpose(input, output, compression_factor):
        raise NotImplementedError()

    def direct(self, input, reusein=False, reuseout=False):
        if numpy.all(numpy.array(self.factor) == 1):
            if not reusein: return input.copy('a')
            return input
        input, output = self.validate_input_direct(input, reusein, reuseout)
        self.compression_direct(input, output, self.factor)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        if numpy.all(numpy.array(self.factor) == 1):
            if not reusein: return input.copy('a')
            return input
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        self.compression_transpose(input, output, self.factor)
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        if numpy.any(numpy.array(shapein[-1]) % self.factor != 0):
            raise ValidationError('The input timeline size ('+str(shapein[-1])+\
                ') is not an integer times the compression factor (' + \
                str(self.factor)+').')
        return combine_sliced_shape(shapein[0:-1],
                                    numpy.array(shapein[-1]) / self.factor)

    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return None
        return combine_sliced_shape(shapeout[0:-1],
                                    numpy.array(shapeout[-1]) * self.factor)

    def __str__(self):
        return super(Compression, self).__str__()+' (x'+str(self.factor)+')'
       

#-------------------------------------------------------------------------------


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    """

    def __init__(self, compression_factor, shapein=None, description=None):
        Compression.__init__(self, compression_factor, shapein=shapein,
                             description=description)

    def compression_direct(self, input, output, compression_factor):
        tmf.compression_average_direct(input.T, output.T,
            numpy.resize(compression_factor, [len(input.nsamples)]))
        
    def compression_transpose(self, input, output, compression_factor):
        tmf.compression_average_transpose(input.T, output.T,
            numpy.resize(compression_factor, [len(input.nsamples)]))
        
        

#-------------------------------------------------------------------------------


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number
    specified by the compression factor
    """

    def __init__(self, compression_factor, shapein=None, description=None):
        Compression.__init__(self, compression_factor, shapein=shapein,
                             description=description)

    def compression_direct(self, input, output, compression_factor):
        tmf.downsampling_direct(input.T, output.T,
            numpy.resize(compression_factor, [len(input.nsamples)]))
        
    def compression_transpose(self, input, output, compression_factor):
        tmf.downsampling_transpose(input.T, output.T,
            numpy.resize(compression_factor, [len(input.nsamples)]))
      

#-------------------------------------------------------------------------------


class Identity(Symmetric):
    """
    Identity class.
    """

    def __init__(self, shapein=None, description=None):
        Symmetric.__init__(self, shapein=shapein, description=description)

    def direct(self, input, reusein=False, reuseout=False):
        return self.validate_input(input, reusein and reuseout)
       

#-------------------------------------------------------------------------------


class Scalar(Diagonal):
    """
    Class for scalar multiplication
    """

    def __init__(self, value, shapein=None, description=None):
        if not numpy.iscomplex(value):
            value = numpy.real(value)
        Diagonal.__init__(self, value, shapein=shapein, description=description)
       
    def __str__(self):
        return super(self.__class__, self).__str__()+' (' + \
               str(self.diagonal) + ')'
    

#-------------------------------------------------------------------------------


class Masking(Symmetric):
    """
    Class for mask application
    Applying a boolean (or int8) mask sets to zero values whose mask is
    True (non-null). Otherwise, the input is multiplied by the mask.
    """

    def __init__(self, mask, shapein=None, description=None):
        if mask is None:
            print('Warning: input mask is None.')
        else:
            mask = numpy.asanyarray(mask)
            if not _is_scientific_dtype(mask.dtype):
                raise TypeError("Invalid mask type: '" + str(mask.dtype) + "'.")
            if numpy.rank(mask) == 0:
                raise TypeError('The input mask should not be scalar.')
        if mask is not None and numpy.iscomplexobj(mask):
            dtype = mask.dtype
        else:
            dtype = var.FLOAT_DTYPE

        Symmetric.__init__(self, description=description, dtype=dtype,
                           shapein=shapein)
        self.mask = mask

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        if self.mask is None:
            return output
        if self.mask.dtype.itemsize == 1:
            status = tmf.masking(output.T, self.mask.view(numpy.int8).T)
            if status != 0: raise RuntimeError()
        else:
            output *= self.mask
        return output


#-------------------------------------------------------------------------------


class Unpacking(AcquisitionModelLinear):
    """
    Convert 1d map into an nd array, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, field=0., description=None):
        mask = numpy.array(mask, numpy.bool8, copy=False)
        AcquisitionModelLinear.__init__(self,
                                        shapein=numpy.sum(mask == 0),
                                        shapeout=mask.shape,
                                        description=description)
        self.mask = mask
        self.field = field

    def direct(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_direct(input, reusein, reuseout)
        tmf.unpack_direct(input.T, self.mask.view(numpy.int8).T, output.T,
                          self.field)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        tmf.unpack_transpose(input.T, self.mask.view(numpy.int8).T, output.T)
        return output


#-------------------------------------------------------------------------------


class Reshaping(SquareBase):
    """
    Reshape arrays
    """

    def __init__(self, shapein, shapeout, description=None):
        if shapein is None or shapeout is None:
            raise ValueError('The shapes are not defined.')
        if numpy.product(flatten_sliced_shape(shapein)) != \
           numpy.product(flatten_sliced_shape(shapeout)):
            raise ValueError('The number of elements of the input and output o'\
                             'f the Reshaping operator are incompatible.')
        AcquisitionModelLinear.__init__(self, description=description,
                                        shapein=shapein, shapeout=shapeout)

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input_direct(input, reusein and reuseout)
        output = _smart_reshape(output, self.shapeout)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input_transpose(input, reusein and reuseout)
        output = _smart_reshape(output, self.shapein)
        return output


#-------------------------------------------------------------------------------


class ResponseTruncatedExponential(Square):
    """
    ResponseTruncatedExponential(tau)

    Apply a truncated exponential response to the signal

    Parameters
    ==========
    
    tau: number
        Time constant divided by the signal sampling time
    """
    
    def __init__(self, tau, shapein=None, description=None):
        """
        """
        Square.__init__(self, description=description, shapein=shapein,
                        typein=Tod)
        if hasattr(tau, 'SI'):
            tau = tau.SI
            if tau.unit != '':
                raise ValueError('The time constant must be dimensionless.')
        self.tau = numpy.array(tau, dtype=var.FLOAT_DTYPE, ndmin=1)

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        tmf.convolution_trexp_direct(output.T, numpy.array(output.nsamples),
                                     self.tau)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        tmf.convolution_trexp_transpose(output.T, numpy.array(output.nsamples),
                                        self.tau)
        return output


#-------------------------------------------------------------------------------


class Padding(AcquisitionModelLinear):
    "Pads before and after a Tod."

    def __init__(self, left=0, right=0, value=0., shapein=None,
                 description=None):
        if shapein is not None:
            shapeout = self.validate_shapein(shapein)
        else:
            shapeout = None
        AcquisitionModelLinear.__init__(self, description=description,
                                        shapein=shapein, shapeout=shapeout,
                                        typein=Tod)
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
   
    def direct(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_direct(input, reusein, reuseout)
        dest = 0
        dest_padded = 0
        for islice in range(len(input.nsamples)):
            nsamples = input.nsamples[islice]
            left = self.left[islice if len(self.left) > 1 else 0]
            output[...,dest_padded:dest_padded+left] = self.value
            output[...,dest_padded+left:dest_padded+left+nsamples] = \
                input[...,dest:dest+nsamples]
            output[...,dest_padded+left+nsamples:dest_padded+output.nsamples[islice]] = self.value
            dest += nsamples
            dest_padded += output.nsamples[islice]
        return output
   
    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        dest = 0
        dest_padded = 0
        for islice in range(len(input.nsamples)):
            nsamples = output.nsamples[islice]
            left  = self.left [islice if len(self.left)  > 1 else 0]
            output[...,dest:dest+nsamples] = \
                input[...,dest_padded+left:dest_padded+left+nsamples]
            dest += nsamples
            dest_padded += input.nsamples[islice]
        return output

    def validate_input_direct(self, input, reusein, reuseout):
        input, output = super(Padding, self).validate_input_direct(input, reusein, reuseout)
        if len(self.left) != 1 and len(self.left) != len(input.nsamples):
            raise ValueError("The input Tod has a number of slices '" + \
                             str(len(input.nsamples)) + \
                             "' incompatible with the specified padding.")
        return input, output
       
    def validate_input_transpose(self, input, reusein, reuseout):
        input, output = super(Padding, self).validate_input_transpose(input, reusein, reuseout)
        if len(self.left) != 1 and len(self.left) != len(input.nsamples):
            raise ValueError("The input Tod has a number of slices '" + \
                             str(len(input.nsamples)) +
                             "' incompatible with the specified padding.")
        return input, output

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        return combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) + \
                                    self.left + self.right)
       
    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return None
        return combine_sliced_shape(shapeout[0:-1], numpy.array(shapeout[-1]) -\
                                    self.left - self.right)
       

#-------------------------------------------------------------------------------


class CircularShift(Square):
    
    def __init__(self, n, axes=None, shapein=None, description=None):
        Square.__init__(self, shapein=shapein, description=description)
        if _my_isscalar(n):
            n = (n,)
        if axes is None:
            axes = tuple(numpy.arange(-len(n), 0))
        elif _my_isscalar(axes):
            axes = (axes,)
        self.n = tuple(map(int, n))
        self.axes = tuple(map(int, axes))

    def direct(self, input, reusein=False, reuseout=False):
        for axis, n in zip(self.axes, self.n):
            input = numpy.roll(input, -n, axis=axis)
        return input

    def transpose(self, input, reusein=False, reuseout=False):
        for axis, n in zip(self.axes, self.n):
            input = numpy.roll(input, n, axis=axis)
        return input


#-------------------------------------------------------------------------------


class Fft(Square):
    """
    Performs complex fft
    """

    def __init__(self, shape, axes=None, flags=['estimate'], description=None):
        Square.__init__(self, description=description, shapein=shape,
                        dtype=var.COMPLEX_DTYPE)
        if fftw3.planning.lib_threads is None:
            nthreads = 1
        else:
            nthreads = tmf.info_nthreads()
        self.n = numpy.product(shape)
        self.axes = axes
        self._in  = numpy.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self._out = numpy.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self.forward_plan = fftw3.Plan(self._in, self._out, direction='forward',
                                       flags=flags, nthreads=nthreads)
        self.backward_plan= fftw3.Plan(self._in, self._out,direction='backward',
                                       flags=flags, nthreads=nthreads)

    def direct(self, input, reusein=False, reuseout=False):
        self._in[:] = input
        fftw3.execute(self.forward_plan)
        return Map(self._out)

    def transpose(self, input, reusein=False, reuseout=False):
        self._in[:] = input
        fftw3.execute(self.backward_plan)
        return Map(self._out / self.n, copy=False)


#-------------------------------------------------------------------------------


class FftHalfComplex(Square):
    """
    Performs real-to-half-complex fft
    """

    def __init__(self, nsamples, shapein=None, description=None):
        Square.__init__(self, description=description, shapein=shapein,
                        typein=Tod)
        self.nsamples_array = numpy.array(nsamples, ndmin=1, dtype='int64')
        self.nsamples = tuple(self.nsamples_array)
        self.nsamples_tot = numpy.sum(self.nsamples_array)
        self.forward_plan = numpy.empty(self.nsamples_array.size, dtype='int64')
        self.backward_plan = numpy.empty(self.nsamples_array.size,dtype='int64')
        for i, n in enumerate(self.nsamples):
            tarray = numpy.empty(n, dtype=var.FLOAT_DTYPE)
            farray = numpy.empty(n, dtype=var.FLOAT_DTYPE)
            self.forward_plan[i] = \
                fftw3.Plan(tarray, farray, direction='forward',
                           flags=['measure'], realtypes=['halfcomplex r2c'],
                           nthreads=1)._get_parameter()
            self.backward_plan[i] = \
                fftw3.Plan(farray, tarray, direction='backward',
                           flags=['measure'], realtypes=['halfcomplex c2r'],
                           nthreads=1)._get_parameter()

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        output_ = _smart_reshape(output, (numpy.product(input.shape[:-1]),
                                 input.shape[-1]))
        tmf.fft_plan(output_.T, self.nsamples_array, self.forward_plan)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        output_ = _smart_reshape(output, (numpy.product(input.shape[:-1]), 
                                 input.shape[-1]))
        tmf.fft_plan(output_.T, self.nsamples_array, self.backward_plan)
        dest = 0
        for n in self.nsamples:
            output_[:,dest:dest+n] /= n
            dest += n
        return output

    def validate_shapein(self, shape):
        if shape is None:
            return None
        nsamples = shape[-1]
        if nsamples != self.nsamples and nsamples != self.nsamples_tot:
            raise ValidationError("Invalid FFT size '" + str(nsamples) + \
                                  "' instead of '"+str(self.nsamples)+"'.")
        return combine_sliced_shape(shape[0:-1], self.nsamples)


#-------------------------------------------------------------------------------


class Convolution(Symmetric):

    def __init__(self, kernel, shapein=None, description=None):
        Symmetric.__init__(self, shapein=shapein, description=description)
        self.kernel = numpy.asanyarray(kernel)

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        output[:] = scipy.signal.fftconvolve(input, self.kernel, mode='same')
        return output


#-------------------------------------------------------------------------------


class InvNtt(Diagonal):

    def __init__(self, nsamples, filter, description=None):
        nsamples = numpy.asarray(nsamples)
        ndetectors = filter.shape[-2]
        ncorrelations = filter.shape[-1] - 1
        nslices = nsamples.size
        if numpy.rank(filter) == 2:
            filter = numpy.resize(filter,(nslices, ndetectors, ncorrelations+1))
        tod_filter, status = \
            tmf.fft_filter_uncorrelated(filter.T, numpy.asarray(nsamples, 
                                        dtype='int32'), numpy.sum(nsamples))
        if status != 0: raise RuntimeError()
        Diagonal.__init__(self, tod_filter.T, description=description,
                          shapein=tod_filter.T.shape)
        self.ncorrelations = ncorrelations
        self.diagonal /= var.mpi_comm.allreduce(numpy.max(self.diagonal),
                                                op=MPI.MAX)


#-------------------------------------------------------------------------------


class InterpolationLinear(Square):

    def __init__(self, mask, shapein=None, description=None):
        Square.__init__(self, attrin={'mask' : mask}, shapein=shapein,
                        description=description)
        self.mask = mask

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        return interpolate_linear(output)

    def transpose(self, input, reusein=False, reuseout=False):
        raise NotImplementedError()


#-------------------------------------------------------------------------------


class AllReduce(Square):

    def direct(self, input, reusein=False, reuseout=False, op=MPI.SUM):
        output = self.validate_input(input, reusein and reuseout)
        output[:] = var.mpi_comm.allreduce(output, op=MPI.SUM)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        return output


#-------------------------------------------------------------------------------


def asacquisitionmodel(operator, description=None):
    if isinstance(operator, AcquisitionModel):
        return operator
    if _my_isscalar(operator):
        return Scalar(operator)
    if isinstance(operator, LinearOperator):
        direct    = lambda input, reusein=False, reuseout=False: \
                        operator.matvec(input)
        transpose = lambda input, reusein=False, reuseout=False: \
                        operator.rmatvec(input)
        model = AcquisitionModelLinear(direct=direct,
                                       transpose=transpose,
                                       shapein=operator.shape[1],
                                       shapeout=operator.shape[0],
                                       dtype=operator.dtype,
                                       description=description)
        return model
    return asacquisitionmodel(scipy.sparse.linalg.aslinearoperator(operator),
                              description=description)


#-------------------------------------------------------------------------------


def _toacquisitionmodel(model, cls, description=None):
    import copy
    if model.__class__ == cls:
        return model
    if description is None:
        description = cls.__name__
    model2 = copy.copy(model)
    model.__class__ = cls
    for key in model.__dict__:
        model.__dict__[key] = None
    model.description = description
    model.blocks = [model2]


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


def _propagate_attributes(input, output):
    """Copy over attributes form input to output"""

    # get common base class
    cls = input.__class__
    while not issubclass(type(output), cls):
        cls = cls.__base__
    # if the arguments do not have the same shape, only copy the units
    if input.shape != output.shape:
        if cls == numpy.ndarray:
            return
        output._unit = input._unit
        output._derived_units = input._derived_units
        return
    # copy over slots
    while cls != numpy.ndarray:            
        for a in input.view(cls).__slots__:
            setattr(output, a, getattr(input, a))
        cls = cls.__base__

#-------------------------------------------------------------------------------


def _smart_reshape(input, shape):
    curr = input
    shape = flatten_sliced_shape(shape)
    while True:
        if curr.shape == shape:
            return curr
        if curr.base is None or curr.base.dtype != input.dtype or \
           curr.base.__class__ != input.__class__:
            return curr.reshape(shape)
        curr = curr.base


#-------------------------------------------------------------------------------


def _validate_input_unit(input, expected):
    if len(expected) == 0 or not hasattr(input, '_unit') or \
       len(input._unit) == 0:
        return
    for u,v in expected.items():
        if u not in input._unit or input._unit[u] != v:
            raise ValidationError("The input unit '" + input.unit + \
                                  "' is incompatible with the required unit '"+\
                                  Quantity(1, expected).unit + "'.")
    return


#-------------------------------------------------------------------------------


def _validate_output_unit(input, output, unitin, unitout):
    if not hasattr(output, '_unit'):
        return
    if len(unitout) == 0:
        return
    if len(output._unit) == 0:
        output._unit = unitout
        return
    output._unit = _divide_unit(output._unit, unitin)
    output._unit = _multiply_unit(output._unit, unitout)
