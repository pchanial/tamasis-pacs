# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
try:
    import fftw3
except:
    print('Warning: Library PyFFTW3 is not installed.')

import copy
import gc
import multiprocessing
import numpy as np
import scipy.signal
import scipy.sparse.linalg
import tamasisfortran as tmf

from mpi4py import MPI
from scipy.sparse.linalg.interface import LinearOperator
from . import var
from .datatypes import Map, Tod, combine_sliced_shape, flatten_sliced_shape, \
                       validate_sliced_shape
from .numpyutils import _my_isscalar
from .processing import interpolate_linear
from .quantity import Quantity, UnitError, _divide_unit, _multiply_unit
from .utils import diff, diffT, diffTdiff, shift
from .mpiutils import split_shape, split_work

__all__ = [
    'AcquisitionModel',
    'AcquisitionModelLinear',
    'DistributionGlobal',
    'DistributionLocal',
    'CircularShift',
    'Clip',
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
    'Maximum',
    'Minimum',
    'Offset',
    'Packing',
    'Padding',
    'Projection',
    'Reshaping',
    'ResponseTruncatedExponential',
    'Rounding',
    'Scalar',
    'Shift',
    'SqrtInvNtt',
    'Unpacking',
    'acquisitionmodel_factory',
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

    The direct method must not rely on input attributes (except Tod's nsamples
    if the acquisition model is unconstrained) since this method is supposed to
    work on bare ndarrays (in which case the acquisition model must be
    constrained: nsamples is extracted from the shapein property). Attribute
    handling must be dealt with in the AcquisitionModel's __init__ method via
    the attrin and attrout keywords.
    """

    def __init__(self, direct=None, cache=False, dtype=None, description=None,
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
        shapein = validate_sliced_shape(shapein)
        self.shapein = shapein
        shapeout = validate_sliced_shape(shapeout or shapein)
        self.shapeout = shapeout
        self.typein = typein
        self.typeout = typeout or typein
        self.unitin = Quantity(1., unitin )._unit
        self.unitout = Quantity(1., unitout or unitin)._unit

        self.cache = cache
        
        if isinstance(self, AcquisitionModelTranspose):
            return

        if cache:
            # store the input of the direct model. Its memory allocation
            # is re-used as the output of the transpose model
            self.cachein = None
            # store the input of the transpose model. Its memory allocation
            # is re-used as the output of the direct model
            self.cacheout = None
        else:
            if typein != (typeout or typein):
                raise ValueError('Inplace handling requires same input and ou' \
                                 'tput type (' + str(typein) + ',' + \
                                 str(typeout or typein) + ').')
            if self.attrin != self.attrout:
                raise ValueError('Inplace handling requires same input and ou' \
                                 'tput attributes.')

        if shapein and type(shapein[-1]) is tuple:
            if not issubclass(typein, Tod):
                raise TypeError('The input type should be a Tod.')
        if shapeout and type(shapeout[-1]) is tuple:
            if not issubclass(typeout or typein, Tod):
                raise TypeError('The output type should be a Tod.')

    def __call__(self, input, inplace=False, cachein=False, cacheout=False):
        return self.direct(input, inplace, cachein, cacheout)

    def direct(self, input, inplace, cachein, cacheout):
        raise NotImplementedError()

    @property
    def shape(self):
        shape = (np.product(flatten_sliced_shape(self.shapeout)),
                 np.product(flatten_sliced_shape(self.shapein)))
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
        self._dtype = np.dtype(dtype)
    
    def validate_shapein(self, shapein):
        """
        Validate input shape and return the output shape of the direct model
        """
        selfshapein = self.shapein
        if shapein is None or shapein == selfshapein:
            return self.shapeout
        if selfshapein is None:
            return shapein
        if flatten_sliced_shape(shapein) == flatten_sliced_shape(selfshapein):
            return self.shapeout
        raise ValidationError('The input of ' + self.description + ' has an i' \
            'ncompatible shape ' + str(shapein) + '. Expected shape is ' + \
            str(self.shapein) + '.')

    def validate_shapeout(self, shapeout):
        """
        Validate input shape and return the output shape of the transpose model
        """
        selfshapeout = self.shapeout
        if shapeout is None or shapeout == selfshapeout:
            return self.shapein
        if selfshapeout is None:
            return shapeout
        if flatten_sliced_shape(shapeout) == flatten_sliced_shape(selfshapeout):
            return self.shapein
        raise ValidationError("The input of '" + self.description + ".T' has " \
            'an incompatible shape ' + str(shapeout) + '. Expected shape is ' +\
            str(self.shapeout) + '.')

    def validate_input_inplace(self, input, inplace):

        input = np.asanyarray(input)
        try:
            input = np.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        except:
            raise TypeError("The input of '" + self.description + \
                "' has a non-numeric type.")

        shape = self.validate_shapein(validate_sliced_shape(input.shape, 
            getattr(input, 'nsamples', None) or (self.shapein[-1] \
            if self.shapein and type(self.shapein[-1]) is tuple else None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + \
                                  self.description+' is not known.')

        typein = self.typein
        if shape and type(shape[-1]) is tuple and (typein is None or \
           not issubclass(typein, Tod)):
            typein = Tod
        if typein is not None and not isinstance(input, typein):
            input = input.view(typein)
        if shape and type(shape[-1]) is tuple:
            input.nsamples = shape[-1]

        if not inplace:
            if var.verbose:
                print('Info: Allocating ' + input.dtype.type.__name__ + \
                      str(input.shape).replace(' ','') +  ' = ' + \
                      str(input.dtype.itemsize * input.size / 2.**20) + \
                      ' MiB in ' + self.description + '.')
            try:
                input = input.copy()
            except MemoryError:
                gc.collect()
                input = input.copy()

        for k,v in self.attrin:
            setattr(input, k, v)

        return input

    def validate_input_direct(self, input, cachein, cacheout):

        def set_cachein(cache):
            self.cachein = cache
        def set_cacheout(cache):
            self.cacheout = cache

        return self.validate_input_cache(input, self.description,
            cachein, cacheout,
            self.attrin, self.attrout,
            self.cachein, self.cacheout,
            set_cachein, set_cacheout,
            self.shapein, self.shapeout,
            self.typein, self.typeout,
            self.unitin, self.unitout,
            lambda shape: self.validate_shapein(shape))

    def validate_input_transpose(self, input, cachein, cacheout):

        def set_cachein(cache):
            self.cachein = cache
        def set_cacheout(cache):
            self.cacheout = cache

        return self.validate_input_cache(input, self.description + '.T',
            cachein, cacheout,
            self.attrout, self.attrin,
            self.cacheout, self.cachein,
            set_cacheout, set_cachein,
            self.shapeout, self.shapein,
            self.typeout, self.typein,
            self.unitout, self.unitin,
            lambda shape: self.validate_shapeout(shape))

    def validate_input_cache(self, input, description,
                             do_cachein, do_cacheout,
                             attrin, attrout,
                             cachein, cacheout,
                             set_cachein, set_cacheout,
                             shapein, shapeout,
                             typein, typeout,
                             unitin, unitout,
                             validate_input_shape):

        input = np.array(input, ndmin=1, subok=True, copy=False)
        try:
            input = np.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        except:
            raise TypeError("The input of '" + description + "' has a non-num" \
                            'eric type.')

        shapein = validate_sliced_shape(input.shape, getattr(input, 'nsamples',\
            None) or (shapein[-1] if shapein and type(shapein[-1]) is tuple \
            else None))
        shapeout = validate_input_shape(shapein)

        if type(shapein[-1]) is tuple and (typein is None or \
           not issubclass(typein, Tod)):
            typein = Tod
        if typein is not None and not isinstance(input, typein):
            input = input.view(typein)
        if shapein and type(shapein[-1]) is tuple:
            input.nsamples = shapein[-1]

        if do_cachein and id(input) != id(cachein):
            # validate input before storing it (nsamples is not enforced)
            _validate_input_unit(input, unitin)
            for k,v in attrin.items():
                setattr(input, k, v)
            
            # store it
            set_cachein(input)

        # get output from the cache
        shapeout_flat = flatten_sliced_shape(shapeout)
        if do_cacheout and cacheout is not None and \
           shapeout_flat == cacheout.shape and cacheout.dtype == input.dtype:
            output = cacheout

        else:

            # allocate output
            if var.verbose:
                reason = 'cache not requested' if not do_cacheout else \
                         'empty cache' if cacheout is None else \
                         'type mismatch' if cacheout.dtype != input.dtype else \
                         'shape mismatch'
                print('Info: Allocating ' + self.dtype.type.__name__ + \
                      str(shapeout_flat).replace(' ','') +  ' = ' + \
                      str(input.dtype.itemsize * np.product(shapeout_flat) \
                      / 2.**20) + ' MiB in ' + description + ' (' + reason + \
                      ').')
            if typeout is None:
                typeout = input.__class__
            if type(shapeout[-1]) is tuple and not issubclass(typeout, Tod):
                typeout = Tod
            if hasattr(typeout, 'empty'):
                try:
                    output = typeout.empty(shapeout, dtype=self.dtype)
                except MemoryError:
                    gc.collect()
                    output = typeout.empty(shapeout, dtype=self.dtype)
            else:
                try:
                    output = np.empty(shapeout_flat, self.dtype)
                except MemoryError:
                    gc.collect()
                    output = np.empty(shapeout_flat, self.dtype)

            # store output
            if do_cacheout:
                set_cacheout(output)

        if type(shapeout[-1]) is tuple:
            output = Tod(output, nsamples=shapeout[-1], copy=False)

        _propagate_attributes(input, output, unitin, unitout, attrout)

        return input, output

    def __mul__(self, other):
        if isinstance(other, np.ndarray):
            return self.matvec(other)
        return Composition([self, other])

    def __rmul__(self, other):
        if not _my_isscalar(other):
            raise NotImplementedError("It is not possible to multiply '" + \
                str(type(other)) + "' with an AcquisitionModel.")
        return Composition([other, self])

    def __imul__(self, other):
        _tocompositemodel(self, Composition, [copy.copy(self), other])
        return self

    def __add__(self, other):
        return Addition([self, other])

    def __radd__(self, other):
        return Addition([other, self])

    def __iadd__(self, other):
        _tocompositemodel(self, Addition, [copy.copy(self), other])
        return self

    def __sub__(self, other):
        return Addition([self, -other])

    def __rsub__(self, other):
        return Addition([other, -self])

    def __isub__(self, other):
        _tocompositemodel(self, Addition, [copy.copy(self), -other])
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
    def __init__(self, transpose=None, **keywords):
        AcquisitionModel.__init__(self, **keywords)
        if transpose is not None:
            if not hasattr(transpose, '__call__'):
                raise TypeError('The input transpose method is not callable.')
            self.transpose = transpose

    def transpose(self, input, inplace, cachein, cacheout):
        raise NotImplementedError()

    def matvec(self, v, inplace=False, cachein=False, cacheout=False):
        v = v.reshape(flatten_sliced_shape(self.shapein))
        return self.direct(v, inplace, cachein, cacheout).ravel()

    def rmatvec(self, v, inplace=False, cachein=False, cacheout=False):
        v = v.reshape(flatten_sliced_shape(self.shapeout))
        return self.transpose(v, inplace, cachein, cacheout).ravel()

    def dense(self):
        d = np.ndarray(self.shape, dtype=self.dtype)
        v = np.zeros(self.shape[1], dtype=var.FLOAT_DTYPE)
        for i in range(self.shape[1]):
            v[:] = 0
            v[i] = 1
            d[:,i] = self.matvec(v, inplace=True, cachein=True, cacheout=True)
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
                                        cache=model.cache,
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
    def cachein(self):
        return self.model.cacheout

    @cachein.setter
    def cachein(self, cache):
        self.model.cacheout = cache
    
    @property
    def cacheout(self):
        return self.model.cachein

    @cacheout.setter
    def cacheout(self, cache):
        self.model.cachein = cache
    
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

    def __init__(self, models):

        models = [ asacquisitionmodel(m) for m in models ]
        self.blocks = []

        for model in models:
            if isinstance(model, self.__class__):
                self.blocks.extend(model.blocks)
            else:
                self.blocks.append(model)

        AcquisitionModel.__init__(self)

        if all([hasattr(m, 'matvec') for m in self.blocks]):
            self.matvec = \
                lambda v, inplace=False, cachein=False, cacheout=False: \
                    self.direct(v.reshape(flatten_sliced_shape(self.shapein)),
                                inplace, cachein, cacheout).ravel()
            def dense():
                d = np.ndarray(self.shape, dtype=self.dtype)
                v = np.zeros(self.shape[1], dtype=var.FLOAT_DTYPE)
                for i in range(self.shape[1]):
                    v[:] = 0
                    v[i] = 1
                    d[:,i] = self.matvec(v, True, True, True)
                return d
            self.dense = dense

        if all([hasattr(m, 'rmatvec') for m in self.blocks]):
            self.rmatvec = \
                lambda v, inplace=False, cachein=False, cacheout=False: \
                    self.transpose(v.reshape(flatten_sliced_shape(
                        self.shapeout)), inplace, cachein, cacheout).ravel()

    @property
    def dtype(self):
        for block in self.blocks:
            if block.dtype.type in (np.complex64, np.complex128, np.complex256):
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
        return np.ndarray

    @typein.setter
    def typein(self, value):
        pass

    @property
    def typeout(self):
        for model in reversed(self.blocks):
            if model.typeout is not None:
                return model.typeout
        return np.ndarray

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

    def validate_input(self, input, shape):
        input = np.array(input, ndmin=1, subok=True, copy=False)
        if shape is not None and type(shape[-1]) is tuple:
            input = Tod(input, nsamples=shape[-1], copy=False)
        return input

    def __str__(self):
        result = AcquisitionModel.__str__(self) + ':'
        components = []
        for block in self.blocks:
            components.extend(str(block).split('\n'))
        result += '\n    '+'\n    '.join(components)
        return result


#-------------------------------------------------------------------------------


class Addition(Composite):
    """
    Class for acquisition models addition

    If at least one of the input already is the result of an addition,
    a flattened list of operators is created by associativity, in order to
    benefit from the AcquisitionModel's caching mechanism.
    """

    def direct(self, input, inplace, cachein, cacheout):
        input = self.validate_input(input, self.shapein)
        output = self.blocks[0].direct(input, False, False, False)
        for i, model in enumerate(self.blocks[1:]):
            last = i == len(self.blocks) - 2
            tmf.add_inplace(np.array(output, ndmin=1, copy=False).T,
                            np.array(model.direct(input, inplace and last,
                            cachein, cacheout), ndmin=1, copy=False).T)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input = self.validate_input(input, self.shapeout)
        output = self.blocks[0].transpose(input, False, False, False)
        for i, model in enumerate(self.blocks[1:]):
            last = i == len(self.blocks) - 2
            tmf.add_inplace(np.array(output, ndmin=1, copy=False).T,
                            np.array(model.transpose(input, inplace and last,
                            cachein, cacheout), ndmin=1, copy=False).T)
        return output

    @property
    def shapein(self):
        shapein = None
        for model in self.blocks:
            shapein_ = model.shapein
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
            shapeout_ = model.shapeout
            if shapeout_ is None:
                continue
            if shapeout is None or type(shapeout_[-1]) is tuple:
                shapeout = shapeout_
                continue
            if flatten_sliced_shape(shapeout) != \
               flatten_sliced_shape(shapeout_):
                raise ValidationError("Incompatible shape in operands: '" + \
                          str(shapeout) +"' and '" + str(shapeout_) + "'.")
        return shapeout

    @shapeout.setter
    def shapeout(self, value):
        pass

    @property
    def T(self):
        return Addition([model.T for model in self.blocks])

    def __iadd__(self, other):
        oldblocks = self.blocks
        if isinstance(other, Addition):
            self.blocks.extend(other.blocks)
        else:
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


#-------------------------------------------------------------------------------


class Composition(Composite):
    """
    Class for acquisition models composition

    If at least one of the input already is the result of a composition,
    a flattened list of operators is created by associativity, in order to
    benefit from the AcquisitionModel's caching mechanism.
    """

    def direct(self, input, inplace, cachein, cacheout):
        input = self.validate_input(input, self.shapein)
        caches = [m.cache for m in self.blocks]
        if any(caches):
            first_cache = caches.index(True)
            last_cache = len(self.blocks) - caches.index(True) - 1
        else:
            first_cache = len(self.blocks)
            last_cache = -1
        for i, model in enumerate(reversed(self.blocks)):
            input = model.direct(input, inplace or i != 0,
                                 cachein or i > first_cache,
                                 cacheout or i < last_cache)
        return input

    def transpose(self, input, inplace, cachein, cacheout):
        input = self.validate_input(input, self.shapeout)
        caches = [m.cache for m in self.blocks]
        first_cache = len(self.blocks) - caches.index(True) - 1
        last_cache = caches.index(True)
        for i, model in enumerate(self.blocks):
            input = model.transpose(input, inplace or i != 0,
                                    cachein or i > first_cache,
                                    cacheout or i < last_cache)
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

    @property
    def T(self):
        return Composition([model.T for model in reversed(self.blocks)])

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


#-------------------------------------------------------------------------------


class Square(object):
    """
    Square operator
    
    The input and output must have the same shape
    This operator does not implement the cache mechanism, but operation on
    the input can be done inplace or on a copy.
    """

    def __init__(self, shapein=None, **keywords):
        self.shapeout = shapein
        self.validate_shapeout = self.validate_shapein
    

#-------------------------------------------------------------------------------


class Symmetric(AcquisitionModelLinear, Square):
    """Symmetric operator"""

    def __init__(self, **keywords):
        AcquisitionModelLinear.__init__(self, **keywords)
        Square.__init__(self, **keywords)
        self.transpose = self.direct
        self.rmatvec = self.matvec

    @property
    def T(self):
        return self


#-------------------------------------------------------------------------------


class Diagonal(Symmetric):
    """
    Diagonal operator.

    Multiply by a diagonal matrix. The input of a Diagonal instance can be of
    rank greater than the specified diagonal array, in which case the latter
    is broadcast along the fast dimensions.
    """

    def __init__(self, diagonal, **keywords):
        diagonal = np.array(diagonal, dtype=var.get_default_dtype(diagonal),
                            order='c')
        Symmetric.__init__(self, dtype=diagonal.dtype, **keywords)
        self.isscalar = diagonal.ndim == 0
        self.data = np.array(diagonal, ndmin=1, copy=False)

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        if self.dtype == input.dtype == var.FLOAT_DTYPE:
            tmf.multiply_inplace(output.T, self.data.T)
        else:
            output.T[:] *= self.data.T
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return self.shapein
        if self.isscalar:
            return shapein
        if flatten_sliced_shape(shapein[0:self.data.ndim]) != self.data.shape:
            raise ValueError('The input has an incompatible shape ' + \
                             str(shapein) + '.')
        return shapein

    def matvec(self, v, inplace=False, cachein=False, cacheout=False):
        shape = list(self.data.shape)
        shape.append(-1)
        v = v.reshape(shape)
        return self.direct(v, inplace, cachein, cacheout).ravel()


class Offset(AcquisitionModel, Square):
    """
    Offset operator.

    Add an offset to the input. The input of an Offset instance can be of rank
    greater than the specified diagonal array, in which case the latter is
    broadcast along the fast dimensions. This operator is not linear.
    """
    def __init__(self, offset, **keywords):
        offset = np.array(offset, dtype=var.get_default_dtype(offset),
                          order='c')
        AcquisitionModel.__init__(self, dtype=offset.dtype, **keywords)
        Square.__init__(self, **keywords)
        self.isscalar = offset.ndim == 0
        self.data = np.array(offset, ndmin=1, copy=False)

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        if self.dtype == var.FLOAT_DTYPE:
            tmf.add_inplace(output.T, self.data.T)
        else:
            output.T[:] += self.data.T
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return self.shapein
        if self.isscalar:
            return shapein
        if flatten_sliced_shape(shapein[0:self.data.ndim]) != self.data.shape:
            raise ValueError('The input has an incompatible shape ' + \
                             str(shapein) + '.')
        return shapein

    def matvec(self, v, inplace=False, cachein=False, cacheout=False):
        shape = list(self.data.shape)
        shape.append(-1)
        v = v.reshape(shape)
        return self.direct(v, inplace, cachein, cacheout).ravel()


#-------------------------------------------------------------------------------


class Rounding(AcquisitionModel, Square):
    """Rounding operator.
    
    The rounding method may be one of the following:
        - rtz : round towards zero (truncation)
        - rti : round towards infinity
        - rtmi : round towards minus infinity (floor)
        - rtpi : round towards positive infinity (ceil)
        - rhtz : round half towards zero
        - rhti : round half towards infinity (numpy's round, fortran's nint)
        - rhtmi : round half towards minus infinity
        - rhtpi : round half towards positive infinity
        - rhs : round half stochastically
    """

    def __init__(self, method='rhti', **keywords):
        AcquisitionModel.__init__(self, **keywords)
        Square.__init__(self, **keywords)
        method = method.lower()
        table = {'rtz'   : tmf.round_rtz,
                 'rti'   : tmf.round_rti,
                 'rtmi'  : tmf.round_rtmi,
                 'rtpi'  : tmf.round_rtpi,
                 'rhtz'  : tmf.round_rhtz,
                 'rhti'  : tmf.round_rhti,
                 'rhtmi' : tmf.round_rhtmi,
                 'rhtpi' : tmf.round_rhtpi,
                 'rhs'   : tmf.round_rhs}
        if method not in table:
            raise ValueError('The rounding method must be one of the following'\
                ': ' + ','.join("'" + k + "'" for k in table.keys()) + '.')
        self.round = table[method]

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        self.round(output.T)
        return output


#-------------------------------------------------------------------------------


class DiscreteDifference(AcquisitionModelLinear, Square):
    """
    Discrete difference operator.

    Calculate the nth order discrete difference along given axis.
    """

    def __init__(self, axis=0, n=1, comm=None, **keywords):
        AcquisitionModelLinear.__init__(self, **keywords)
        Square.__init__(self, **keywords)
        self.n = n
        self.axis = axis
        self.comm = comm or var.comm_map

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        for i in range(self.n):
            diff(output, self.axis, comm=self.comm)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        for i in range(self.n):
            diffT(output, self.axis, comm=self.comm)
        return output        


#-------------------------------------------------------------------------------


class DdTdd(Symmetric):
    """Calculate operator dX.T dX along a given axis."""

    def __init__(self, axis=0, scalar=1., description=None, comm=None,
                 **keywords):
        if description is None and scalar != 1.:
            description = str(scalar) + ' DdTdd'
        Symmetric.__init__(self, **keywords)
        self.axis = axis
        self.scalar = scalar
        self.comm = comm or var.comm_map

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        diffTdiff(output, self.axis, self.scalar, comm=self.comm)
        return output

   
#-------------------------------------------------------------------------------


class Projection(AcquisitionModelLinear):
    """
    This class handles  operations by the pointing matrix

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
        - npixels_per_sample: maximum number of sky map pixels that can be
          intercepted by a detector
    """

    def __init__(self, observation, method=None, header=None, resolution=None,
                 npixels_per_sample=0, oversampling=True, comm_map=None,
                 packed=False, description=None):

        self.comm_map = comm_map or var.comm_map
        self.comm_tod = observation.comm_tod

        self.method, pmatrix, self.header, self.ndetectors, nsamples, \
        self.npixels_per_sample, (unitout, unitin), (duout, duin) = \
            observation.get_pointing_matrix(header,
                                            resolution,
                                            npixels_per_sample,
                                            method=method,
                                            oversampling=oversampling)

        self.nsamples_tot = int(np.sum(nsamples))
        self._pmatrix = pmatrix
        if self.npixels_per_sample == 0:
            pmatrix = np.empty(0, dtype=np.int64)
        self.pmatrix = pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]) \
                              .view(np.recarray)
        self.pmatrix.shape = (self.ndetectors, np.sum(nsamples),
                              self.npixels_per_sample)

        shapein = tuple([self.header['NAXIS'+str(i+1)] for i in \
                         range(self.header['NAXIS'])])[::-1]
        mask = Map.empty(shapein, dtype=np.bool8, header=self.header)
        tmf.pointing_matrix_mask(self._pmatrix, mask.view(np.int8).T, 
            self.npixels_per_sample, self.nsamples_tot, self.ndetectors)

        ismapdistributed = self.comm_map.Get_size() > 1
        istoddistributed = self.comm_tod.Get_size() > 1
        self.ispacked = packed or ismapdistributed
        if self.ispacked:
            tmf.pointing_matrix_pack(self._pmatrix, mask.view(np.int8).T,
                self.npixels_per_sample, self.nsamples_tot, self.ndetectors)
            shapein = (int(np.sum(~mask)))

        attrin = {'header' : self.header}
        if duin is not None:
            attrin['derived_units'] = duin
        attrout = {}
        if duout is not None:
            attrout['derived_units'] = duout
        shapeout = combine_sliced_shape(self.ndetectors, nsamples)
        AcquisitionModelLinear.__init__(self,
                                        cache=True,
                                        description=description,
                                        attrin=attrin,
                                        attrout=attrout,
                                        shapein=shapein,
                                        shapeout=shapeout,
                                        typein=Map,
                                        typeout=Tod,
                                        unitin=unitin,
                                        unitout=unitout)
        self.mask = mask
        if not self.ispacked and not istoddistributed:
            return

        if self.ispacked:
            if ismapdistributed:
                self *= DistributionLocal(self.mask)
            else:
                self *= Packing(self.mask)
        elif istoddistributed:
            self *= DistributionGlobal(self.shapein, share=True,
                                       comm=self.comm_tod)
        s = self.blocks[0]
        self.header = s.header
        self.mask = s.mask
        self.method = s.method
        self.ndetectors = s.ndetectors
        self.npixels_per_sample = s.npixels_per_sample
        self.nsamples_tot = s.nsamples_tot
        self.pmatrix = s.pmatrix

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        tmf.pointing_matrix_direct(self._pmatrix, input.T, output.T,
                                   self.npixels_per_sample)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        tmf.pointing_matrix_transpose(self._pmatrix, input.T, output.T, 
                                      self.npixels_per_sample)
        return output

    def get_ptp(self):
        npixels = np.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample,
            self.nsamples_tot, self.ndetectors, npixels).T


#-------------------------------------------------------------------------------


class DistributionGlobal(AcquisitionModelLinear):
    """
    Distribute a global map to different MPI processes.
    By default, they are locally distributed, in the sense that an MPI process
    will only handle a subset of the global map.
    """

    def __init__(self, shape, share=False, comm=None, **keywords):

        if comm is None:
            comm = var.comm_map
        self.comm = comm

        # if share is true, the maps are not distributed
        if share:
            def direct(input, inplace, cachein, cacheout):
                return self.validate_input_inplace(input, inplace)
            def transpose(input, inplace, cachein, cacheout):
                output = self.validate_input_inplace(input, inplace)
                if self.comm.Get_size() > 1:
                    self.comm.Allreduce(MPI.IN_PLACE, [output, MPI.DOUBLE],
                                        op=MPI.SUM)
                return output
            AcquisitionModelLinear.__init__(self, shapein=shape, typein=Map,
                                            direct=direct, transpose=transpose,
                                            cache=False, **keywords)
            return

        shapeout = split_shape(shape, comm)
        self.counts = []
        self.offsets = [0]
        for rank in range(comm.Get_size()):
            s = split_work(shape[0], rank=rank, comm=comm)
            n = (s.stop - s.start) * np.product(shape[1:])
            self.counts.append(n)
            self.offsets.append(self.offsets[-1] + n)
        self.offsets.pop()
        attrin = { 'comm':MPI.COMM_SELF, 'shape_global':shape }
        attrout = { 'comm':self.comm, 'shape_global':shape }
        AcquisitionModelLinear.__init__(self, cache=True, shapein=shape,
            shapeout=shapeout, attrin=attrin, attrout=attrout, typein=Map,
            **keywords)

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        s = split_work(self.shapein[0], comm=self.comm)
        n = s.stop - s.start
        output[0:n] = input[s.start:s.stop]
        output[n:] = 0
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        s = split_work(self.shapein[0], comm=self.comm)
        n = s.stop - s.start
        self.comm.Allgatherv([input[0:n], MPI.DOUBLE], [output, (self.counts,
                             self.offsets), MPI.DOUBLE])
        return output


#-------------------------------------------------------------------------------


class DistributionLocal(AcquisitionModelLinear):
    """
    Scatter a distributed map to different MPI processes under the control of a
    local non-distributed mask.
    """

    def __init__(self, mask, operator=MPI.SUM, comm=None, **keywords):
        if comm is None:
            comm = var.comm_map
        shapeout = (int(np.sum(~mask)),)
        shapein = split_shape(mask.shape, comm)
        attrin = { 'comm':comm, 'shape_global':mask.shape }
        attrout = { 'comm':MPI.COMM_SELF, 'shape_global':shapeout}
        AcquisitionModelLinear.__init__(self, cache=True, typein=Map,
            shapein=shapein, shapeout=shapeout, attrin=attrin, **keywords)
        self.comm = comm
        self.mask = mask
        self.operator = operator

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        status = tmf.mpi_allscatterlocal(input.T, self.mask.view(np.int8).T,
            output.T, self.comm.py2f())
        if status == 0: return output
        if status < 0:
            raise RuntimeError('Incompatible sizes.')
        raise MPI.Exception(status)

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        status = tmf.mpi_allreducelocal(input.T, self.mask.view(np.int8).T,
            output.T, self.operator.py2f(), self.comm.py2f())
        if status == 0: return output
        if status < 0:
            raise RuntimeError('Incompatible mask.')
        raise MPI.Exception(status)


#-------------------------------------------------------------------------------


class Compression(AcquisitionModelLinear):
    """
    Abstract class for compressing the input signal.
    """

    def __init__(self, compression_factor, shapein=None, **keywords):
        if _my_isscalar(compression_factor):
            compression_factor = (compression_factor,)
        self.factor = np.array(compression_factor, int)

        if np.all(self.factor == 1):
            def direct(input, inplace, cachein, cacheout):
                return self.validate_input_inplace(input, inplace)
            transpose = direct
            cache = False
        else:
            direct = None
            transpose = None
            cache = True
        shapeout = self.validate_shapein(shapein)
        AcquisitionModelLinear.__init__(self,
                                        direct=direct,
                                        transpose=transpose,
                                        cache=cache,
                                        shapein=shapein,
                                        shapeout=shapeout,
                                        typein=Tod,
                                        **keywords)

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        if np.any(np.array(shapein[-1]) % self.factor != 0):
            raise ValidationError('The input timeline size ('+str(shapein[-1])+\
                ') is not an integer times the compression factor (' + \
                str(self.factor)+').')
        return combine_sliced_shape(shapein[0:-1],
                                    np.array(shapein[-1]) / self.factor)

    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return None
        if self.shapeout is not None and flatten_sliced_shape(shapeout) == \
           flatten_sliced_shape(self.shapeout):
            return self.shapein
        return combine_sliced_shape(shapeout[0:-1],
                                    np.array(shapeout[-1]) * self.factor)

    def __str__(self):
        return super(Compression, self).__str__()+' (x'+str(self.factor)+')'
       

#-------------------------------------------------------------------------------


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    """

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        tmf.compression_average_direct(input.T, output.T,
            np.array(input.nsamples, np.int32), self.factor.astype(np.int32))
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        tmf.compression_average_transpose(input.T, output.T,
            np.array(input.nsamples, np.int32), self.factor.astype(np.int32))
        return output


#-------------------------------------------------------------------------------


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number
    specified by the compression factor
    """

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        tmf.downsampling_direct(input.T, output.T,
            np.array(input.nsamples, np.int32), self.factor.astype(np.int32))
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        tmf.downsampling_transpose(input.T, output.T,
            np.array(input.nsamples, np.int32), self.factor.astype(np.int32))
        return output
      

#-------------------------------------------------------------------------------


class Identity(Symmetric):
    """
    Identity class.
    """

    def __init__(self, **keywords):
        Symmetric.__init__(self, **keywords)

    def direct(self, input, inplace, cachein, cacheout):
        return self.validate_input_inplace(input, inplace)
       

#-------------------------------------------------------------------------------


class Scalar(Diagonal):
    """
    Class for scalar multiplication
    """

    def __init__(self, value, **keywords):
        if not np.iscomplex(value):
            value = np.real(value)
        Diagonal.__init__(self, value, **keywords)
       
    def __str__(self):
        return super(self.__class__, self).__str__()+' (' + \
               str(self.data) + ')'
    

#-------------------------------------------------------------------------------


class Masking(Symmetric):
    """
    Mask operator.

    Sets to zero values whose mask is True (non-null). The input of a Masking
    instance can be of rank greater than the speficied mask, in which case the
    latter is broadcast along the fast dimensions.
    """

    def __init__(self, mask, **keywords):
        Symmetric.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        if mask is None:
            print('Warning: input mask is None.')
            mask = False
        mask = np.array(mask, order='c', dtype=np.bool8)
        self.isscalar = mask.ndim == 0
        self.mask = np.array(mask, ndmin=1, copy=False)

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        status = tmf.masking(output.T, self.mask.view(np.int8).T)
        if status != 0: raise RuntimeError()
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return self.shapein
        if self.isscalar:
            return shapein
        if flatten_sliced_shape(shapein[0:self.mask.ndim]) != self.mask.shape:
            raise ValueError('The input has shape ' + str(shapein) + ' incomp' \
                'atible with that of the mask ' + str(self.mask.shape) + '.')
        return shapein


#-------------------------------------------------------------------------------


class Unpacking(AcquisitionModelLinear):
    """
    Convert 1d map into an nd array, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, field=0., **keywords):
        mask = np.array(mask, np.bool8)
        AcquisitionModelLinear.__init__(self,
                                        cache=True,
                                        shapein=np.sum(mask == 0),
                                        shapeout=mask.shape,
                                        **keywords)
        self.mask = mask
        self.field = field

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        tmf.unpack_direct(input.T, self.mask.view(np.int8).T, output.T,
                          self.field)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        tmf.unpack_transpose(input.T, self.mask.view(np.int8).T, output.T)
        return output


#-------------------------------------------------------------------------------


class Packing(AcquisitionModelLinear):
    """
    Convert an nd array in a 1d map, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, field=0., **keywords):
        mask = np.array(mask, np.bool8)
        AcquisitionModelLinear.__init__(self,
                                        cache=True,
                                        shapein=mask.shape,
                                        shapeout=np.sum(mask == 0),
                                        **keywords)
        self.mask = mask
        self.field = field

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        tmf.unpack_transpose(input.T, self.mask.view(np.int8).T, output.T)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        tmf.unpack_direct(input.T, self.mask.view(np.int8).T, output.T,
                          self.field)
        return output


#-------------------------------------------------------------------------------


class Reshaping(AcquisitionModelLinear):
    """
    Reshape arrays
    """

    def __init__(self, shapein, shapeout, **keywords):
        if shapein is None or shapeout is None:
            raise ValueError('The shapes are not defined.')
        if np.product(flatten_sliced_shape(shapein)) != \
           np.product(flatten_sliced_shape(shapeout)):
            raise ValueError('The number of elements of the input and output o'\
                             'f the Reshaping operator are incompatible.')
        AcquisitionModelLinear.__init__(self, shapein=shapein,
                                        shapeout=shapeout, **keywords)

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_direct(input, inplace)
        output = _smart_reshape(output, self.shapeout)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        output = self.validate_input_transpose(input, inplace)
        output = _smart_reshape(output, self.shapein)
        return output

    def validate_input_direct(self, input, inplace):
        input = np.array(input, ndmin=1, copy=not inplace, subok=True)
        shapeout = self.validate_shapein(input.shape)
        return input

    def validate_input_transpose(self, input, inplace):
        input = np.array(input, ndmin=1, copy=not inplace, subok=True)
        shapeout = self.validate_shapeout(input.shape)
        return input


#-------------------------------------------------------------------------------


class ResponseTruncatedExponential(AcquisitionModelLinear, Square):
    """
    ResponseTruncatedExponential(tau)

    Apply a truncated exponential response to the signal

    Parameters
    ==========
    
    tau: number
        Time constant divided by the signal sampling period
    """
    
    def __init__(self, tau, **keywords):
        """
        """
        AcquisitionModelLinear.__init__(self, typein=Tod, **keywords)
        Square.__init__(self, **keywords)
        if hasattr(tau, 'SI'):
            tau = tau.SI
            if tau.unit != '':
                raise ValueError('The time constant must be dimensionless.')
        self.tau = np.array(tau, dtype=var.FLOAT_DTYPE, ndmin=1)

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        tmf.convolution_trexp_direct(output.T, np.array(output.nsamples),
                                     self.tau)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        tmf.convolution_trexp_transpose(output.T, np.array(output.nsamples),
                                        self.tau)
        return output


#-------------------------------------------------------------------------------


class Padding(AcquisitionModelLinear):
    "Pads before and after a Tod."

    def __init__(self, left=0, right=0, value=0., shapein=None, **keywords):
        if shapein is not None:
            shapeout = self.validate_shapein(shapein)
        else:
            shapeout = None
        AcquisitionModelLinear.__init__(self,
                                        cache=True,
                                        shapein=shapein,
                                        shapeout=shapeout,
                                        typein=Tod,
                                        **keywords)
        left = np.array(left, ndmin=1, dtype=int)
        right = np.array(right, ndmin=1, dtype=int)
        if np.any(left < 0):
            raise ValueError('Left padding is not positive.')
        if np.any(right < 0):
            raise ValueError('Right padding is not positive.')
        if np.rank(left) != 1 or np.rank(right) != 1:
            raise ValueError('Padding must be scalar or a vector.')
        self.left  = tuple(left)
        self.right = tuple(right)
        self.value = value
   
    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        dest = 0
        dest_padded = 0
        for islice in range(len(input.nsamples)):
            nsamples = input.nsamples[islice]
            left = self.left[islice if len(self.left) > 1 else 0]
            output[...,dest_padded:dest_padded+left] = self.value
            output[...,dest_padded+left:dest_padded+left+nsamples] = \
                input[...,dest:dest+nsamples]
            output[...,dest_padded+left+nsamples:dest_padded+ \
                output.nsamples[islice]] = self.value
            dest += nsamples
            dest_padded += output.nsamples[islice]
        return output
   
    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
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

    def validate_input_direct(self, input, cachein, cacheout):
        input, output = super(Padding, self).validate_input_direct(input,
            cachein, cacheout)
        if len(self.left) != 1 and len(self.left) != len(input.nsamples):
            raise ValueError("The input Tod has a number of slices '" + \
                             str(len(input.nsamples)) + \
                             "' incompatible with the specified padding.")
        return input, output
       
    def validate_input_transpose(self, input, cachein, cacheout):
        input, output = super(Padding, self).validate_input_transpose(input,
            cachein, cacheout)
        if len(self.left) != 1 and len(self.left) != len(input.nsamples):
            raise ValueError("The input Tod has a number of slices '" + \
                             str(len(input.nsamples)) +
                             "' incompatible with the specified padding.")
        return input, output

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        return combine_sliced_shape(shapein[0:-1], np.array(shapein[-1]) + \
                                    self.left + self.right)
       
    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return None
        return combine_sliced_shape(shapeout[0:-1], np.array(shapeout[-1]) -\
                                    self.left - self.right)


#-------------------------------------------------------------------------------


class Shift(AcquisitionModelLinear, Square):

    def __init__(self, n, axis=None, **keywords):
        AcquisitionModelLinear.__init__(self, **keywords)
        Square.__init__(self, **keywords)
        if axis is None:
            if not isinstance(n, (list, tuple, np.ndarray)):
                n = (n,)
            axis = tuple(np.arange(-len(n), 0))
        elif not isinstance(axis, (list, tuple, np.ndarray)):
            n = (n,)
            axis = (axis,)
        elif not isinstance(n, (list, tuple, np.ndarray)) or \
             len(n) !=  len(axis):
            n = len(axis) * (n,)
        self.n = [np.array(v, dtype=int) for v in n]
        self.axis = [int(v) for v in axis]

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        for n, axis in zip(self.n, self.axis):
            shift(output, n, axis)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        for n, axis in zip(self.n, self.axis):
            shift(output, -n, axis)
        return output        


#-------------------------------------------------------------------------------


class CircularShift(AcquisitionModelLinear, Square):
    
    def __init__(self, n, axis=None, **keywords):
        AcquisitionModelLinear.__init__(self, **keywords)
        Square.__init__(self, **keywords)
        if _my_isscalar(n):
            n = (n,)
        if axis is None:
            axis = tuple(np.arange(-len(n), 0))
        elif _my_isscalar(axis):
            axis = (axis,)
        self.n = tuple(map(int, n))
        self.axis = tuple(map(int, axis))

    def direct(self, input, inplace, cachein, cacheout):
        for axis, n in zip(self.axis, self.n):
            input = np.roll(input, -n, axis=axis)
        return input

    def transpose(self, input, inplace, cachein, cacheout):
        for axis, n in zip(self.axis, self.n):
            input = np.roll(input, n, axis=axis)
        return input


#-------------------------------------------------------------------------------


class Fft(AcquisitionModelLinear, Square):
    """
    Performs complex fft
    """

    def __init__(self, shape, flags=['estimate'], **keywords):
        AcquisitionModelLinear.__init__(self, shapein=shape,
            dtype=var.COMPLEX_DTYPE, **keywords)
        Square.__init__(self, **keywords)
        if fftw3.planning.lib_threads is None:
            nthreads = 1
        else:
            nthreads = tmf.info_nthreads()
        self.n = np.product(shape)
        self._in  = np.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self._out = np.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self.forward_plan = fftw3.Plan(self._in, self._out, direction='forward',
                                       flags=flags, nthreads=nthreads)
        self.backward_plan= fftw3.Plan(self._in, self._out,direction='backward',
                                       flags=flags, nthreads=nthreads)

    def direct(self, input, inplace, cachein, cacheout):
        self._in[:] = input
        fftw3.execute(self.forward_plan)
        return Map(self._out)

    def transpose(self, input, inplace, cachein, cacheout):
        self._in[:] = input
        fftw3.execute(self.backward_plan)
        return Map(self._out / self.n, copy=False)


#-------------------------------------------------------------------------------


class FftHalfComplex(AcquisitionModelLinear, Square):
    """
    Performs real-to-half-complex fft
    """

    def __init__(self, nsamples, **keywords):
        AcquisitionModelLinear.__init__(self, typein=Tod, **keywords)
        Square.__init__(self, **keywords)
        self.nsamples = tuple(np.array(nsamples, ndmin=1, dtype=int))
        self.forward_plan = np.empty(len(self.nsamples), dtype=int)
        self.backward_plan = np.empty(len(self.nsamples), dtype=int)
        for i, n in enumerate(self.nsamples):
            tarray = np.empty(n, dtype=var.FLOAT_DTYPE)
            farray = np.empty(n, dtype=var.FLOAT_DTYPE)
            self.forward_plan[i] = \
                fftw3.Plan(tarray, farray, direction='forward',
                           flags=['measure'], realtypes=['halfcomplex r2c'],
                           nthreads=1)._get_parameter()
            self.backward_plan[i] = \
                fftw3.Plan(farray, tarray, direction='backward',
                           flags=['measure'], realtypes=['halfcomplex c2r'],
                           nthreads=1)._get_parameter()

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        output_ = _smart_reshape(output, (np.product(input.shape[:-1]),
                                 input.shape[-1]))
        tmf.fft_plan(output_.T, np.array(self.nsamples), self.forward_plan)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        output_ = _smart_reshape(output, (np.product(input.shape[:-1]), 
                                 input.shape[-1]))
        tmf.fft_plan(output_.T, np.array(self.nsamples), self.backward_plan)
        dest = 0
        for n in self.nsamples:
            output_[:,dest:dest+n] /= n
            dest += n
        return output

    def validate_shapein(self, shape):
        if shape is None:
            return None
        nsamples = shape[-1]
        if nsamples != self.nsamples and nsamples != sum(self.nsamples):
            raise ValidationError("Invalid FFT size '" + str(nsamples) + \
                                  "' instead of '"+str(self.nsamples)+"'.")
        return combine_sliced_shape(shape[0:-1], self.nsamples)


#-------------------------------------------------------------------------------


class Convolution(Symmetric):

    def __init__(self, kernel, **keywords):
        Symmetric.__init__(self, **keywords)
        self.kernel = np.asanyarray(kernel)

    def direct(self, input, inplace, cachein, cacheout):
        output = self.validate_input_inplace(input, inplace)
        output[:] = scipy.signal.fftconvolve(input, self.kernel, mode='same')
        return output


#-------------------------------------------------------------------------------


class InvNtt(AcquisitionModelLinear):

    def __init__(self, obs, filename=None, **keywords):
        nsamples = obs.get_nsamples()
        length = np.asarray(2**np.ceil(np.log2(np.array(nsamples) + 200)), int)
        invntt = self._get_diagonal(length, obs.get_filter_uncorrelated(
                                    filename=filename, **keywords))
        fft = FftHalfComplex(length)
        padding = Padding(left=invntt.ncorrelations, right=length - nsamples - \
                          invntt.ncorrelations)
        _tocompositemodel(self, Composition,
                          [ padding.T, fft.T, invntt, fft, padding ])

    def _get_diagonal(self, nsamples, filter, **keywords):
        nsamples = np.asarray(nsamples)
        ndetectors = filter.shape[-2]
        ncorrelations = filter.shape[-1] - 1
        nslices = nsamples.size
        if np.rank(filter) == 2:
            filter = np.resize(filter,(nslices, ndetectors, ncorrelations+1))
        tod_filter, status = \
            tmf.fft_filter_uncorrelated(filter.T, np.asarray(nsamples, 
                                        dtype=np.int32), np.sum(nsamples))
        if status != 0: raise RuntimeError()
        d = Diagonal(tod_filter.T, shapein=tod_filter.T.shape, **keywords)
        d = np.maximum(d, 0)
        d.data /= var.comm_tod.allreduce(np.max(d.data), op=MPI.MAX)
        d.ncorrelations = ncorrelations
        return d


#-------------------------------------------------------------------------------


class SqrtInvNtt(InvNtt):
    def __init__(self, *args, **kw):
        invntt = InvNtt(*args, **kw)
        _tocompositemodel(self, Composition, invntt.blocks)
        data = self.blocks[2].data
        data[:] = np.sqrt(data)
        #np.sqrt(data, out=data) does not work with numpy 1.5


#-------------------------------------------------------------------------------


class InterpolationLinear(AcquisitionModelLinear, Square):

    def __init__(self, mask, **keywords):
        AcquisitionModelLinear.__init__(self, attrin={'mask':mask}, **keywords)
        Square.__init__(self, **keywords)
        self.mask = mask

    def direct(self, input, inplace, cachein, cacheout):
        return interpolate_linear(output)

    def transpose(self, input, inplace, cachein, cacheout):
        raise NotImplementedError()


#-------------------------------------------------------------------------------


def acquisitionmodel_factory(direct, transpose=None, description=None, **keywords):
    """Creates an AcquisitionModel from a function"""
    description = description or direct.__name__
    if description == '<lambda>':
        description = ''

    if transpose is None:
        a = AcquisitionModel(description=description, **keywords)
    else:
        a = AcquisitionModelLinear(description=description, **keywords)

    def d(input, inplace, cachein, cacheout):
        output = a.validate_input_inplace(input, inplace)
        return direct(output)
    a.direct = d
    if transpose is None:
        return a

    def t(input, inplace, cachein, cacheout):
        output = a.validate_input_inplace(input, inplace)
        return transpose(output)
    a.transpose = t
    return a

def Clip(vmin, vmax, description=None, **keywords):
    description = description or 'Clip'
    return acquisitionmodel_factory(lambda x: np.clip(x, vmin, vmax, out=x),
                                    description=description, **keywords)

def Maximum(value, description=None, **keywords):
    description = description or 'Maximum'
    return acquisitionmodel_factory(lambda x: np.maximum(x, value, x),
                                    description=description, **keywords)

def Minimum(value, description=None, **keywords):
    description = description or 'Minimum'
    return acquisitionmodel_factory(lambda x: np.minimum(x, value, x),
                                    description=description, **keywords)


#-------------------------------------------------------------------------------


def asacquisitionmodel(operator, shapein=None, shapeout=None, description=None):
    if isinstance(operator, AcquisitionModel):
        if shapein and operator.shapein and shapein != operator.shapein:
            raise ValueError('The input shapein ' + str(shapein) + ' is incom' \
                'patible with that of the input ' + str(operator.shapein) + '.')
        if shapeout and operator.shapeout and shapeout != operator.shapeout:
            raise ValueError('The input shapeout ' + str(shapeout) + ' is inco'\
                'mpatible with that of the input ' + str(operator.shapeout) +  \
                '.')
        if shapein and not operator.shapein or \
           shapeout and not operator.shapeout:
            operator = copy.copy(operator)
            operator.shapein = shapein
            operator.shapeout = shapeout
        return operator
    if _my_isscalar(operator):
        return Scalar(operator)
    if isinstance(operator, LinearOperator):
        direct    = lambda input, inplace, cachein, cacheout: \
                        operator.matvec(input)
        transpose = lambda input, inplace, cachein, cacheout: \
                        operator.rmatvec(input)
        model = AcquisitionModelLinear(direct=direct,
                                       transpose=transpose,
                                       shapein=shapein or operator.shape[1],
                                       shapeout=shapeout or operator.shape[0],
                                       dtype=operator.dtype,
                                       description=description)
        return model
    return asacquisitionmodel(scipy.sparse.linalg.aslinearoperator(operator),
                              description=description)


#-------------------------------------------------------------------------------


def _tocompositemodel(model, cls, models):
    if model.__class__ == cls:
        return model
    model.__class__ = cls
    model.__init__(models)


#-------------------------------------------------------------------------------


def _get_dtype(type1, type2):
    t1 = type1.type()
    t2 = type2.type()
    t = t1 * t2
    return t.dtype


#-------------------------------------------------------------------------------


def _is_scientific_dtype(dtype):
    """Return true if the data type is """
    return issubclass(dtype.type, np.number) or dtype.type == np.bool8


#-------------------------------------------------------------------------------


def _propagate_attributes(input, output, unitin, unitout, attrout):
    """Copy over attributes from input to output"""

    # if the arguments do not have the same shape, only copy the units
    if input.shape != output.shape:
        try:
            setattr(output, '_unit', input._unit)
        except:
            pass
        try:
            setattr(output, '_derived_units', input._derived_units)
        except:
            pass

    elif hasattr(input, '__dict__'):

        # copy over input's attributes
        for k, v in input.__dict__.items():
            setattr(output, k, v)

    _validate_output_unit(input, output, unitin, unitout)

    # copy over operator's attributes
    for k,v in attrout.items():
        setattr(output, k, v)


#-------------------------------------------------------------------------------


def _smart_reshape(input, shape):
    curr = input
    shape = flatten_sliced_shape(shape)
    while True:
        if curr.shape == shape:
            return curr
        base = curr.base
        if base is None or base.dtype != input.dtype or \
           base.__class__ != input.__class__ or base.size != input.size:
            return curr.reshape(shape)
        curr = base


#-------------------------------------------------------------------------------


def _validate_input_unit(input, expected):
    if len(expected) == 0 or not hasattr(input, '_unit') or \
       len(input._unit) == 0:
        return
    for u,v in expected.items():
        if u not in input._unit or input._unit[u] != v:
            unit = Quantity(1, expected).unit
            raise ValidationError("The input unit '" + input.unit + "' is inco"\
                "mpatible with the required unit '" + unit + "'.")


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
