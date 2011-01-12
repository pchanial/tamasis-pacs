try:
    import fftw3
except:
    print('Warning: Library PyFFTW3 is not installed.')

import multiprocessing
import numpy
import scipy.sparse.linalg
import tamasisfortran as tmf

from mpi4py import MPI
from . import var
from .datatypes import Map, Tod, combine_sliced_shape, flatten_sliced_shape, validate_sliced_shape
from .numpyutils import _my_isscalar
from .processing import interpolate_linear
from .quantity import Quantity, UnitError, _divide_unit, _multiply_unit

__all__ = ['AcquisitionModel', 'Diagonal', 'Scalar', 'Identity', 'DiscreteDifference',
           'Projection', 'CompressionAverage', 'DownSampling', 'Masking',
           'Unpacking', 'Reshaping', 'Padding', 'Fft', 'FftHalfComplex',
           'InvNtt', 'InterpolationLinear', 'AllReduce', 'CircularShift', 
           'asacquisitionmodel']


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

    def __init__(self, description=None, types=None, units=None, derived_units=None):

        if description is None:
            description = self.__class__.__name__

        self.description = description

        # input and output shape are unconstrained
        self._shapein  = None
        self._shapeout = None
        
        # store type information
        if type(types) not in (list,tuple):
            types = (types,types)
        self.types = types

        # store unit information
        if type(units) not in (list,tuple):
            units = (units,units)
        self.units = (Quantity(1., units[0])._unit, Quantity(1., units[1])._unit)

        # store derived units information
        if type(derived_units) not in (list,tuple):
            self.derived_units = (derived_units, derived_units)
        else:
            self.derived_units = tuple(derived_units)

        # store the input of the transpose model. Its memory allocation is re-used as the output of the direct model
        self._output_direct = None

        # store the input of the direct model. Its memory allocation is re-used as the output of the transpose model
        self._output_transpose = None

        # dtype
        self._dtype = None

        # Acquisition model components
        self.blocks = []
    
    def __call__(self, input, reusein=False, reuseout=False):
        return self.direct(input, reusein, reuseout)

    def direct(self, input, reusein=False, reuseout=False):
        raise NotImplementedError()

    def transpose(self, input, reusein=False, reuseout=False):
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
        return var.FLOAT_DTYPE

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
    
    def validate_shape_direct(self, shapein):
        """
        Validate input shape and return the output shape of the direct model
        """
        if shapein is None or self.shapein is None:
            if self.shapeout is not None:
               return self.shapeout
            else:
               return shapein
        if flatten_sliced_shape(shapein) != flatten_sliced_shape(self.shapein):
            raise ValidationError('The input of '+self.__class__.__name__+' has an incompatible shape '+str(shapein)+' instead of '+str(self.shapein)+'.')
        return self.shapeout

    def validate_shape_transpose(self, shapeout):
        """
        Validate input shape and return the output shape of the transpose model
        """
        if shapeout is None or self.shapeout is None:
            if self.shapein is not None:
                return self.shapein
            else:
                return shapeout
        if flatten_sliced_shape(shapeout) != flatten_sliced_shape(self.shapeout):
            raise ValidationError("The input of '" + type(self).__name__ + ".T' has an incompatible shape " + str(shapeout) + ' instead of ' + str(self.shapeout) + '.')
        return self.shapein

    def validate_input_direct(self, input, reusein, reuseout, **options):

        # validate input
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + "' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        if self.types[1] is not None and not isinstance(input, self.types[1]):
            input = input.view(self.types[1])
        _validate_input_unit(input, self.units[1])

        # store input for the transpose's output
        if reusein:
            self._output_transpose = input

        # validate output
        shape = self.validate_shape_direct(validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + type(self).__name__+' is not known.')
        shape_flat = flatten_sliced_shape(shape)
        typeout = self.types[0] or input.__class__
        if self._output_direct is None or shape_flat != self._output_direct.shape or \
           self._output_direct.dtype != input.dtype or typeout != type(self._output_transpose):
            self._output_direct = None
            if var.verbose: print('Info: Allocating '+str(input.dtype.itemsize*numpy.product(shape_flat)/2.**20)+' MiB for the output of ' + type(self).__name__+'.')
            if typeout == numpy.ndarray:
                output = numpy.empty(shape, input.dtype, **options)
            else:
                output = typeout.empty(shape, dtype=input.dtype, **options)
            if reuseout:
                self._output_direct = output
        else:
            output = self._output_direct
            if not reuseout:
                self._output_direct = None
        _propagate_attributes(input, output)
        _validate_output_unit(input, output, self.units[1], self.units[0], self.derived_units[0])

        return input, output
        

    def validate_input_transpose(self, input, reusein, reuseout, **options):

        # validate input
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + ".T' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        if self.types[0] is not None and not isinstance(input, self.types[0]):
            input = input.view(self.types[0])
        _validate_input_unit(input, self.units[0])

        # store input for the direct's output
        if reusein:
            self._output_direct = input

        # validate output
        shape = self.validate_shape_transpose(validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + type(self).__name__+' is not known.')
        shape_flat = flatten_sliced_shape(shape)
        typeout = self.types[1] or input.__class__
        if self._output_transpose is None or shape_flat != self._output_transpose.shape or \
           self._output_transpose.dtype != input.dtype or typeout != type(self._output_transpose):
            self._output_transpose = None
            if var.verbose: print('Info: Allocating '+str(input.dtype.itemsize*numpy.product(shape_flat)/2.**20)+' MiB for the output of ' + type(self).__name__+'.T.')
            if typeout == numpy.ndarray:
                output = numpy.empty(shape_flat, input.dtype, **options)
            else:
                output = typeout.empty(shape_flat, dtype=input.dtype, **options)
            if reuseout:
                self._output_transpose = output
        else:
            output = self._output_transpose
            if not reuseout:
                self._output_transpose = None
        _propagate_attributes(input, output)
        _validate_output_unit(input, output, self.units[0], self.units[1], self.derived_units[1])

        return input, output

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

    def direct(self, input, reusein=False, reuseout=False):
        return self.model.transpose(input, reusein=reusein, reuseout=reuseout)

    def transpose(self, input, reusein=False, reuseout=False):
        return self.model.direct(input, reusein=reusein, reuseout=reuseout)

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

    def validate_shape_direct(self, shapein):
        return self.model.validate_shape_transpose(shapein)

    def validate_shape_transpose(self, shapeout):
        return self.model.validate_shape_direct(shapeout)

    def validate_input_direct(self, input, reusein, reuseout):
        raise NotImplementedError()

    def validate_input_transpose(self, input, reusein, reuseout):
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
   
    def direct(self, input, reusein=False, reuseout=False):
        output = self.blocks[0].direct(input, reusein=False, reuseout=False)
        for i, model in enumerate(self.blocks[1:]):
            reusein_ = reusein and i == len(self.blocks)-1
            output += model.direct(input, reusein=reusein_, reuseout=True)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.blocks[0].transpose(input, reusein=False, reuseout=False)
        for i, model in enumerate(self.blocks[1:]):
            reusein_ = reusein and i == len(self.blocks)-1
            output += model.transpose(input, reusein=reusein_, reuseout=True)
        return output

    @property
    def T(self):
        return Addition([model.T for model in self.blocks])

    @property
    def shapein(self):
        shapein = None
        for model in reversed(self.blocks):
            shapein_ = model.validate_shape_transpose(None)
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
            shapeout_ = model.validate_shape_direct(None)
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

    def direct(self, input, reusein=False, reuseout=False):
        for i, model in enumerate(self.blocks):
            reusein_  = reusein  or i != 0
            reuseout_ = reuseout or i != len(self.blocks)-1
            input = model.direct(input, reusein=reusein_, reuseout=reuseout_)
        return input

    def transpose(self, input, reusein=False, reuseout=False):
        for i, model in enumerate(reversed(self.blocks)):
            reusein_  = reusein  or i != 0
            reuseout_ = reuseout or i != len(self.blocks)-1
            input = model.transpose(input, reusein=reusein_, reuseout=reuseout_)
        return input

    @property
    def T(self):
        return Composition([model.T for model in self.blocks])

    @property
    def shapein(self):
        shapein = None
        for model in reversed(self.blocks):
            shapein = model.validate_shape_transpose(shapein)
        return shapein

    @property
    def shapeout(self):
        shapeout = None
        for model in self.blocks:
            shapeout = model.validate_shape_direct(shapeout)
        return shapeout

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


class SquareBase(AcquisitionModel):
    """
    Square operator
    
    The input and output must have the same number of elements.
    This operator does not implement the cache mechanism, but operation on the input
    can be done inplace or on a copy.
    """

    def validate_input_direct(self, input, inplace):
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + "' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        if self.types[1] is not None and not isinstance(input, self.types[1]):
            input = input.view(self.types[1])

        shape = self.validate_shape_direct(validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + type(self).__name__+' is not known.')

        if inplace:
            return input

        if var.verbose:
            print('Info: Allocating ' + str(input.dtype.itemsize * numpy.product(input.shape)/2.**20) + ' MiB in ' + \
                  type(self).__name__ + '.')
        return input.copy('a')

    def validate_input_transpose(self, input, inplace):
        input = numpy.asanyarray(input)
        if not _is_scientific_dtype(input.dtype):
            raise TypeError("The input of '" + type(self).__name__ + ".T' has a non-numeric type '" + str(input.dtype) + "'.")
        input = numpy.asanyarray(input, _get_dtype(self.dtype, input.dtype))
        if self.types[0] is not None and not isinstance(input, self.types[0]):
            input = input.view(self.types[0])

        shape = self.validate_shape_transpose(validate_sliced_shape(input.shape, getattr(input, 'nsamples', None)))
        if shape is None:
            raise ValidationError('The shape of the output of ' + type(self).__name__+' is not known.')

        if inplace:
            return input

        if var.verbose:
            print('Info: Allocating ' + str(input.dtype.itemsize * numpy.product(input.shape)/2.**20) + ' MiB in ' + \
                  type(self).__name__ + '.T.')
        return input.copy('a')


#-------------------------------------------------------------------------------


class Square(SquareBase):
    """
    Square operator
    
    The input and output must have the same shape
    This operator does not implement the cache mechanism, but operation on the input
    can be done inplace or on a copy.
    """

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

    def validate_shape(self, shape):
        return super(Square, self).validate_shape_direct(shape)

    def validate_shape_direct(self, shape):
        return self.validate_shape(shape)

    def validate_shape_transpose(self, shape):
        return self.validate_shape(shape)

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
    def __init__(self, diagonal, nsamples=None, description=None):
        AcquisitionModel.__init__(self, description)
        self.diagonal = numpy.asarray(diagonal, dtype=var.get_default_dtype(diagonal))
        if not _my_isscalar(self.diagonal):
            self.shapein = validate_sliced_shape(self.diagonal.shape, nsamples)

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
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

    def direct(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_direct(input, reusein, reuseout)
        output[:] = numpy.diff(input, self.n, self.axis)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        shape = list(input.shape)
        shape[self.axis] = 1
        border = numpy.zeros(shape)
        out = numpy.concatenate((border, input, border), axis=self.axis)
        output[:] = -numpy.diff(out, axis=self.axis)
        return output

    def validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        shapeout = list(shapein)
        shapeout[self.axis] -= self.n
        return tuple(shapeout)

    def validate_shape_transpose(self, shapeout):
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
    def __init__(self, observation, method=None, header=None, resolution=None, npixels_per_sample=0, oversampling=True, description=None):
        self._pmatrix, self.header, ndetectors, nsamples, self.npixels_per_sample, units, derived_units = \
            observation.get_pointing_matrix(header, resolution, npixels_per_sample, method=method, oversampling=oversampling)
        AcquisitionModel.__init__(self, description, types=(Tod,Map), units=units, derived_units=derived_units)

        self.pmatrix = self._pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]).view(numpy.recarray)
        self.pmatrix.resize((ndetectors, numpy.sum(nsamples), self.npixels_per_sample))
        self.shapein = tuple([self.header['naxis'+str(i+1)] for i in reversed(list(range(self.header['naxis'])))])
        self.shapeout = combine_sliced_shape(ndetectors, nsamples)

    def direct(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_direct(input, reusein, reuseout)
        output.nsamples = self.shapeout[-1]
        output.derived_units = self.derived_units[0]
        tmf.pointing_matrix_direct(self._pmatrix, input.T, output.T, self.npixels_per_sample)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        output.header = self.header
        output.derived_units = self.derived_units[1]
        tmf.pointing_matrix_transpose(self._pmatrix, input.T, output.T,  self.npixels_per_sample)
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
        AcquisitionModel.__init__(self, description, types=Tod)
        if hasattr(compression_factor, 'compression_factor'):
            compression_factor = compression_factor.compression_factor
        elif hasattr(compression_factor, 'slice'):
            compression_factor = compression_factor.slice.compression_factor
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
        output.nsamples = tuple(numpy.divide(input.nsamples, self.factor))
        self.compression_direct(input, output, self.factor)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        if numpy.all(numpy.array(self.factor) == 1):
            if not reusein: return input.copy('a')
            return input
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        output.nsamples = tuple(numpy.multiply(input.nsamples, self.factor))
        self.compression_transpose(input, output, self.factor)
        return output

    def validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        if numpy.any(numpy.array(shapein[-1]) % self.factor != 0):
            raise ValidationError('The input timeline size ('+str(shapein[-1])+') is not an integer times the compression factor ('+str(self.factor)+').')
        return combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) / self.factor)

    def validate_shape_transpose(self, shapeout):
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

    def compression_direct(self, input, output, compression_factor):
        tmf.compression_average_direct(input.T, output.T, numpy.resize(compression_factor, [len(input.nsamples)]))
        output.mask = None
        
    def compression_transpose(self, input, output, compression_factor):
        tmf.compression_average_transpose(input.T, output.T, numpy.resize(compression_factor, [len(input.nsamples)]))
        output.mask = None
        
        

#-------------------------------------------------------------------------------


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number specified by the compression factor
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, compression_factor, description)

    def compression_direct(self, input, output, compression_factor):
        tmf.downsampling_direct(input.T, output.T, numpy.resize(compression_factor, [len(input.nsamples)]))
        output.mask = None
        
    def compression_transpose(self, input, output, compression_factor):
        tmf.downsampling_transpose(input.T, output.T, numpy.resize(compression_factor, [len(input.nsamples)]))
        output.mask = None
      

#-------------------------------------------------------------------------------


class Identity(Symmetric):
    """
    Identity class.
    """
    def __init__(self, description=None):
        AcquisitionModel.__init__(self, description)

    def direct(self, input, reusein=False, reuseout=False):
        return self.validate_input(input, reusein and reuseout)
       

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
        return var.FLOAT_DTYPE


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

    def direct(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_direct(input, reusein, reuseout)
        tmf.unpack_direct(input.T, self.mask.view(numpy.int8).T, output.T, self.field)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        tmf.unpack_transpose(input.T, self.mask.view(numpy.int8).T, output.T)
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

    def direct(self, input, reusein=False, reuseout=False):
        output = numpy.zeros(self.shapeout) + self.field
        output[self.mask] = input
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        return input[self.mask]


#-------------------------------------------------------------------------------


class Reshaping(SquareBase):
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

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input_direct(input, reusein and reuseout)
        output = _smart_reshape(output, self.shapeout)
        if type(self.shapeout[-1]) is tuple:
            output = Tod(output, nsamples=self.shapeout[-1], copy=False)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input_transpose(input, reusein and reuseout)
        output = _smart_reshape(output, self.shapein)
        if type(self.shapein[-1]) is tuple:
            output = Tod(output, nsamples=self.shapein[-1], copy=False)
        return output


#-------------------------------------------------------------------------------


class Padding(AcquisitionModel):
    "Pads before and after a Tod."
    def __init__(self, left=0, right=0, value=0., description=None):
        AcquisitionModel.__init__(self, description, types=Tod)
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
        output.nsamples = tuple(numpy.array(input.nsamples) + self.left + self.right)
        dest = 0
        dest_padded = 0
        for islice in range(len(input.nsamples)):
            nsamples = input.nsamples[islice]
            left = self.left[islice if len(self.left) > 1 else 0]
            output[...,dest_padded:dest_padded+left] = self.value
            output[...,dest_padded+left:dest_padded+left+nsamples] = input[...,dest:dest+nsamples]
            output[...,dest_padded+left+nsamples:] = self.value
            dest += nsamples
            dest_padded += output.nsamples[islice]
        return output
   
    def transpose(self, input, reusein=False, reuseout=False):
        input, output = self.validate_input_transpose(input, reusein, reuseout)
        output.nsamples = tuple(numpy.array(input.nsamples) - self.left - self.right)
        dest = 0
        dest_padded = 0
        for islice in range(len(input.nsamples)):
            nsamples = output.nsamples[islice]
            left  = self.left [islice if len(self.left)  > 1 else 0]
            output[...,dest:dest+nsamples] = input[...,dest_padded+left:dest_padded+left+nsamples]
            dest += nsamples
            dest_padded += input.nsamples[islice]
        return output

    def validate_input_direct(self, input, reusein, reuseout):
        input, output = super(Padding, self).validate_input_direct(input, reusein, reuseout)
        if len(self.left) != 1 and len(self.left) != len(input.nsamples):
            raise ValueError("The input Tod has a number of slices '" + str(len(input.nsamples)) +
                             "' incompatible with the specified padding.")
        return input, output
       
    def validate_input_transpose(self, input, reusein, reuseout):
        input, output = super(Padding, self).validate_input_transpose(input, reusein, reuseout)
        if len(self.left) != 1 and len(self.left) != len(input.nsamples):
            raise ValueError("The input Tod has a number of slices '" + str(len(input.nsamples)) +
                             "' incompatible with the specified padding.")
        return input, output

    def validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        return combine_sliced_shape(shapein[0:-1], numpy.array(shapein[-1]) + self.left + self.right)
       
    def validate_shape_transpose(self, shapeout):
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
        AcquisitionModel.__init__(self, description)
        if fftw3.planning.lib_threads is None:
            nthreads = 1
        else:
            nthreads = tmf.info_nthreads()
        self.shapein = shape
        self.n = numpy.product(shape)
        self.axes = axes
        self._in  = numpy.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self._out = numpy.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self.forward_plan = fftw3.Plan(self._in, self._out, direction='forward', flags=flags, nthreads=nthreads)
        self.backward_plan = fftw3.Plan(self._in, self._out, direction='backward', flags=flags, nthreads=nthreads)
        self.dtype = var.COMPLEX_DTYPE

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
    def __init__(self, nsamples, description=None):
        AcquisitionModel.__init__(self, description, types=Tod)
        self.nsamples_array = numpy.array(nsamples, ndmin=1, dtype='int64')
        self.nsamples = tuple(self.nsamples_array)
        self.nsamples_tot = numpy.sum(self.nsamples_array)
        self.forward_plan = numpy.empty(self.nsamples_array.size, dtype='int64')
        self.backward_plan = numpy.empty(self.nsamples_array.size, dtype='int64')
        for i, n in enumerate(self.nsamples):
            tarray = numpy.empty(n, dtype=var.FLOAT_DTYPE)
            farray = numpy.empty(n, dtype=var.FLOAT_DTYPE)
            self.forward_plan[i] = fftw3.Plan(tarray, farray, direction='forward', flags=['measure'], realtypes=['halfcomplex r2c'], nthreads=1)._get_parameter()
            self.backward_plan[i] = fftw3.Plan(farray, tarray, direction='backward', flags=['measure'], realtypes=['halfcomplex c2r'], nthreads=1)._get_parameter()
            

    def direct(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        output_ = _smart_reshape(output, (numpy.product(input.shape[:-1]), input.shape[-1]))
        tmf.fft_plan(output_.T, self.nsamples_array, self.forward_plan)
        return output

    def transpose(self, input, reusein=False, reuseout=False):
        output = self.validate_input(input, reusein and reuseout)
        output_ = _smart_reshape(output, (numpy.product(input.shape[:-1]), input.shape[-1]))
        tmf.fft_plan(output_.T, self.nsamples_array, self.backward_plan)
        dest = 0
        for n in self.nsamples:
            output_[:,dest:dest+n] /= n
            dest += n
        return output

    def validate_shape(self, shape):
        if shape is None:
            return None
        nsamples = shape[-1]
        if nsamples != self.nsamples and nsamples != self.nsamples_tot:
            raise ValidationError("Invalid FFT size '"+str(nsamples)+"' instead of '"+str(self.nsamples)+"'.")
        return combine_sliced_shape(shape[0:-1], self.nsamples)


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
        self.diagonal /= var.mpi_comm.allreduce(numpy.max(self.diagonal), op=MPI.MAX)


#-------------------------------------------------------------------------------


class InterpolationLinear(Square):
    def __init__(self, mask, description=None):
        Square.__init__(self, description)
        self.mask = mask
    def direct(self, input, reusein=False, reuseout=False):
        input.mask = self.mask
        return interpolate_linear(input)
    def transpose(self, input, reusein=False, reuseout=False):
        raise NotImplementedError()


#-------------------------------------------------------------------------------


class AllReduce(Square):

    def direct(self, input, reusein=False, reuseout=False, op=MPI.SUM):
        input = self.validate_input(input, reusein and reuseout)
        input[:] = var.mpi_comm.allreduce(input, op=MPI.SUM)
        return input
    def transpose(self, input, reusein=False, reuseout=False):
        input = self.validate_input(input, reusein and reuseout)
        return input


#-------------------------------------------------------------------------------


def asacquisitionmodel(operator, description=None):
    if isinstance(operator, AcquisitionModel):
        return operator
    if _my_isscalar(operator):
        return Scalar(operator)
    if isinstance(operator, scipy.sparse.linalg.interface.LinearOperator):
        model = AcquisitionModel(description)
        model.direct    = lambda input, reusein=False, reuseout=False: operator.matvec(input)
        model.transpose = lambda input, reusein=False, reuseout=False: operator.rmatvec(input)
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
    cls = input.__class__
    while not issubclass(type(output), cls):
        cls = cls.__base__
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
        if curr.base is None or curr.base.dtype != input.dtype or curr.base.__class__ != input.__class__:
            return curr.reshape(shape)
        curr = curr.base


#-------------------------------------------------------------------------------


def _validate_input_unit(input, expected):
    if expected is None or input._unit is None:
        return
    for u,v in expected.items():
        if u not in input._unit or input._unit[u] != v:
            raise ValidationError("The input unit '" + input.unit + "' is incompatible with the required unit '" + \
                                  Quantity(1, expected).unit + "'.")
    return True


#-------------------------------------------------------------------------------


def _validate_output_unit(input, output, unitin, unitout, duout):
    if not hasattr(output, '_unit'):
        return
    if hasattr(input, '_unit'):
        output._unit = input._unit or unitin
        duout = duout or input.derived_units
    output.derived_units = duout
    if unitout is None:
        return
    if output._unit is None:
        output._unit = unitout
        return
    output._unit = _divide_unit(output._unit, unitin)
    output._unit = _multiply_unit(output._unit, unitout)
