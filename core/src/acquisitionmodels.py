from __future__ import division

import functools
import inspect
import numpy as np
import operator
import pyoperators

import tamasisfortran as tmf

from pyoperators import (Operator, IdentityOperator, DiagonalOperator,
                         BlockColumnOperator, BlockDiagonalOperator,
                         CompositionOperator, DistributionIdentityOperator,
                         MaskOperator, NumexprOperator)
from pyoperators.decorators import (linear, orthogonal, real, square, symmetric,
                                    unitary, inplace)
from pyoperators.utils import isscalar, tointtuple, openmp_num_threads
from pyoperators.utils.mpi import MPI, distribute_shape

from . import var
from .datatypes import Map, Tod
from .quantities import Quantity, _divide_unit, _multiply_unit
from .utils import diff, diffT, diffTdiff, shift
from .wcsutils import str2fitsheader

try:
    import fftw3
    MAX_FFTW_NUM_THREADS = 1 if fftw3.planning.lib_threads is None \
        else openmp_num_threads()
except:
    print('Warning: Library PyFFTW3 is not installed.')

__all__ = [
    'BlackBodyOperator',
    'CompressionAverage', # obsolete
    'CompressionAverageOperator',
    'ConvolutionOperator',
    'DistributionLocalOperator',
    'DdTdd', # Obsolete
    'DdTddOperator',
    'DiscreteDifference', # obsolete
    'DiscreteDifferenceOperator',
    'DownSamplingOperator',
    'FftOperator',
    'FftHalfComplexOperator',
    'InvNtt', # obsolete
    'InvNttOperator',
    'Masking', # obsolete
    'PackOperator',
    'PadOperator',
    'Projection', # obsolete
    'ProjectionOperator',
    'ResponseTruncatedExponential', # obsolete
    'ConvolutionTruncatedExponentialOperator',
    'RollOperator',
    'ShiftOperator',
    'SqrtInvNttOperator',
    'Unpacking', # obsolete
    'UnpackOperator',
    'PointingMatrix',
]

def block_diagonal(*partition_args, **keywords):
    """
    Class decorator that partitions an Operator along a specified axis.
    It adds a 'partition' keyed argument to the class constructor.

    Subclasses can also benefit from a decorated class, as long as they do not
    define a __init__ or __new__ method.
    
    Parameters
    ----------
    axisin: int
        Input partition axis.

    *partition_args: string varargs
       Specify the class' arguments of the __init__ method that  can be a
       sequence to be dispatched to the partitioned operators.

    Example
    -------
    >>> @block_diagonal('value', axisin=0)
    >>> @linear
    >>> class MyOp(Operator):
    >>>     def __init__(self, value, shapein=None):
    >>>         self.value = value
    >>>         Operator.__init__(self, lambda i,o: np.multiply(i,value,o),
    >>>                           shapein=shapein)
    >>>
    >>> op = MyOp([1,2,3], shapein=3, partitionin=(1,1,1))
    >>> op.todense()
    array([[ 1.,  0.,  0.],
           [ 0.,  2.,  0.],
           [ 0.,  0.,  3.]])

    >>> [o.value for o in op.operands]
    [1, 2, 3]
        
    """
    # the following way to extract keywords is unnecessary in Python3 (PEP 3102)
    if 'axisin' not in keywords:
        raise TypeError("Missing 'axisin' keyword.")

    axisin = keywords['axisin']
    del keywords['axisin']

    if len(keywords) > 0:
        raise TypeError('Invalid keyed argument.')

    def func(cls):

        @functools.wraps(cls.__new__)
        def partition_new(cls_, *args, **keywords):

            # get number of blocks through the partitionin keyword
            n1 = 0
            partitionin = None
            if 'partitionin' in keywords:
                if keywords['partitionin'] is not None:
                    partitionin = tointtuple(keywords['partitionin'])
                    n1 = len(partitionin)
                del keywords['partitionin']

            # get __init__ argument names, except self
            class_args, jk, jk, jk = inspect.getargspec(cls.__init_original__)
            class_args.pop(0)

            # get number of blocks through the arguments
            ns = [0 if isscalar(a) else len(a) for i,a in enumerate(args) \
                  if class_args[i] in partition_args] + \
                 [0 if isscalar(v) else len(v) for k,v in keywords.items() \
                  if k in partition_args]

            # required for copy.copy, which calls __new__ without arguments
            if n1 == 0 and len(ns) == 0:
                result = cls.__new_original__(cls_)
                if len(class_args) == 0:
                    cls.__init_original__(result)
                return result

            n2 = 0 if len(ns) == 0 else max(ns)
            if any(n not in (0, n2) for n in ns):
                raise ValueError('The partition variables do not have the same '
                                 'number of elements.')

            # check the two methods are compatible
            if n1 != 0 and n2 != 0 and n1 != n2:
                raise ValueError('The specified partitioning is incompatible wi'
                    'th the number of elements in the partition variables.')
            n = max(n1, n2, 1)

            # Implicit partition
            if partitionin is None:
                partitionin = n * (None,)

            # dispatch arguments
            n = len(partitionin)
            argss = tuple(tuple(a[i] if class_args[j] in partition_args and \
                not isscalar(a) else a for j,a in enumerate(args)) \
                for i in range(n))
            keyss = tuple(dict((k, v[i]) if k in partition_args and \
                not isscalar(v) else (k,v) for k,v in keywords.items()) \
                for i in range(n))
            
            # the input shapein/out describe the BlockDiagonalOperator
            def _reshape(s, p, a):
                s = list(tointtuple(s))
                s[a] = p
                return tuple(s)
            if 'shapein' in keywords:
                shapein = keywords['shapein']
                for keys, p in zip(keyss, partitionin):
                    keys['shapein'] = _reshape(shapein, p, axisin)
            if 'shapeout' in keywords:
                shapeout = keywords['shapeout']
                for keys, p in zip(keyss, partitionin):
                    keys['shapeout'] = _reshape(shapeout, p, axisin)

            # instantiate the partitioned operators
            ops = [cls.__new_original__(cls_, *a, **k) \
                   for a,k in zip(argss, keyss)]
            for o, a, k in zip(ops, argss, keyss):
                if isinstance(o, cls):
                    cls.__init_original__(o, *a, **k)

            return BlockDiagonalOperator(ops, partitionin=partitionin,
                                         axisin=axisin)

        @functools.wraps(cls.__init__)
        def partition_init(self, *args, **keys):
            # opearator has already been initialised
            pass

        cls.__new_original__ = staticmethod(cls.__new__)
        cls.__new__ = staticmethod(partition_new)
        cls.__init_original__ = cls.__init__
        cls.__init__ = partition_init
        return cls

    return func


class BlackBodyOperator(Operator):
    """
    Instanceless class whose __new__ method specialises in
    BlackBodyFixedTemperatureOperator and BlackBodyFreeTemperatureOperator.
    """
    def __new__(cls, wavelength=None, wavelength0=None, temperature=None,
                beta=0.):
        if temperature is not None:
            return BlackBodyFixedTemperatureOperator(temperature, wavelength,
                                                     wavelength0, beta)
        raise NotImplementedError()


@symmetric
class BlackBodyFixedTemperatureOperator(NumexprOperator):
    """
    Diagonal operator whose normalised values follow the Planck equation,
    optionnally modified by a power-law emissivity of given slope.
    """
    def __new__(cls, temperature=None, wavelength=None, wavelength0=None,
                beta=0.):
        if not isscalar(wavelength):
            return BlockColumnOperator([BlackBodyFixedTemperatureOperator(
                temperature, w, wavelength0, beta) for w in wavelength],
                new_axisout=0)
        return NumexprOperator.__new__(cls)

    def __init__(self, temperature, wavelength, wavelength0=None, beta=0.):
        """
        Parameters
        ----------
        temperature : float
            Black body temperature, in Kelvin.
        wavelength : float
            Wavelength, in meters, at which black body values will be
            computed.
        wavelength0 : float
            Reference wavelength, in meters, for which the operator returns 1.
        beta : float
            Slope of the emissivity (the spectrum is multiplied by nu**beta)
            The default value is 0 (non-modified black body).
        """
        self.temperature = np.array(temperature, float, copy=False)
        self.wavelength = float(wavelength)
        self.wavelength0 = float(wavelength0)
        self.beta = float(beta)
        c = 2.99792458e8
        h = 6.626068e-34
        k = 1.380658e-23
        if self.temperature.size == 1:
            coef = (self.wavelength0/self.wavelength)**(3+self.beta) * \
                   (np.exp(h*c/(self.wavelength0*k*self.temperature))-1) / \
                   (np.exp(h*c/(self.wavelength*k*self.temperature))-1)
            expr = 'coef * input'
            global_dict = {'coef':coef}
        else:
            coef1 = (self.wavelength0/self.wavelength)**(3+self.beta)
            coef2 = h*c/(self.wavelength0*k)
            coef3 = h*c/(self.wavelength*k)
            expr = 'coef1 * (exp(coef2/T) - 1) / (exp(coef3/T) - 1)'
            global_dict = {'coef1':coef1, 'coef2':coef2, 'coef3':coef3,
                           'T':self.temperature}
        NumexprOperator.__init__(self, expr, global_dict, dtype=float)


@block_diagonal('factor', axisin=-1)
@real
@linear
class CompressionOperator(Operator):
    """
    Abstract class for compressing the input signal.
    Compression is operated on the fastest axis.
    """

    def __new__(cls, factor=None, shapein=None, **keywords):
        if factor == 1:
            op = IdentityOperator(shapein=shapein, dtype=float, **keywords)
            op.factor = 1
            return op
        return Operator.__new__(cls)

    def __init__(self, factor, shapein=None, shapeout=None, **keywords):
        self.factor = int(factor)
        Operator.__init__(self, shapein=shapein, shapeout=shapeout, **keywords)

    def reshapein(self, shapein):
        return shapein[:-1] + (shapein[-1] // self.factor,)

    def reshapeout(self, shapeout):
        return shapeout[:-1] + (shapeout[-1] * self.factor,)

    def validatein(self, shapein):
        if shapein[-1] % self.factor != 0:
            raise ValueError("The input timeline size '{0}' is incompatible wit"
                "h the compression factor '{1}'.".format(shapein[-1],
                self.factor))

    def __str__(self):
        return super(CompressionOperator, self).__str__() + ' (x{0})'.format(
               self.factor)


class CompressionAverageOperator(CompressionOperator):
    """
    Compress the input signal by averaging blocks of specified size.
    """

    def direct(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.compression_average_direct(input_, ishape[0], ishape[1], istride,
            output_, oshape[1], ostride, self.factor)

    def transpose(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.compression_average_transpose(input_, ishape[0], ishape[1],
            istride, output_, oshape[1], ostride, self.factor)


class DownSamplingOperator(CompressionOperator):
    """
    Downsample the input signal by picking up one sample out of a number
    specified by the compression factor
    """

    def direct(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.downsampling_direct(input_, ishape[0], ishape[1], istride, output_,
            oshape[1], ostride, self.factor)

    def transpose(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.downsampling_transpose(input_, ishape[0], ishape[1], istride,
            output_, oshape[1], ostride, self.factor)
      

@real
@symmetric
class InvNttOperator(Operator):

    def __new__(cls, obs=None, method='uncorrelated', filename=None,
                **keywords):
        if obs is None:
            return Operator.__new__(cls)
        nsamples = obs.get_nsamples()
        method = method.lower()
        if method not in ('uncorrelated', 'uncorrelated python'):
            raise ValueError("Invalid method '{0}'.".format(method))

        filter = obs.get_filter_uncorrelated(filename=filename, **keywords)
        if filter.ndim == 2:
            filter = filter[np.newaxis,...]
        nfilters = filter.shape[0]
        if nfilters != 1 and nfilters != len(nsamples):
            raise ValueError("Incompatible number of filters '{0}'. Expected nu"
                             "mber is '{1}'".format(nfilters, len(nsamples)))
        ncorrelations = filter.shape[-1] - 1
        filter_length = np.asarray(2**np.ceil(np.log2(np.array(nsamples) + \
                                   2 * ncorrelations)), int)
        fft_filters = []
        for i, n in enumerate(nsamples):
            i = min(i, nfilters-1)
            fft_filter, status = tmf.fft_filter_uncorrelated(filter[i].T,
                                                             filter_length)
            if status != 0: raise RuntimeError()
            np.maximum(fft_filter, 0, fft_filter)
            fft_filters.append(fft_filter.T)
            norm = var.comm_tod.allreduce(max([np.max(f) for f in fft_filters]),
                                          op=MPI.MAX)
        for f in fft_filters:
            np.divide(f, norm, f)
            
        if method == 'uncorrelated':
            cls = InvNttUncorrelatedOperator
        else:
            cls = InvNttUncorrelatedPythonOperator
        #XXX should generate BlockDiagonalOperator with no duplicates...
        return BlockDiagonalOperator([cls(f, ncorrelations, n) \
                   for f, n in zip(fft_filters, nsamples)], axisin=-1)


@real
@symmetric
class SqrtInvNttOperator(InvNttOperator):
    def __init__(self, *args, **kw):
        invntt = InvNttOperator(*args, **kw)
        self._tocomposite(op.Composition, invntt.blocks)
        data = self.blocks[2].data
        np.sqrt(data, data)


@real
@symmetric
@inplace
class InvNttUncorrelatedOperator(Operator):

    def __init__(self, fft_filter, ncorrelations, nsamples,
                 fftw_flags=['measure', 'unaligned']):
        filter_length = fft_filter.shape[-1]
        array = np.empty(filter_length, dtype=var.FLOAT_DTYPE)
        self.fft_filter = fft_filter
        self.filter_length = filter_length
        self.left = filter_length - nsamples - ncorrelations
        self.right = ncorrelations
        self.fftw_flags = fftw_flags
        self.fplan = fftw3.Plan(array, direction='forward', flags=fftw_flags,
            realtypes=['halfcomplex r2c'], nthreads=1)
        self.bplan = fftw3.Plan(array, direction='backward', flags=fftw_flags,
            realtypes=['halfcomplex c2r'], nthreads=1)
        Operator.__init__(self, shapein=(fft_filter.shape[0], nsamples))
        
    def direct(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.invntt_uncorrelated(input_, ishape[1], istride, output_, ostride, self.fft_filter.T, self.fplan._get_parameter(), self.bplan._get_parameter(), self.left, self.right)
       

@real
@symmetric
class InvNttUncorrelatedPythonOperator(Operator):

    def __new__(cls, fft_filter=None, ncorrelations=None, nsamples=None):
        if fft_filter is None:
            return Operator.__new__(cls)
        filter_length = fft_filter.shape[-1]
        invntt = DiagonalOperator(fft_filter)
        fft = FftHalfComplexOperator(filter_length)
        padding = PadOperator(left=ncorrelations, right=filter_length - \
                              nsamples - ncorrelations)
        return padding.T * fft.T * invntt * fft * padding


@block_diagonal('left', 'right', axisin=-1)
@real
@linear
class PadOperator(Operator):
    """Pads before and after along the fast dimension of an ndarray."""

    def __init__(self, left=0, right=0., **keywords):
        left = int(left)
        right = int(right)
        if left < 0:
            raise ValueError('Left padding is not positive.')
        if right < 0:
            raise ValueError('Right padding is not positive.')
        self.left = left
        self.right = right
        Operator.__init__(self, **keywords)        
   
    def direct(self, input, output):
        if hasattr(output, 'mask'):
            output.mask = None
        right = -self.right if self.right != 0 else output.shape[-1]
        output[...,:self.left] = 0
        output[...,self.left:right] = input
        output[...,right:] = 0
   
    def transpose(self, input, output):
        if hasattr(output, 'mask'):
            output.mask = None
        right = -self.right if self.right != 0 else input.shape[-1]
        output[...] = input[...,self.left:right]

    def reshapein(self, shapein):
        return shapein[:-1] + (shapein[-1] + self.left + self.right,)
       
    def reshapeout(self, shapeout):
        return shapeout[:-1] + (shapeout[-1] - self.left - self.right,)


class PointingMatrix(np.ndarray):
    DTYPE = np.dtype([('value', np.float32), ('index', np.int32)])
    def __new__(cls, array, info=None, copy=True, ndmin=0):
        result = np.array(array, copy=copy, ndmin=ndmin).view(cls)
        result.info = info
        return result

    def __array_finalize__(self, obj):
        self.info = getattr(obj, 'info', None)

    def __getattr__(self, name):
        if self.dtype.names is not None and name in self.dtype.names:
            return self[name].view(np.ndarray)
        return super(PointingMatrix, self).__getattribute__(name)

    @classmethod
    def empty(cls, shape, info=None, verbose=True):
        if verbose:
            n = np.product(shape)
            if n > 0:
                print('Info: Allocating ' + str(np.product(shape) / 2**17) + \
                  ' MiB for the pointing matrix.')
        return PointingMatrix(np.empty(shape, dtype=cls.DTYPE), info=info,
                              copy=False)

    @classmethod
    def zeros(cls, shape, info=None, verbose=True):
        if verbose:
            n = np.product(shape)
            if n > 0:
                print('Info: Allocating ' + str(np.product(shape) / 2**17) + \
                  ' MiB for the pointing matrix.')
        result = PointingMatrix(np.zeros(shape, dtype=cls.DTYPE), info=info,
                                copy=False)
        result.index = -1
        return result

@real
@linear
class ProjectionOperator(Operator):
    """
    This class handles operations by the pointing matrix

    It is assumed that each pointing matrix row has at most 'npixels_per_sample'
    non-null elements, or in other words, an output element can only intersect
    a fixed number of input elements.

    Given an input vector x and an output vector y, y = P(x) translates into:
        y[i] = sum(P.matrix[i,:].value * x[P.matrix[i,:].index])

    If the input is not MPI-distributed unlike the output, the projection
    operator is automatically multiplied by the operator DistributionIdentity-
    Operator, to enable MPI reductions.

    If the input is MPI-distributed, this operator is automatically packed (see
    below) and multiplied by the operator DistributionLocalOperator, which
    takes the local input as argument.

    Arguments
    ---------
    input : pointing matrix or Observation
        It is recommended to use a pointing matrix or a list of pointing
        matrices instead of an observation:
            proj = ProjectionOperator(obs.get_pointing_matrix(...))
    packed : boolean
        If true, the input indices are renumbered to discard the input elements
        which are not projected onto the output. This feature is used to handle
        MPI-distributed input arrays.
    commin: mpi4py.Comm
        The input MPI communicator.
    commout: mpi4py.Comm
        The output MPI communicator.

    Attributes
    ----------
    header : pyfits.Header
        The map FITS header
    matrix : composite dtype ('value', 'f4'), ('index', 'i4')
        The pointing matrix. 
    """

    def __new__(cls, input=None, method=None, header=None, resolution=None,
                npixels_per_sample=0, units=None, derived_units=None,
                downsampling=False, packed=False, shapein=None, shapeout=None,
                commin=MPI.COMM_WORLD, commout=MPI.COMM_WORLD, **keywords):

        # check if there is only one pointing matrix
        isobservation = hasattr(input, 'get_pointing_matrix')
        if isobservation or input is None:
            if hasattr(input, 'slice'):
                nmatrices = len(input.slice)
            else:
                nmatrices = 1
        else:
            if isinstance(input, PointingMatrix):
                input = (input,)
            if any(not isinstance(i, PointingMatrix) for i in input):
                raise TypeError('The input is not a PointingMatrix, (nor a sequ'
                                'ence of).')
            nmatrices = len(input)

        ismapdistributed = commin is not None and commin.size > 1
        istoddistributed = commout is not None and commout.size > 1

        # the output class will be a ProjectionOperator
        if nmatrices == 1 and not ismapdistributed and not istoddistributed \
           and not packed:
            instance = Operator.__new__(cls)
            return instance

        # the input has multiple pointing matrices or will be multiplied by an
        # operator
        if isobservation:
            if header is None:
                if not hasattr(input, 'get_map_header'):
                    raise AttributeError("No map header has been specified and "
                        "the observation has no 'get_map_header' method.")
                header = input.get_map_header(resolution=resolution,
                                              downsampling=downsampling)
        else:
            headers = [i.info['header'] for i in input]
            shapein_headers = [tuple([h['NAXIS' + str(i+1)] for i in \
                               range(h['NAXIS'])])[::-1] for h in headers]
            if any(s != shapein_headers[0] for s in shapein_headers):
                raise ValueError('The pointing matrices do not have the same in'
                    "put shape: {0}.".format(', '.join(str(shapein_headers))))
            if header is None:
                header = headers[0]

        if isinstance(header, str):
            header = str2fitsheader(header)

        shapein = shapein if shapein is not None else tuple([header['NAXIS' + \
            str(i+1)] for i in range(header['NAXIS'])])[::-1]

        # get the block column operator
        if isobservation:
            input = input.get_pointing_matrix(header, npixels_per_sample,
                        method=method, downsampling=downsampling)
        if nmatrices == 1:
            input = (input,)

        operands = [ProjectionOperator(p, packed=False, shapein=shapein,
                    commin=None, commout=None, **keywords) for p in input]
        result = BlockColumnOperator(operands, axisout=-1)

        if nmatrices > 1:
            # get the global mask
            def get_mask(output=None):
                if output is None:
                    #FIXME shapein_global = collect_shape(shapein)
                    shapein_global = shapein
                    output = Map.ones(shapein_global, dtype=np.bool8,
                                      header=header)
                for p in result.operands:
                    p.get_mask(output)
                return output
            result.header = header
            result.get_mask = get_mask

        # handle packed or distributed input and output
        return cls.distribute(result, commin, commout, packed)

    def __init__(self, input, method=None, header=None, resolution=None,
                 npixels_per_sample=0, units=None, derived_units=None,
                 downsampling=False, packed=False, shapein=None, shapeout=None,
                 commin=MPI.COMM_WORLD, commout=MPI.COMM_WORLD, **keywords):

        if hasattr(input, 'get_pointing_matrix'):
            if header is None:
                if not hasattr(input, 'get_map_header'):
                    raise AttributeError("No map header is specified and the ob"
                        "servation has no 'get_map_header' method.")
                header = input.get_map_header(resolution=resolution,
                                              downsampling=downsampling)
            elif isinstance(header, str):
                header = str2fitsheader(header)

            commout = input.comm_tod
            matrix = input.get_pointing_matrix(header, npixels_per_sample,
                method=method, downsampling=downsampling)
        else:
            if not isinstance(input, PointingMatrix):
                input = input[0]
            if input.ndim != 3:
                raise ValueError('The input pointing matrix has not 3 dimension'
                                 's.')
            matrix = input

        if header is None:
            header = matrix.info.get('header', None)
        if units is None:
            units = matrix.info.get('units', None)
        if derived_units is None:
            derived_units = matrix.info.get('derived_units', None)

        if units is None:
            units = ('', '')
        if derived_units is None:
            derived_units = ({}, {})

        unit = _divide_unit(Quantity(1, units[0])._unit,
                            Quantity(1, units[1])._unit)

        ndetectors, nsamples, npixels_per_sample = matrix.shape

        # validate shapein
        shapein_expected = None if header is None else tuple([header['NAXIS' + \
            str(i+1)] for i in range(header['NAXIS'])])[::-1]
        if shapein is None and shapein_expected is None:
            raise ValueError("The pointing matrix itself does not contain infor"
                "mation about the input shape. Use the 'header' or 'shapein' ke"
                "ywords to specify it.")
        if shapein is None:
            shapein = shapein_expected
        elif shapein_expected is not None and shapein != shapein_expected:
            raise ValueError("The specified input shape '{0}' is incompatible w"
                "ith that of the header '{1}'.".format(shapein,
                shapein_expected))

        # validate shapeout
        shapeout_expected = (ndetectors, nsamples)
        if shapeout is None:
            shapeout = shapeout_expected
        elif shapeout != shapeout_expected:
            raise ValueError("The specified output shape '{0}' is incompatible "
                "with that of the pointing matrix '{1}'.".format(shapeout,
                shapeout_expected))

        # attribute setting functions
        def set_attrin(attr):
            if header is not None:
                attr['header'] = self.header.copy()
            unitout = attr['_unit'] if '_unit' in attr else {}
            if unitout:
                attr['_unit'] = _divide_unit(unitout, unit)
            if derived_units[1] is not None:
                attr['_derived_units'] = derived_units[1]

        def set_attrout(attr):
            attr['header'] = None
            unitin = attr['_unit'] if '_unit' in attr else {}
            if unitin:
                attr['_unit'] = _multiply_unit(unitin, unit)
            if derived_units[0] is not None:
                attr['_derived_units'] = derived_units[0]

        # f2py doesn't accept zero-sized opaque arguments
        if matrix.size == 0:
            _matrix = np.empty(1, np.int64)
        else:
            _matrix = matrix.ravel().view(np.int64)

        # initialise operator
        Operator.__init__(self, shapein=shapein, shapeout=shapeout,
                          attrin=set_attrin, attrout=set_attrout,
                          classin=Map, classout=Tod, dtype=var.FLOAT_DTYPE,
                          **keywords)
        self.matrix = matrix
        self._matrix = _matrix
        self.commin = commin
        self.commout = commout
        self.header = header
        self.ndetectors = ndetectors
        self.nsamples = nsamples
        self.npixels_per_sample = npixels_per_sample
        self.unit = unit
        self.duout, self.duin = derived_units

    @staticmethod
    def distribute(self, commin, commout, packed, **keywords):

        ismapdistributed = commin is not None and commin.size > 1
        istoddistributed = commout is not None and commout.size > 1

        if not ismapdistributed and not istoddistributed and not packed:
            return self

        if packed or ismapdistributed:
            mask = self.get_mask()
            shapein = int(np.sum(~mask))
            ops = (self,) if isinstance(self, ProjectionOperator) else \
                  self.operands
            for op in ops:
                tmf.pointing_matrix_pack(self._matrix, mask.view(np.int8).T,
                    self.npixels_per_sample, self.nsamples * self.ndetectors)
                Operator.__init__(op, shapein=shapein, shapeout=self.shapeout,
                                  attrin=op.attrin, attrout=op.attrout,
                                  classin=op.classin, classout=op.classout,
                                  dtype=self.dtype, **keywords)
                op.validatein(tointtuple(shapein))
                op.packed = True

            if ismapdistributed:
                self *= DistributionLocalOperator(mask)
            else:
                self *= PackOperator(mask)
                self.operands[0].validatein(tointtuple(shapein))

        else:
            self *= DistributionIdentityOperator(commout=self.commout)

        s = self.operands[0]
        self.header = s.header
        self.ndetectors = s.ndetectors
        self.npixels_per_sample = s.npixels_per_sample

        return self

    def direct(self, input, output):
        if not output.flags.contiguous:
            if pyoperators.memory.verbose:
                print 'Optimize me: Projection output is not contiguous.'
            output_ = np.empty_like(output)
        else:
            output_ = output
        tmf.pointing_matrix_direct(self._matrix, input.T, output_.T,
                                   self.npixels_per_sample)
        if not output.flags.contiguous:
            output[...] = output_

    def transpose(self, input, output, operation=None):
        if not input.flags.contiguous:
            if pyoperators.memory.verbose:
                print 'Optimize me: Projection.T input is not contiguous.'
            input_ = np.ascontiguousarray(input)
        else:
            input_ = input
        if operation is None:
            output[...] = 0
        elif operation is not operator.__iadd__:
            raise ValueError('Invalid operation.')
        tmf.pointing_matrix_transpose(self._matrix, input_.T, output.T,
                                      self.npixels_per_sample)

    def get_mask(self, output=None):
        shapein_global = self.shapein
        #XXX FIXME: shapein_global = collect_shape(self.shapein)
        if output is None:
            output = Map.ones(shapein_global, dtype=np.bool8,
                              header=self.header)
        elif output.dtype != bool or output.shape != shapein_global:
            raise ValueError('The argument to store the mask is incompatible.')
        tmf.pointing_matrix_mask(self._matrix, output.view(np.int8).T, 
            self.npixels_per_sample, self.nsamples * self.ndetectors)
        return output

    def get_ptp(self):
        npixels = np.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._matrix, self.npixels_per_sample,
            self.nsamples_tot, self.ndetectors, npixels).T

    def intersects(self, fitscoords, axis=None, f=None):
        """
        Return True if a map pixel is seen by the pointing matrix.

        """

        naxes = self.header['NAXIS']
        axes = [1] + [self.header['NAXIS'+str(i+1)] for i in range(naxes-1)]
        index = np.array(fitscoords) - 1
        index = int(np.sum(index * np.cumproduct(axes)))

        shape = self.matrix.shape
        if f is not None:
            out = f(self._matrix, index, shape[-1], shape[-2], shape[-3])
        elif axis is None:
            out = tmf.operators.pmatrix_intersects(self._matrix, index,
                      shape[-1], shape[-2], shape[-3])
        elif axis == 0 or axis == -3:
            out = tmf.operators.pmatrix_intersects_axis3(self._matrix, index,
                      shape[-1], shape[-2], shape[-3])
        elif axis == 1 or axis == -2:
            out = tmf.operators.pmatrix_intersects_axis2(self._matrix, index,
                      shape[-1], shape[-2], shape[-3])
        else:
            raise ValueError('Invalid axis.')

        if isinstance(out, np.ndarray):
            out = out.view(bool)
        else:
            out = bool(out)
        return out


@real
@linear
@square
@inplace
class ConvolutionTruncatedExponentialOperator(Operator):
    """
    Apply a truncated exponential response to the signal

    Parameters
    ==========
    
    tau: number
        Time constant divided by the signal sampling period.
    """
    
    def __init__(self, tau, **keywords):
        """
        """
        if hasattr(tau, 'SI'):
            tau = tau.SI
            if tau.unit != '':
                raise ValueError('The time constant must be dimensionless.')
        self.tau = np.array(tau, dtype=var.FLOAT_DTYPE, ndmin=1)
        Operator.__init__(self, **keywords)

    def direct(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.convolution_trexp_direct(input_, ishape[0], ishape[1], istride,
            output_, oshape[1], ostride, self.tau)

    def transpose(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.convolution_trexp_transpose(input_, ishape[0], ishape[1], istride,
            output_, oshape[1], ostride, self.tau)


@real
@linear
@square
@inplace
class DiscreteDifferenceOperator(Operator):
    """
    Discrete difference operator.

    Calculate the nth order discrete difference along given axis.
    """

    def __init__(self, axis=-1, comm=None, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        self.axis = axis
        self.comm = comm or var.comm_map
        self.set_rule('.T.', self._rule_ddtdd, CompositionOperator)

    def direct(self, input, output):
        diff(input, output, self.axis, comm=self.comm)

    def transpose(self, input, output):
        diffT(input, output, self.axis, comm=self.comm)

    @staticmethod
    def _rule_ddtdd(dT, self):
        return DdTddOperator(self.axis, comm=self.comm)


@real
@symmetric
@inplace
class DdTddOperator(Operator):
    """Calculate operator dX.T dX along a given axis."""

    def __init__(self, axis=-1, scalar=1., comm=None, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        self.__name__ = (str(scalar) + ' ' if scalar != 1 else '') + \
            self.__name__
        self.axis = axis
        self.scalar = scalar
        self.comm = comm or var.comm_map        
        self.set_rule('{HomothetyOperator}.', lambda o,s: DdTddOperator(
                      self.axis, o.data * s.scalar), CompositionOperator)

    def direct(self, input, output):
        diffTdiff(input, output, self.axis, self.scalar, comm=self.comm)


@real
@linear
class DistributionLocalOperator(Operator):
    """
    Scatter a distributed map to different MPI processes under the control of a
    local non-distributed mask.
    """

    def __init__(self, mask, operator=MPI.SUM, comm=None, **keywords):
        if comm is None:
            comm = var.comm_map
        shapeout = (int(np.sum(~mask)),)
        shapein = distribute_shape(mask.shape, comm)
        attrin = { 'comm':comm, 'shape_global':mask.shape }
        attrout = { 'comm':MPI.COMM_SELF, 'shape_global':shapeout}
        Operator.__init__(self, classin=Map, shapein=shapein,
                                shapeout=shapeout, attrin=attrin,
                                attrout=attrout, **keywords)
        self.comm = comm
        self.mask = mask
        self.operator = operator

    def direct(self, input, output):
        status = tmf.operators_mpi.allscatterlocal(input.T, self.mask.view(
            np.int8).T, output.T, self.comm.py2f())
        if status == 0: return
        if status < 0:
            raise RuntimeError('Incompatible sizes.')
        raise MPI.Exception(status)

    def transpose(self, input, output):
        status = tmf.operators_mpi.allreducelocal(input.T, self.mask.view(
            np.int8).T, output.T, self.operator.py2f(), self.comm.py2f())
        if status == 0: return
        if status < 0:
            raise RuntimeError('Incompatible mask.')
        raise MPI.Exception(status)


@real
@linear
@inplace
class PackOperator(Operator):
    """
    Convert an nd array in a 1d map, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, **keywords):
        mask = np.array(mask, np.bool8)
        Operator.__init__(self, shapein=mask.shape, shapeout=np.sum(mask == 0),
                          dtype=var.FLOAT_DTYPE, **keywords)
        self.mask = mask

    def direct(self, input, output):
        tmf.packing(input.ravel(), self.mask.view(np.int8).ravel(),
                    output.ravel())

    def transpose(self, input, output):
        tmf.unpacking(input.ravel(), self.mask.view(np.int8).ravel(),
                      output.ravel())


@real
@linear
@inplace
class UnpackOperator(Operator):
    """
    Convert 1d map into an nd array, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, **keywords):
        mask = np.array(mask, np.bool8)
        Operator.__init__(self, shapein=np.sum(mask == 0), shapeout=mask.shape,
                          dtype=var.FLOAT_DTYPE, **keywords)
        self.mask = mask

    def direct(self, input, output):
        tmf.unpacking(input.ravel(), self.mask.view(np.int8).ravel(),
                      output.ravel())

    def transpose(self, input, output):
        tmf.packing(input.ravel(), self.mask.view(np.int8).ravel(),
                    output.ravel())


@real
@linear
@square
@inplace
class ShiftOperator(Operator):

    def __init__(self, n, axis=None, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        if axis is None:
            axis = -1 if isscalar(n) else range(-len(n), 0)
        self.axis = (axis,) if isscalar(axis) else tuple(axis)
        self.n = (n,) * len(self.axis) if isscalar(n) else \
                 tuple([np.array(m,int) for m in n])
        if len(self.axis) != len(self.n):
            raise ValueError('There is a mismatch between the number of axes an'
                             'd offsets.')

    def direct(self, input, output):
        shift(input, output, self.n[0], self.axis[0])
        for n, axis in zip(self.n, self.axis)[1:]:
            shift(output, output, n, axis)

    def transpose(self, input, output):
        shift(input, output, -self.n[0], self.axis[0])
        for n, axis in zip(self.n, self.axis)[1:]:
            shift(output, output, -n, axis)


@orthogonal
@inplace
class RollOperator(Operator):
    
    def __init__(self, n, axis=None, **keywords):
        Operator.__init__(self, **keywords)
        if axis is None:
            axis = -1 if isscalar(n) else range(-len(n), 0)
        self.axis = (axis,) if isscalar(axis) else tuple(axis)
        self.n = (n,) * len(self.axis) if isscalar(n) else tuple(n)
        if len(self.axis) != len(self.n):
            raise ValueError('There is a mismatch between the number of axes an'
                             'd offsets.')

    def direct(self, input, output):
        output[...] = np.roll(input, self.n[0], axis=self.axis[0])
        for n, axis in zip(self.n, self.axis)[1:]:
            output[...] = np.roll(output, n, axis=axis)

    def transpose(self, input, output):
        output[...] = np.roll(input, -self.n[0], axis=self.axis[0])
        for n, axis in zip(self.n, self.axis)[1:]:
            output[...] = np.roll(output, -n, axis=axis)


@unitary
@inplace
class FftOperator(Operator):
    """
    Performs complex fft
    """

    def __init__(self, shape, flags=['measure'], nthreads=None, **keywords):
        Operator.__init__(self, shapein=shape,
                                dtype=var.COMPLEX_DTYPE, **keywords)
        nthreads = min(nthreads or openmp_num_threads(), MAX_FFTW_NUM_THREADS)
        self.n = np.product(shape)
        self._in  = np.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self._out = np.zeros(shape, dtype=var.COMPLEX_DTYPE)
        self.forward_plan = fftw3.Plan(self._in, self._out, direction='forward',
                                       flags=flags, nthreads=nthreads)
        self.backward_plan= fftw3.Plan(self._in, self._out,direction='backward',
                                       flags=flags, nthreads=nthreads)

    def direct(self, input, output):
        self._in[...] = input
        fftw3.execute(self.forward_plan)
        output[...] = self._out

    def transpose(self, input, output):
        self._in[...] = input
        fftw3.execute(self.backward_plan)
        output[...] = self._out / self.n

@real
@orthogonal
@inplace
class FftHalfComplexOperator(Operator):
    """
    Performs real-to-half-complex fft
    """

    def __init__(self, size, fftw_flags=['measure', 'unaligned'], **keywords):
        size = int(size)
        array1 = np.empty(size, dtype=var.FLOAT_DTYPE)
        array2 = np.empty(size, dtype=var.FLOAT_DTYPE)
        fplan = fftw3.Plan(array1, array2, direction='forward', flags= \
            fftw_flags, realtypes=['halfcomplex r2c'], nthreads=1)
        bplan = fftw3.Plan(array1, array2, direction='backward', flags=\
            fftw_flags, realtypes=['halfcomplex c2r'], nthreads=1)
        ifplan = fftw3.Plan(array1, direction='forward', flags=fftw_flags,
            realtypes=['halfcomplex r2c'], nthreads=1)
        ibplan = fftw3.Plan(array1, direction='backward', flags=fftw_flags,
            realtypes=['halfcomplex c2r'], nthreads=1)
        self.size = size
        self.fftw_flags = fftw_flags
        self.fplan = fplan
        self.bplan = bplan
        self.ifplan = ifplan
        self.ibplan = ibplan
        Operator.__init__(self, **keywords)

    def direct(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        if self.same_data(input, output):
            tmf.fft_plan_inplace(input_, ishape[0], ishape[1], istride,
                                 self.ifplan._get_parameter())
            return
        output_, oshape, ostride = _ravel_strided(output)
        tmf.fft_plan_outplace(input_, ishape[0], ishape[1], istride, output_,
                     oshape[1], ostride, self.fplan._get_parameter())

    def transpose(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        if self.same_data(input, output):
            tmf.fft_plan_inplace(input_, ishape[0], ishape[1], istride,
                                 self.ibplan._get_parameter())
        else:
            output_, oshape, ostride = _ravel_strided(output)
            tmf.fft_plan_outplace(input_, ishape[0], ishape[1], istride,
                output_, oshape[1], ostride, self.bplan._get_parameter())
        output /= self.size

    def validatein(self, shapein):
        if shapein[-1] != self.size:
            raise ValueError("Invalid input dimension '{0}'. Expected dimension"
                             " is '{1}'".format(shapein[-1], self.size))

@linear
@square
@real
@inplace
class ConvolutionOperator(Operator):
    def __init__(self, shape, kernel, flags=['measure'], nthreads=None,
                 dtype=None, **keywords):

        dtype = dtype or var.get_default_dtype(kernel)
        kernel = np.array(kernel, dtype, copy=False)
        Operator.__init__(self, shapein=shape, dtype=dtype, **keywords)

        ndim = len(self.shapein)
        if ndim != kernel.ndim:
            raise ValueError("The kernel dimension '" + str(kernel.ndim) + "' i"
                "s incompatible with the specified shape '" + str(ndim) + "'.")

        # if the kernel is larger than the image, we don't crop it since it
        # might affect normalisation of the kernel
        if any([ks > s for ks,s in zip(kernel.shape, shape)]):
            raise ValueError('The kernel must not be larger than the input.')

        nthreads = min(nthreads or openmp_num_threads(), MAX_FFTW_NUM_THREADS)

        ker_origin = (np.array(kernel.shape)-1) / 2
        if any([int(o) != o for o in ker_origin]):
            raise ValueError('Kernel with even dimension is not yet handled.')

        # pad kernel with zeros
        ker_slice = [ slice(0,s) for s in kernel.shape ]
        kernel, kernel[ker_slice] = np.zeros(shape, kernel.dtype), kernel
        for axis, o in enumerate(ker_origin):
            kernel = np.roll(kernel, int(-o), axis=axis)
        
        # set up fftw plans
        self._in = np.empty(shape, kernel.dtype)
        self._out = np.empty(shape, complex)
        self.forward_plan = fftw3.Plan(self._in, self._out,
                                       direction='forward',
                                       flags=flags,
                                       nthreads=nthreads)
        self.backward_plan = fftw3.Plan(self._out, self._in,
                                        direction='backward',
                                        flags=flags,
                                        nthreads=nthreads)

        # FT kernel
        self._in[:] = kernel
        fftw3.execute(self.forward_plan)
        self.kernel = self._out / kernel.size
        
    def direct(self, input, output):
        self._in[...] = input
        fftw3.execute(self.forward_plan)
        self._out *= self.kernel
        fftw3.execute(self.backward_plan)
        output[...] = self._in
        
    def transpose(self, input, output):
        self._in[...] = input
        fftw3.execute(self.forward_plan)
        # avoid temp array in 'self._out *= self.kernel.conjugate()'
        tmf.multiply_conjugate_complex(self._out.ravel(), self.kernel.ravel(),
                                       self._out.ravel())
        fftw3.execute(self.backward_plan)
        output[...] = self._in


def _ravel_strided(array):
    # array_ = array.reshape((-1,array.shape[-1]))
    #XXX bug in numpy 1.5.1: reshape with -1 leads to bogus strides
    array_ = array.reshape(reduce(lambda x,y:x*y, array.shape[0:-1], 1),
                           array.shape[-1])
    n2, n1 = array_.shape
    if array_.flags.c_contiguous:
        return array_.ravel(), array_.shape, n1
    if not array_[0,:].flags.c_contiguous:
        raise RuntimeError()
    s2, s1 = array_.strides
    s2 //= s1
    base = array
    while not base.flags.c_contiguous:
        base = base.base
    start = (array.__array_interface__['data'][0] - \
             base.__array_interface__['data'][0]) // s1
    stop = start + (n2 - 1) * s2 + n1
    flat = base.ravel()[start:stop]
    return flat, array_.shape, s2

def obsolete(cls):
    @functools.wraps(cls.__base__.__init__)
    def init(self, *args, **keywords):
        print("The class '{0}' has been renamed. Use class '{1}' instead." \
              .format(cls.__name__, cls.__base__.__name__))
        cls.__init_original__(self, *args, **keywords)
    cls.__init_original__ = cls.__base__.__init__
    cls.__init__ = init
    return cls

        
@obsolete
class CompressionAverage(CompressionAverageOperator):
    pass
@obsolete
class DiscreteDifference(DiscreteDifferenceOperator):
    pass
@obsolete
class DdTdd(DdTddOperator):
    pass
@obsolete
class InvNtt(InvNttOperator):
    pass
@obsolete
class Masking(MaskOperator):
    pass
@obsolete
class Packing(PackOperator):
    pass
@obsolete
class Projection(ProjectionOperator):
    pass
@obsolete
class ResponseTruncatedExponential(ConvolutionTruncatedExponentialOperator):
    pass
@obsolete
class Unpacking(UnpackOperator):
    pass
