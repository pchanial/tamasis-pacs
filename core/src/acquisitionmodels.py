from __future__ import division

import functools
import inspect
import numpy as np
import operator

import tamasisfortran as tmf

try:
    import fftw3
except:
    print('Warning: Library PyFFTW3 is not installed.')

from mpi4py import MPI

from operators import (Operator, IdentityOperator, DiagonalOperator,
                       ExpansionOperator, PartitionOperator)
from operators.decorators import (idempotent, linear, orthogonal, real, square,
                                  symmetric, unitary, inplace)
from operators.utils import isscalar, tointtuple, openmp_num_threads

from . import var

from .datatypes import Map, Tod
from .observations import Observation
from .utils import diff, diffT, diffTdiff, shift
from .mpiutils import split_shape, split_work
from .wcsutils import str2fitsheader

__all__ = [
    'CompressionAverage',
    'Convolution',
    'DistributionGlobal',
    'DistributionLocal',
    'DdTdd',
    'DiscreteDifference',
    'DownSampling',
    'Fft',
    'FftHalfComplex',
    'InvNtt',
    'Masking',
    'Packing',
    'Padding',
    'Projection',
    'Projection_old',
    'ResponseTruncatedExponential',
    'RollOperator',
    'ShiftOperator',
    'SqrtInvNtt',
    'Unpacking',
]

def partitioned(*partition_args, **keywords):
    """
    Class decorator that partitions an Operator along a specified axis.
    It adds a 'partition' keyed argument to the class constructor.
    
    Parameters
    ----------
    axis: int
        Partition axis of the input (default is 0)

    *partition_args: string varargs
       Specify the class' arguments of the __init__ method that  can be a
       sequence to be dispatched to the partitioned operators.

    Example
    -------
    >>> @partitioned('value')
    >>> @linear
    >>> class MyOp(Operator):
    >>>     def __init__(self, value, shapein=None):
    >>>         self.value = value
    >>>         Operator.__init__(self, lambda i,o: np.multiply(i,value,o),
    >>>                           shapein=shapein)
    >>>
    >>> op = MyOp([1,2,3], shapein=3, partition=(1,1,1))
    >>> op.todense()
    array([[ 1.,  0.,  0.],
           [ 0.,  2.,  0.],
           [ 0.,  0.,  3.]])

    >>> [o.value for o in op.operands]
    [1, 2, 3]
        
    """
    # the following way to extract keywords is unnecessary in Python3 (PEP 3102)
    if 'axis' in keywords:
        axisin = keywords['axis']
        del keywords['axis']
    else:
        axisin = 0
    if len(keywords) > 0:
        raise TypeError('Invalid keyed argument.')

    def func(cls):

        @functools.wraps(cls.__new__)
        def partition_new(cls_, *args, **keywords):
            if 'partition' in keywords and keywords['partition'] is not None:
                partitionin = tointtuple(keywords['partition'])
                del keywords['partition']
            else:
                inst = cls.__new_original__(cls_, *args, **keywords)
                cls_.__init__ = partition_init
                return inst

            # dispatch arguments
            n = len(partitionin)
            class_args, jk, jk, jk = inspect.getargspec(cls_.__init_original__)
            class_args.pop(0)
            argss = tuple(tuple(a[i] if class_args[j] in partition_args and \
                not isscalar(a) else a for j,a in enumerate(args)) \
                for i in range(n))
            keyss = tuple(dict((k, v[i]) if k in partition_args and \
                not isscalar(v) else (k,v) for k,v in keywords.iteritems()) \
                for i in range(n))
            
            # the input shapein/out describe the PartitionOperator
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
            for o,a,k in zip(ops, argss, keyss):
                if isinstance(o, cls):
                    o.__init_original__(*a, **k)
            if len(ops) == 1:
                return ops[0]

            return PartitionOperator(ops, partitionin=partitionin,
                                     axisin=axisin)

        @functools.wraps(cls.__init__)
        def partition_init(self, *args, **keywords):
            if 'partition' in keywords:
                partitionin = tointtuple(keywords['partition'])
                if partitionin is not None and len(partitionin) == 1:
                    return
                del keywords['partition']
            cls.__init_original__(self, *args, **keywords)

        cls.__new_original__ = staticmethod(cls.__new__)
        cls.__new__ = staticmethod(partition_new)
        cls.__init_original__ = cls.__init__
        cls.__init__ = partition_init
        return cls

    return func


@partitioned('factor', axis=-1)
@real
@linear
class Compression(Operator):

    """
    Abstract class for compressing the input signal.
    Compression is operated on the fastest axis.
    """

    def __new__(cls, factor, shapein=None, **keywords):
        if factor == 1:
            op = IdentityOperator(shapein=shapein, **keywords)
            op.factor = 1
            return op
        return super(Compression, cls).__new__(cls, factor, shapein, **keywords)

    def __init__(self, factor, shapein=None, shapeout=None, **keywords):
        self.factor = int(factor)
        Operator.__init__(self, shapein=shapein, shapeout=shapeout, **keywords)

    def reshapein(self, shapein):
        if shapein is None:
            return None
        if shapein[-1] % self.factor != 0:
            raise ValueError("The input timeline size '{0}') is incompatible wi"
                "th the compression factor '{1}'".format(shapein[-1],
                self.factor))
        shapeout = list(shapein)
        shapeout[-1] = shapein[-1] // self.factor
        return shapeout

    def reshapeout(self, shapeout):
        if shapeout is None:
            return None
        shapeout = list(shapeout)
        shapeout[-1] = shapeout[-1] * self.factor
        return  shapeout

    def __str__(self):
        return super(Compression, self).__str__() + ' (x{})'.format(self.factor)


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    """

    def __str__(self):
        return super(CompressionAverage, self).__str__() + ' (x{})'.format(self.factor)

    def direct(self, input, output):
        if hasattr(output, 'mask'):
            output.mask = None
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.compression_average_direct(input_, ishape[0], ishape[1], istride,
            output_, oshape[1], ostride, self.factor)

    def transpose(self, input, output):
        if hasattr(output, 'mask'):
            output.mask = None
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.compression_average_transpose(input_, ishape[0], ishape[1],
            istride, output_, oshape[1], ostride, self.factor)


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number
    specified by the compression factor
    """

    def direct(self, input, output):
        if hasattr(output, 'mask'):
            output.mask = None
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.downsampling_direct(input_, ishape[0], ishape[1], istride, output_,
            oshape[1], ostride, self.factor)

    def transpose(self, input, output):
        if hasattr(output, 'mask'):
            output.mask = None
        input_, ishape, istride = _ravel_strided(input)
        output_, oshape, ostride = _ravel_strided(output)
        tmf.downsampling_transpose(input_, ishape[0], ishape[1], istride,
            output_, oshape[1], ostride, self.factor)
      

@real
@symmetric
class InvNtt(Operator):

    def __new__(cls, obs, method='uncorrelated', filename=None, **keywords):
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
            cls = InvNttUncorrelated
        else:
            cls = InvNttUncorrelatedPython
        #XXX should generate partitionoperator with no duplicates...
        return PartitionOperator([cls(f, ncorrelations, n) \
                   for f, n in zip(fft_filters, nsamples)], axisin=-1)


@real
@symmetric
class SqrtInvNtt(InvNtt):
    def __init__(self, *args, **kw):
        invntt = InvNtt(*args, **kw)
        self._tocomposite(op.Composition, invntt.blocks)
        data = self.blocks[2].data
        np.sqrt(data, data)


@real
@symmetric
@inplace
class InvNttUncorrelated(Operator):

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
class InvNttUncorrelatedPython(Operator):

    def __new__(cls, fft_filter, ncorrelations, nsamples):
        filter_length = fft_filter.shape[-1]
        invntt = DiagonalOperator(fft_filter)
        fft = FftHalfComplex(filter_length)
        padding = Padding(left=ncorrelations, right=filter_length - nsamples - \
                              ncorrelations)
        return padding.T * fft.T * invntt * fft * padding


@partitioned('left', 'right', axis=-1)
@real
@linear
class Padding(Operator):
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
        if shapein is None:
            return None
        shapeout = list(shapein)
        shapeout[-1] += self.left + self.right
        return shapeout
       
    def reshapeout(self, shapeout):
        if shapeout is None:
            return None
        shapein = list(shapeout)
        shapein[-1] -= self.left + self.right
        return shapein


#-------------------------------------------------------------------------------


@real
@linear
class Projection_old(Operator):
    """
    This class handles operations by the pointing matrix

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

        if header is None:
            header = observation.get_map_header(resolution=resolution,
                                                oversampling=oversampling)
        elif isinstance(header, str):
            header = str2fitsheader(header)
        self.header = header

        pmatrix, self.method, self.ndetectors, nsamples, \
        self.npixels_per_sample, (self.unitout, self.unitin), \
        (self.duout, self.duin) = \
            observation.get_pointing_matrix(header,
                                            npixels_per_sample,
                                            method=method,
                                            oversampling=oversampling)

        self.nsamples = nsamples
        self._pmatrix = pmatrix
        if self.npixels_per_sample == 0:
            pmatrix = np.empty(0, dtype=np.int64)
        self.pmatrix = pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]) \
                              .view(np.recarray)
        self.pmatrix.shape = (self.ndetectors, self.nsamples,
                              self.npixels_per_sample)

        shapein = tuple([self.header['NAXIS'+str(i+1)] for i in \
                         range(self.header['NAXIS'])])[::-1]
        mask = Map.empty(shapein, dtype=np.bool8, header=self.header)
        tmf.pointing_matrix_mask(self._pmatrix, mask.view(np.int8).T, 
            self.npixels_per_sample, self.nsamples, self.ndetectors)

        ismapdistributed = self.comm_map.Get_size() > 1
        istoddistributed = self.comm_tod.Get_size() > 1
        self.ispacked = packed or ismapdistributed
        if self.ispacked:
            tmf.pointing_matrix_pack(self._pmatrix, mask.view(np.int8).T,
                self.npixels_per_sample, self.nsamples, self.ndetectors)
            shapein = (int(np.sum(~mask)))

        if self.unitin:
            self.unitin = 'Jy ' + self.unitin
        if self.unitout:
            self.unitout = 'Jy ' + self.unitout

        shapeout = (self.ndetectors, np.sum(nsamples))
        Operator.__init__(self, shapein=shapein, shapeout=shapeout,
                          dtype=var.FLOAT_DTYPE)
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
        self.nsamples = s.nsamples
        self.pmatrix = s.pmatrix

    def direct(self, input, output):
#        output.__class__ = Tod
#        output.__array_finalize__(None)
#        output.unit = self.unitout
#        output.derived_units = self.duout
        output.header = None
        tmf.pointing_matrix_direct(self._pmatrix, input.T, output.T,
                                   self.npixels_per_sample)

    def transpose(self, input, output):
#        output.__class__ = Map
#        output.__array_finalize__(None)
#        output.header = self.header
#        output.unit = self.unitin
#        output.derived_units = self.duin
        output[...] = 0
        tmf.pointing_matrix_transpose(self._pmatrix, input.T, output.T, 
                                      self.npixels_per_sample)

    def get_ptp(self):
        npixels = np.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample,
            self.nsamples, self.ndetectors, npixels).T


@real
@linear
class Projection(Operator):
    """
    This class handles operations by the pointing matrix

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

    def __new__(cls, input, method=None, header=None, resolution=None,
                npixels_per_sample=0, oversampling=True, **keywords):
        if not isinstance(input, Observation) or len(input.slice) == 1:
            instance = super(Projection, cls).__new__(cls)
            return instance
        obs = input
        if header is None:
            header = obs.get_map_header(resolution=resolution, oversampling=oversampling)
        elif isinstance(header, str):
            header = str2fitsheader(header)
        operands = []
        for islice in range(len(obs.slice)):
            pmatrix, method, ndetectors, nsamples, npixels_per_sample, units, \
                derived_units = obs.get_pointing_matrix(header,
                npixels_per_sample, method=method, oversampling=oversampling,
                islice=islice)
            p = Projection(pmatrix, method=method, header=header, ndetectors=ndetectors, nsamples=nsamples,
                           npixels_per_sample=npixels_per_sample, units=units, derived_units=derived_units, comm_tod=obs.comm_tod)
            operands.append(p)
        return ExpansionOperator(operands, partitionout=obs.get_nsamples(),
                                 axisout=-1)
        
    def __init__(self, input, method=None, header=None, ndetectors=None,
                 nsamples=None, npixels_per_sample=0, units=None,
                 derived_units=None, resolution=None, oversampling=True,
                 comm_map=None, comm_tod=None, packed=False):

        self.comm_map = comm_map or var.comm_map

        if isinstance(input, Observation):
            if header is None:
                header = input.get_map_header(resolution=resolution,
                                              oversampling=oversampling)
            elif isinstance(header, str):
                header = str2fitsheader(header)

            comm_tod = input.comm_tod
            pmatrix, method, ndetectors, nsamples, npixels_per_sample, units, \
                derived_units = input.get_pointing_matrix(
                    header,
                    npixels_per_sample,
                    method=method,
                    oversampling=oversampling)
        else:
            pmatrix = input

        self.comm_tod = comm_tod
        self.header = header
        self.ndetectors = ndetectors
        self.nsamples = nsamples
        self.npixels_per_sample = npixels_per_sample
        self.unitout, self.unitin = units
        self.duout, self.duin = derived_units
        self._pmatrix = pmatrix
        if self.npixels_per_sample == 0:
            pmatrix = np.empty(0, dtype=np.int64)
        self.pmatrix = pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]) \
                              .view(np.recarray)
        self.pmatrix.shape = (self.ndetectors, self.nsamples,
                              self.npixels_per_sample)

        shapein = tuple([self.header['NAXIS'+str(i+1)] for i in \
                         range(self.header['NAXIS'])])[::-1]
        mask = Map.empty(shapein, dtype=np.bool8, header=self.header)
        tmf.pointing_matrix_mask(self._pmatrix, mask.view(np.int8).T, 
            self.npixels_per_sample, self.nsamples, self.ndetectors)

        ismapdistributed = self.comm_map.Get_size() > 1
        istoddistributed = self.comm_tod.Get_size() > 1
        self.ispacked = packed or ismapdistributed
        if self.ispacked:
            tmf.pointing_matrix_pack(self._pmatrix, mask.view(np.int8).T,
                self.npixels_per_sample, self.nsamples, self.ndetectors)
            shapein = (int(np.sum(~mask)))

        if self.unitin:
            self.unitin = 'Jy ' + self.unitin
        if self.unitout:
            self.unitout = 'Jy ' + self.unitout

        shapeout = (self.ndetectors, self.nsamples)
        Operator.__init__(self, shapein=shapein, shapeout=shapeout,
                          dtype=var.FLOAT_DTYPE)
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
        self.pmatrix = s.pmatrix

    def direct(self, input, output):
        #output.unit = self.unitout
        #output.derived_units = self.duout
        #output.header = None
        if not output.flags.contiguous:
            print 'Warning, Projection output is not contiguous.'
            output_ = np.empty_like(output)
        else:
            output_ = output
        tmf.pointing_matrix_direct(self._pmatrix, input.T, output_.T,
                                   self.npixels_per_sample)
        if not output.flags.contiguous:
            output[...] = output_
        #output.__class__ = Tod


    def transpose(self, input, output, operation=None):
        #output.header = self.header
        #output.unit = self.unitin
        #output.derived_units = self.duin
        if not input.flags.contiguous:
            print 'Warning, Projection input is not contiguous.'
            input_ = np.ascontiguousarray(input)
        else:
            input_ = input
        if operation is None:
            output[...] = 0
        elif operation is not operator.__iadd__:
            raise ValueError('Invalid operation.')
        tmf.pointing_matrix_transpose(self._pmatrix, input_.T, output.T, 
                                      self.npixels_per_sample)
        #output.__class__ = Map

    def get_ptp(self):
        npixels = np.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample,
            self.nsamples_tot, self.ndetectors, npixels).T


@real
@linear
@square
@inplace
class ResponseTruncatedExponential(Operator):
    """
    ResponseTruncatedExponential(tau)

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
class DiscreteDifference(Operator):
    """
    Discrete difference operator.

    Calculate the nth order discrete difference along given axis.
    """

    def __init__(self, axis=-1, comm=None, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        self.axis = axis
        self.comm = comm or var.comm_map
        self.add_rule('.T.', self._rule_ddtdd)

    def direct(self, input, output):
        diff(input, output, self.axis, comm=self.comm)

    def transpose(self, input, output):
        diffT(input, output, self.axis, comm=self.comm)

    def _rule_ddtdd(self, dT):
        return DdTdd(self.axis, comm=self.comm)


@real
@symmetric
@inplace
class DdTdd(Operator):
    """Calculate operator dX.T dX along a given axis."""

    def __init__(self, axis=-1, scalar=1., comm=None, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        self.__name__ = (str(scalar) + ' ' if scalar != 1 else '') + \
            self.__name__
        self.axis = axis
        self.scalar = scalar
        self.comm = comm or var.comm_map        
        self.add_rule('{ScalarOperator}.',lambda s: DdTdd(self.axis,
            s.data * self.scalar))

    def direct(self, input, output):
        diffTdiff(input, output, self.axis, self.scalar, comm=self.comm)


@real
@linear
class DistributionGlobal(Operator):
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
            def direct(input, output):
                if input.data is output.data:
                    pass
                output[:] = input
            def transpose(input, output):
                if input.data is not output.data:
                    output[:] = input
                if self.comm.Get_size() > 1:
                    self.comm.Allreduce(MPI.IN_PLACE, [output, MPI.DOUBLE],
                                        op=MPI.SUM)
            Operator.__init__(self, shapein=shape, classin=Map,
                                    direct=direct, transpose=transpose,
                                    **keywords)
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
        Operator.__init__(self, shapein=shape, shapeout=shapeout,
                                attrin=attrin, attrout=attrout, classin=Map,
                                **keywords)

    def direct(self, input, output):
        s = split_work(self.shapein[0], comm=self.comm)
        n = s.stop - s.start
        output[0:n] = input[s.start:s.stop]
        output[n:] = 0

    def transpose(self, input, output):
        s = split_work(self.shapein[0], comm=self.comm)
        n = s.stop - s.start
        self.comm.Allgatherv([input[0:n], MPI.DOUBLE], [output, (self.counts,
                             self.offsets), MPI.DOUBLE])


@real
@linear
class DistributionLocal(Operator):
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
        Operator.__init__(self, classin=Map, shapein=shapein,
                                shapeout=shapeout, attrin=attrin,
                                attrout=attrout, **keywords)
        self.comm = comm
        self.mask = mask
        self.operator = operator

    def direct(self, input, output):
        status = tmf.mpi_allscatterlocal(input.T, self.mask.view(np.int8).T,
            output.T, self.comm.py2f())
        if status == 0: return
        if status < 0:
            raise RuntimeError('Incompatible sizes.')
        raise MPI.Exception(status)

    def transpose(self, input, output):
        status = tmf.mpi_allreducelocal(input.T, self.mask.view(np.int8).T,
            output.T, self.operator.py2f(), self.comm.py2f())
        if status == 0: return
        if status < 0:
            raise RuntimeError('Incompatible mask.')
        raise MPI.Exception(status)


@real
@symmetric
@idempotent
@inplace
class Masking(Operator):
    """
    Mask operator.

    Sets to zero values whose mask is True (non-null). The input of a Masking
    instance can be of rank greater than the speficied mask, in which case the
    latter is broadcast along the fast dimensions.
    """

    def __init__(self, mask, **keywords):
        if mask is None:
            print('Warning: input mask is None.')
            mask = False
        mask = np.array(mask, order='c', dtype=np.bool8)
        self.isscalar = mask.ndim == 0
        self.mask = np.array(mask, ndmin=1, copy=False)
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)

    def direct(self, input, output):
        if self.same_data(input, output):
            tmf.masking_inplace(input.ravel(), self.mask.view(np.int8).ravel())
        else:
            tmf.masking_outplace(input.ravel(), self.mask.view(np.int8).ravel(),
                                 output.ravel())

    def validate_shapein(self, shapein):
        if shapein is None:
            return self.shapein
        if self.isscalar:
            return shapein
        if shapein[0:self.mask.ndim] != self.mask.shape:
            raise ValueError('The input has a shape ' + str(shapein) + ' incomp'
                'atible with that of the mask ' + str(self.mask.shape) + '.')
        return shapein


@real
@linear
@inplace
class Packing(Operator):
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
class Unpacking(Operator):
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


class Slicing(Operator):
    "Pads before and after a Tod."

    def __init__(self, shapein, slices, **keywords):
        if isinstance(slices, slice):
            slices = (slice,)

        if len(slices) != len(self.shapein):
            raise ValueError('The dimensions of the slices is incompatible wi'\
                             'th the input shape.')

        if shapein is not None:
            shapeout = self.validate_shapein(shapein)
        else:
            shapeout = None
        Operator.__init__(self, shapein=shapein, shapeout=shapeout,
                                **keywords)

        self.slices = slices
        self.mask = np.ones(shapeout, dtype=np.bool8)
        self.mask[self.slices] = False

    def direct(self, input, output):
        if input.data is output.data or output.base == id(input):
            return
        output[:] = input[self.slices]

    def transpose(self, input, output):
        if input.data is output.data:
            return
        if input.base != id(output):
            output[self.slices] = input
        tmf.masking_inplace(output.ravel(), self.mask.view(np.int8).ravel())


@real
@linear
@square
@inplace
class ShiftOperator(Operator):

    def __init__(self, n, axis=-1, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
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
    
    def __init__(self, n, axis=-1, **keywords):
        Operator.__init__(self, **keywords)
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
class Fft(Operator):
    """
    Performs complex fft
    """

    def __init__(self, shape, flags=['estimate'], **keywords):
        Operator.__init__(self, shapein=shape,
                                dtype=var.COMPLEX_DTYPE, **keywords)
        if fftw3.planning.lib_threads is None:
            nthreads = 1
        else:
            nthreads = openmp_num_threads()
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
class FftHalfComplex(Operator):
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

    def reshapein(self, shapein):
        if shapein is None:
            return None
        if shapein[-1] != self.size:
            raise ValueError("Invalid input dimension '{0}'. Expected dimension"
                             " is '{1}'".format(shapein[-1], self.size))
        return shapein

@real
@linear
@inplace
class Convolution(Operator):
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

        nthreads = nthreads or openmp_num_threads()

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
    array_ = array.view(np.ndarray).reshape((-1, array.shape[-1]))
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
    
