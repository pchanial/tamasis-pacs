from __future__ import division

import functools
import numpy as np

import tamasisfortran as tmf

try:
    import fftw3
except:
    print('Warning: Library PyFFTW3 is not installed.')

from mpi4py import MPI

from operators import (Operator, IdentityOperator, DiagonalOperator,
                       PartitionOperator)
from operators.decorators import (idempotent, linear, orthogonal, real, square,
                                  symmetric, unitary)
from operators.utils import isscalar, tointtuple, openmp_num_threads

from . import var

from .datatypes import Map, Tod
from .utils import diff, diffT, diffTdiff, shift
from .mpiutils import split_shape, split_work

__all__ = [
    'DistributionGlobal',
    'DistributionLocal',
    'RollOperator',
    'Convolution',
    'Diagonal',
    'DdTdd',
    'DiscreteDifference',
    'Fft',
    'FftHalfComplex',
    'InterpolationLinear',
    'ShiftOperator',
    'CompressionAverage',
    'DownSampling',
    'CompressionAverage_old',
    'InvNtt',
    'Padding',
    'Projection',
    'ResponseTruncatedExponential',
    'SqrtInvNtt',
    'DistributionGlobal',
    'DistributionLocal',
    'Convolution',
    'DdTdd',
    'DiscreteDifference',
    'Fft',
    'FftHalfComplex',
    'FftHalfComplex_old',
    'InterpolationLinear',
    'Masking',
]

Diagonal = DiagonalOperator

def partitioned_last_axis(cls):

    @staticmethod
    @functools.wraps(cls.__new__)
    def partitioned_new(cls_, *args, **keywords):
        if 'nsamples' in keywords:
            nsamples = tointtuple(keywords['nsamples'])
            del keywords['nsamples']
        else:
            nsamples = None
            print 'k', keywords
        if nsamples is None or len(nsamples) == 1:
            inst = cls.__new_original__(cls_, *args, **keywords)
            print inst.__class__ 
            cls_.__init__ = partitioned_init
            return inst

        n = len(nsamples)
        argss = tuple(tuple(a if isscalar(a) else a[i] for a in args) \
                      for i in range(n))
        keyss = tuple(dict((k, v if isscalar(v) else v[i]) \
                      for k,v in keywords.iteritems()) for i in range(n))
        operands = [cls.__new_original__(cls_, *a, **k) \
                    for a,k in zip(argss, keyss)]
        print 'args:', argss
        print 'keys:', keyss
        for o,a,k in zip(operands, argss, keyss):
            o.__init_original__(*a, **k)
        return PartitionOperator(operands, partitionin=nsamples, axisin=-1)

    @functools.wraps(cls.__init__)
    def partitioned_init(self, *args, **keywords):
        try:
            del keywords['nsamples']
        except KeyError:
            pass
        cls.__init_original__(self, *args, **keywords)

    cls.__new_original__ = cls.__new__
    cls.__new__ = partitioned_new
    cls.__init_original__ = cls.__init__
    cls.__init__ = partitioned_init
    return cls

@real
@linear
class Compression_old(Operator):

    """
    Abstract class for compressing the input signal.
    Compression is operated on the fastest axis.
    """

    def __new__(cls, factor, shapein=None, **keywords):
        factor = tointtuple(factor)
        if np.all(factor == 1):
            op = IdentityOperator(shapein=shapein, **keywords)
            op.factor = (1,)
            return op
        return super(Compression_old, cls).__new__(cls, factor, shapein, **keywords)

    def __init__(self, factor, shapein=None, shapeout=None, nsamples=None,
                 **keywords):

        factor = tointtuple(factor)
        nsamples = tointtuple(nsamples)

        if nsamples is None:
            if len(factor) > 1:
                raise ValueError("When dealing with mixed compression factors, "
                "the keyed argument 'nsamples' must be set.")
        elif len(factor) > 1 and len(nsamples) != len(factor):
            raise ValueError("The number of specified compression factors '{0}'"
                " is different from the number of slices given by nsamples '{1}"
                "'".format(len(factor), len(nsamples)))
        elif np.any(np.mod(nsamples, factor) != 0):
            raise ValueError("The compression factor is incompatible with the s"
                "slicing given by nsamples.")

        self.factor = factor

        if shapein is not None:
            shapeout_expected = self.reshapein(shapein)
            if shapeout and tointtuple(shapeout) != shapeout_expected:
                raise ValueError('The output shape is incompatible with the com'
                                 'pression factor.')

        self.nsamplesin = nsamples
        if nsamples is None:
            self.nsamplesout = None
        else:
            self.nsamplesout = tointtuple(np.divide(nsamples, factor))

        Operator.__init__(self, shapein=shapein, dtype=var.FLOAT_DTYPE,
                          **keywords)

    def reshapein(self, shapein):
        if self.nsamplesout is None:
            if shapein[-1] % self.factor != 0:
                raise ValueError("The input timeline size '{0}') is incompatibl"
                    "e with the compression factor '{1}'".format(shapein[-1],
                    self.factor))
            expected = shapein[-1] // self.factor
        else:
            if shapein[-1] != np.sum(self.nsamplesin):
                raise ValueError("The input timeline size '{0}' is incompatible"
                    " with the compression factor '{1}' and the slicing given b"
                    "y nsamples '{2}'".format(shapein[-1], self.factor,
                    self.nsamples))
            expected = int(np.sum(self.nsamplesout))
        shapeout = list(shapein)
        shapeout[-1] = expected
        return tuple(shapeout)

    def reshapeout(self, shapeout):
        if self.nsamplesin is None:
            expected = shapeout[-1] * self.factor[0]
        else:
            expected = int(np.sum(self.nsamplesin))
        shapeout = list(shapeout)
        shapeout[-1] = expected
        return  tuple(shapeout)

    def __str__(self):
        if len(self.factor) == 1:
            sfactor = str(self.factor[0])
        else:
            sfactor = str(self.factor)[1:-1].replace(' ', '')
        return super(Compression_old, self).__str__() + ' (x' + sfactor + ')'


class CompressionAverage_old(Compression_old):
    """
    Compress the input signal by averaging blocks of specified size.
    """

    def direct(self, input, output):
        if self.nsamplesout is None:
            nsamples = (input.shape[-1],)
        else:
            nsamples = self.nsamplesin
            output.nsamples = self.nsamplesout
        tmf.compression_average_direct(input.T, output.T,
            np.array(nsamples, np.int32), np.array(self.factor, np.int32))

    def transpose(self, input, output):
        if self.nsamplesin is None:
            nsamples = (input.shape[-1],)
        else:
            nsamples = self.nsamplesout
            output.nsamples = self.nsamplesin
        tmf.compression_average_transpose(input.T, output.T,
            np.array(nsamples, np.int32), np.array(self.factor, np.int32))


@partitioned_last_axis
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

        factor = int(factor)
        if shapein is not None:
            shapeout_expected = self.reshapein(shapein)
            if shapeout and tointtuple(shapeout) != shapeout_expected:
                raise ValueError('The output shape is incompatible with the com'
                                 'pression factor.')

        self.factor = factor
        Operator.__init__(self, shapein=shapein, dtype=var.FLOAT_DTYPE,
                          **keywords)

    def reshapein(self, shapein):
        if shapein[-1] % self.factor != 0:
            raise ValueError("The input timeline size '{0}') is incompatible wi"
                "th the compression factor '{1}'".format(shapein[-1],
                self.factor))
        shapeout = list(shapein)
        shapeout[-1] = shapein[-1] // self.factor
        return tuple(shapeout)

    def reshapeout(self, shapeout):
        shapeout = list(shapeout)
        shapeout[-1] = shapeout[-1] * self.factor[0]
        return  tuple(shapeout)

    def __str__(self):
        return super(Compression, self).__str__() + ' (x{})'.format(self.factor)


@partitioned_last_axis
@real
@linear
class CompressionAverage(Operator):
    """
    Compress the input signal by averaging blocks of specified size.
    """

    def __init__(self, factor, shapein=None, shapeout=None, **keywords):

        factor = int(factor)
        if shapein is not None:
            shapeout_expected = self.reshapein(shapein)
            if shapeout and tointtuple(shapeout) != shapeout_expected:
                raise ValueError('The output shape is incompatible with the com'
                                 'pression factor.')

        self.factor = factor
        Operator.__init__(self, shapein=shapein, dtype=var.FLOAT_DTYPE,
                          **keywords)

    def reshapein(self, shapein):
        if shapein[-1] % self.factor != 0:
            raise ValueError("The input timeline size '{0}') is incompatible wi"
                "th the compression factor '{1}'".format(shapein[-1],
                self.factor))
        shapeout = list(shapein)
        shapeout[-1] = shapein[-1] // self.factor
        return tuple(shapeout)

    def reshapeout(self, shapeout):
        shapeout = list(shapeout)
        shapeout[-1] = shapeout[-1] * self.factor[0]
        return  tuple(shapeout)

    def __str__(self):
        return super(CompressionAverage, self).__str__() + ' (x{})'.format(self.factor)

    def direct(self, input, output):
        nsamples = (input.shape[-1],)
        output_ = np.ascontiguousarray(output)
        tmf.compression_average_direct(input.T, output_.T,
            np.array(nsamples, np.int32), np.array((self.factor,), np.int32))
        output[:] = output_

    def transpose(self, input, output):
        nsamples = (input.shape[-1],)
        output_ = np.ascontiguousarray(output)
        tmf.compression_average_transpose(input.T, output_.T,
            np.array(nsamples, np.int32), np.array(self.factor, np.int32))
        output[:] = output_


class DownSampling(Compression_old):
    """
    Downsample the input signal by picking up one sample out of a number
    specified by the compression factor
    """

    def direct(self, input, output):
        if self.nsamplesout is None:
            nsamples = (input.shape[-1],)
        else:
            nsamples = self.nsamplesin
            output.nsamples = self.nsamplesout
        tmf.downsampling_direct(input.T, output.T,
            np.array(nsamples, np.int32), np.array(self.factor, np.int32))

    def transpose(self, input, output):
        if self.nsamplesin is None:
            nsamples = (input.shape[-1],)
        else:
            nsamples = self.nsamplesout
            output.nsamples = self.nsamplesin
        tmf.downsampling_transpose(input.T, output.T,
            np.array(nsamples, np.int32), np.array(self.factor, np.int32))
      

@real
@symmetric
class InvNtt(Operator):

    def __new__(cls, obs, filename=None, **keywords):
        nsamples = obs.get_nsamples()
        length = np.asarray(2**np.ceil(np.log2(np.array(nsamples) + 200)), int)
        invntt = cls._get_diagonal(length, obs.get_filter_uncorrelated(
                                    filename=filename, **keywords))
        fft = FftHalfComplex(tuple(length), nsamples=tuple(length))
        padding = Padding(left=invntt.ncorrelations, right=length - nsamples - \
                          invntt.ncorrelations, nsamples=nsamples)
        return padding, fft, invntt
        return padding.T * fft.T * invntt * fft * padding

    @staticmethod
    def _get_diagonal(nsamples, filter, **keywords):
        tod_filters = []
        for n in nsamples:
            tod_filter, status = \
                tmf.fft_filter_uncorrelated2(filter.T, n)
            if status != 0: raise RuntimeError()
            np.maximum(tod_filter, 0, tod_filter)
            tod_filters.append(tod_filter)
        norm = var.comm_tod.allreduce(max([np.max(f) for f in tod_filters]),
                                      op=MPI.MAX)
        for f in tod_filters:
            np.divide(f, norm, f)
        d = PartitionOperator([DiagonalOperator(f.T) for f in tod_filters],
                              partitionin=nsamples, axisin=-1)
        d.ncorrelations = filter.shape[-1] - 1
        return d


@real
@symmetric
class SqrtInvNtt(InvNtt):
    def __init__(self, *args, **kw):
        invntt = InvNtt(*args, **kw)
        self._tocomposite(op.Composition, invntt.blocks)
        data = self.blocks[2].data
        np.sqrt(data, data)


@partitioned_last_axis
@real
@linear
class Padding(Operator):
    """Pads before and after along the fast dimension of an ndarray."""

    def __init__(self, left=0, right=0., shapein=None, **keywords):
        left = np.array(left, ndmin=1, dtype=int)
        right = np.array(right, ndmin=1, dtype=int)
        if np.any(left < 0):
            raise ValueError('Left padding is not positive.')
        if np.any(right < 0):
            raise ValueError('Right padding is not positive.')
        if left.ndim != 1 or right.ndim != 1:
            raise ValueError('Padding must be scalar or a vector.')
        if right.size > 1 and left.size > 1 and right.size != left.size:
            raise ValueError('Invalid left and right padding.')

        self.left  = tuple(left)
        if left.size == 1:
            self.left = right.size * self.left
        self.right = tuple(right)
        if right.size == 1:
            self.right = left.size * self.right
        Operator.__init__(self, shapein=shapein, **keywords)
   
    def direct(self, input, output):
        dest = 0
        dest_padded = 0
        for islice in range(1):
            nsamples = input.shape[-1]
            left = self.left[islice if len(self.left) > 1 else 0]
            output[...,dest_padded:dest_padded+left] = 0
            output[...,dest_padded+left:dest_padded+left+nsamples] = \
                input[...,dest:dest+nsamples]
            output[...,dest_padded+left+nsamples:dest_padded+ \
                output.shape[-1]] = 0
            dest += nsamples
            dest_padded += output.shape[-1]
   
    def transpose(self, input, output):
        dest = 0
        dest_padded = 0
        for islice in range(1):
            nsamples = output.shape[-1]
            left  = self.left [islice if len(self.left)  > 1 else 0]
            output[...,dest:dest+nsamples] = \
                input[...,dest_padded+left:dest_padded+left+nsamples]
            dest += nsamples
            dest_padded += input.shape[-1]

    def reshapein(self, shapein):
        shapeout = list(shapein)
        shapeout[-1] += self.left[0] + self.right[0]
        return shapeout
       
    def reshapeout(self, shapeout):
        shapein = list(shapeout)
        shapein[-1] -= self.left[0] + self.right[0]
        return shapein


#-------------------------------------------------------------------------------


@real
@linear
class Projection(Operator):
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
        self.npixels_per_sample, (self.unitout, self.unitin), \
        (self.duout, self.duin) = \
            observation.get_pointing_matrix(header,
                                            resolution,
                                            npixels_per_sample,
                                            method=method,
                                            oversampling=oversampling)

        self.nsamples = nsamples
        self.nsamples_tot = int(np.sum(nsamples))
        self._pmatrix = pmatrix
        if self.npixels_per_sample == 0:
            pmatrix = np.empty(0, dtype=np.int64)
        self.pmatrix = pmatrix.view([('weight', 'f4'), ('pixel', 'i4')]) \
                              .view(np.recarray)
        self.pmatrix.shape = (self.ndetectors, self.nsamples_tot,
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

        if self.unitin:
            self.unitin = 'Jy ' + self.unitin
        if self.unitout:
            self.unitout = 'Jy ' + self.unitout

        shapeout = (self.ndetectors, np.sum(nsamples))
        Operator.__init__(self, shapein=shapein, shapeout=shapeout)
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

    def direct(self, input, output):
        output.__class__ = Tod
        output.__array_finalize__(None)
        output.unit = self.unitout
        output.derived_units = self.duout
        output.header = None
        output.nsamples = self.nsamples
        tmf.pointing_matrix_direct(self._pmatrix, input.T, output.T,
                                   self.npixels_per_sample)

    def transpose(self, input, output):
        output.__class__ = Map
        output.__array_finalize__(None)
        output.header = self.header
        output.unit = self.unitin
        output.derived_units = self.duin
        if hasattr(output, 'nsamples'):
            delattr(output, 'nsamples')
        tmf.pointing_matrix_transpose(self._pmatrix, input.T, output.T, 
                                      self.npixels_per_sample)

    def get_ptp(self):
        npixels = np.product(self.shapein)
        return tmf.pointing_matrix_ptp(self._pmatrix, self.npixels_per_sample,
            self.nsamples_tot, self.ndetectors, npixels).T


@real
@linear
@square
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
        Operator.__init__(self, classin=Tod, **keywords)

    def direct(self, input, output):
        if self.nsamples is not None:
            nsamples = self.nsamples
            output.nsamples = nsamples
        else:
            nsamples = (input.shape[-1],)
        tmf.convolution_trexp_direct(input.ravel(), output.ravel(),
                                     np.array(nsamples), self.tau)

    def transpose(self, input, output):
        if self.nsamples is not None:
            nsamples = self.nsamples
            output.nsamples = nsamples
        else:
            nsamples = (input.shape[-1],)
        tmf.convolution_trexp_transpose(input.ravel(), output.ravel(),
                                        np.array(nsamples), self.tau)


@real
@linear
@square
class DiscreteDifference(Operator):
    """
    Discrete difference operator.

    Calculate the nth order discrete difference along given axis.
    """

    def __init__(self, axis=-1, comm=None, **keywords):
        Operator.__init__(self, **keywords)
        self.axis = axis
        self.comm = comm or var.comm_map

    def direct(self, input, output):
        diff(input, output, self.axis, comm=self.comm)

    def transpose(self, input, output):
        diffT(input, output, self.axis, comm=self.comm)


@real
@symmetric
class DdTdd(Operator):
    """Calculate operator dX.T dX along a given axis."""

    def __init__(self, axis=-1, scalar=1., description=None, comm=None,
                 **keywords):
        if description is None and scalar != 1.:
            description = str(scalar) + ' DdTdd'
        Operator.__init__(self, **keywords)
        self.axis = axis
        self.scalar = scalar
        self.comm = comm or var.comm_map

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
        Operator.__init__(self, **keywords)

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
class Packing(Operator):
    """
    Convert an nd array in a 1d map, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, **keywords):
        mask = np.array(mask, np.bool8)
        Operator.__init__(self, shapein=mask.shape, shapeout=np.sum(mask == 0),
                          **keywords)
        self.mask = mask

    def direct(self, input, output):
        tmf.packing(input.ravel(), self.mask.view(np.int8).ravel(),
                    output.ravel())

    def transpose(self, input, output):
        tmf.unpacking(input.ravel(), self.mask.view(np.int8).ravel(),
                      output.ravel(), 0.)


@real
@linear
class Unpacking(Operator):
    """
    Convert 1d map into an nd array, under the control of a mask.
    The elements for which the mask is True are equal to the field argument.
    """

    def __init__(self, mask, **keywords):
        mask = np.array(mask, np.bool8)
        Operator.__init__(self, shapein=np.sum(mask == 0), shapeout=mask.shape,
                          **keywords)
        self.mask = mask

    def direct(self, input, output):
        tmf.unpacking(input.ravel(), self.mask.view(np.int8).ravel(),
                      output.ravel(), 0.)

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
class ShiftOperator(Operator):

    def __init__(self, n, axis=-1, **keywords):
        Operator.__init__(self, **keywords)
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
class RollOperator(Operator):
    
    def __init__(self, n, axis=-1, **keywords):
        Operator.__init__(self, **keywords)
        self.axis = (axis,) if isscalar(axis) else tuple(axis)
        self.n = (n,) * len(self.axis) if isscalar(n) else tuple(n)
        if len(self.axis) != len(self.n):
            raise ValueError('There is a mismatch between the number of axes an'
                             'd offsets.')

    def direct(self, input, output):
        output[:] = np.roll(input, self.n[0], axis=self.axis[0])
        for n, axis in zip(self.n, self.axis)[1:]:
            output[:] = np.roll(output, n, axis=axis)

    def transpose(self, input, output):
        output[:] = np.roll(input, -self.n[0], axis=self.axis[0])
        for n, axis in zip(self.n, self.axis)[1:]:
            output[:] = np.roll(output, -n, axis=axis)


@unitary
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
        self._in[:] = input
        fftw3.execute(self.forward_plan)
        output[:] = self._out

    def transpose(self, input, output):
        self._in[:] = input
        fftw3.execute(self.backward_plan)
        output[:] = self._out / self.n


@real
@orthogonal
class FftHalfComplex_old(Operator):
    """
    Performs real-to-half-complex fft
    """

    def __init__(self, nsamples, **keywords):
        Operator.__init__(self, **keywords)
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

    def direct(self, input, output):
        print self.same_data(input, output)
        print input.__array_interface__
        print output.__array_interface__
        if not self.same_data(input, output):
            output[:] = input
        output_ = np.reshape(output, (-1, input.shape[-1]))
        tmf.fft_plan(output_.T, np.array(self.nsamples), self.forward_plan)

    def transpose(self, input, output):
        if not self.same_data(input, output):
            output[:] = input
        output_ = np.reshape(output, (-1, input.shape[-1]))
        tmf.fft_plan(output_.T, np.array(self.nsamples), self.backward_plan)
        dest = 0
        for n in self.nsamples:
            output_[:,dest:dest+n] /= n
            dest += n

    def validate_shapein(self, shape):
        if shape is None:
            return None
        nsamples = shape[-1]
        if shape[-1] != sum(self.nsamples):
            raise ValueError("Invalid FFT size '" + str(nsamples) + \
                                  "' instead of '"+str(sum(self.nsamples))+"'.")
        return shape


@partitioned_last_axis
@real
@orthogonal
class FftHalfComplex(Operator):
    """
    Performs real-to-half-complex fft
    """

    def __init__(self, size, **keywords):
        Operator.__init__(self, **keywords)
        self.nsamples = tointtuple(size)
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

    def direct(self, input, output):
        output_ = np.ascontiguousarray(input).reshape((-1,input.shape[-1]))
        tmf.fft_plan(output_.T, np.array(self.nsamples), self.forward_plan)
        output[:] = output_

    def transpose(self, input, output):
        output_ = np.ascontiguousarray(input).reshape((-1,input.shape[-1]))
        tmf.fft_plan(output_.T, np.array(self.nsamples), self.backward_plan)
        dest = 0
        for n in self.nsamples:
            output_[:,dest:dest+n] /= n
            dest += n
        output[:] = output_

    def _reshapein(self, shape):
        if shape is None:
            return None
        nsamples = shape[-1]
        if shape[-1] != sum(self.nsamples):
            raise ValueError("Invalid FFT size '" + str(nsamples) + \
                             "' instead of '"+str(sum(self.nsamples))+"'.")
        return shape

    def __str__(self):
        return self.__class__.__name__ + '(size={})'.format(self.nsamples[0])


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
        self._in[:] = input
        fftw3.execute(self.forward_plan)
        self._out *= self.kernel
        fftw3.execute(self.backward_plan)
        output[:] = self._in
        
    def transpose(self, input, output):
        self._in[:] = input
        fftw3.execute(self.forward_plan)
        # avoid temp array in 'self._out *= self.kernel.conjugate()'
        opi.multiply_conjugate_complex(self._out.ravel(), self.kernel.ravel(),
                                       self._out.ravel())
        fftw3.execute(self.backward_plan)
        output[:] = self._in


class InterpolationLinear(Operator):

    def __init__(self, mask, **keywords):
        Operator.__init__(self, attrin={'mask':mask}, **keywords)
        self.mask = mask

    def direct(self, input, output):
        output[:] = interpolate_linear(output)

    def transpose(self, input, output):
        raise NotImplementedError()

@partitioned_last_axis
class myOp(Operator):
    """ Bla1. """
    def __init__(self, value, **keywords):
        """ Bla2. """
        self.value=value
        Operator.__init__(self, **keywords)
    def direct(self, input, output):
        output[:] = self.value * input

    def __str__(self):
        return self.__class__.__name__ + '(value={})'.format(self.value)
