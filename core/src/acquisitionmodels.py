from __future__ import division

import functools
import numpy as np

from pyoperators import (BlockColumnOperator, BlockDiagonalOperator,
                         CompositionOperator, DiagonalOperator,
                         MPIDistributionIdentityOperator, IdentityOperator,
                         Operator, MaskOperator, ZeroOperator)
from pyoperators.decorators import (linear, orthogonal, real, square, symmetric,
                                    contiguous, inplace, separable)
from pyoperators.utils import (isalias, isscalar, openmp_num_threads,
                               product)
from pyoperators.utils.mpi import MPI, distribute_shape, distribute_shapes
from pysimulators import Map
from pysimulators.acquisitionmodels import (block_diagonal, PointingMatrix,
         ProjectionInMemoryOperator, ProjectionOnFlyOperator)
from pysimulators.wcsutils import fitsheader2shape, str2fitsheader

from . import tamasisfortran as tmf
from . import var
from .utils import diff, diffT, diffTdiff, shift
from .mpiutils import gather_fitsheader_if_needed, scatter_fitsheader

try:
    import fftw3
    MAX_FFTW_NUM_THREADS = 1 if fftw3.planning.lib_threads is None \
        else openmp_num_threads()
except:
    print('Warning: Library PyFFTW3 is not installed.')

__all__ = [
    'CompressionAverage', # obsolete
    'CompressionAverageOperator',
    'MPIDistributionLocalOperator',
    'DdTdd', # Obsolete
    'DdTddOperator',
    'DiscreteDifference', # obsolete
    'DiscreteDifferenceOperator',
    'DownSamplingOperator',
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
    'ShiftOperator',
    'SqrtInvNttOperator',
    'Unpacking', # obsolete
    'UnpackOperator',
]


@block_diagonal('factor', axisin=-1)
@real
@linear
@separable
class CompressionOperator(Operator):
    """
    Abstract class for compressing the input signal.
    Compression is operated on the fastest axis.

    """
    def __init__(self, factor, **keywords):
        if factor == 1:
            self.__class__ = IdentityOperator
            self.__init__(**keywords)
            return
        self.factor = int(factor)
        keywords['dtype'] = float
        Operator.__init__(self, **keywords)

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


@contiguous
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


@contiguous
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

    def __new__(cls, obs=None, method='uncorrelated', **keywords):
        if obs is None:
            return Operator.__new__(cls)
        nsamples = obs.get_nsamples()
        comm_tod = obs.instrument.comm
        method = method.lower()
        if method not in ('uncorrelated', 'uncorrelated python'):
            raise ValueError("Invalid method '{0}'.".format(method))

        filter = obs.get_filter_uncorrelated(**keywords)
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
        for i, n in enumerate(filter_length):
            i = min(i, nfilters-1)
            fft_filter, status = tmf.fft_filter_uncorrelated(filter[i].T, n)
            if status != 0: raise RuntimeError()
            np.maximum(fft_filter, 0, fft_filter)
            fft_filters.append(fft_filter.T)
        
        norm = comm_tod.allreduce(max(np.max(f) for f in fft_filters),
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
        self.left = ncorrelations
        self.right = filter_length - nsamples - ncorrelations
        self.fftw_flags = fftw_flags
        self.fplan = fftw3.Plan(array, direction='forward', flags=fftw_flags,
            realtypes=['halfcomplex r2c'], nthreads=1)
        self.bplan = fftw3.Plan(array, direction='backward', flags=fftw_flags,
            realtypes=['halfcomplex c2r'], nthreads=1)
        Operator.__init__(self, shapein=(fft_filter.shape[0], nsamples))
        
    def direct(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        if isalias(input, output):
            tmf.operators.invntt_uncorrelated_inplace(input_, ishape[1],
                istride, self.fft_filter.T, self.fplan._get_parameter(),
                self.bplan._get_parameter(), self.left, self.right)
            return
        output_, oshape, ostride = _ravel_strided(output)
        tmf.operators.invntt_uncorrelated_outplace(input_, ishape[1], istride,
            output_, ostride, self.fft_filter.T, self.fplan._get_parameter(),
            self.bplan._get_parameter(), self.left, self.right)


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
@separable
class PadOperator(Operator):
    """
    Pads before and after along the fast dimension of an ndarray.

    """
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
        right = -self.right if self.right != 0 else output.shape[-1]
        output[...,:self.left] = 0
        output[...,self.left:right] = input
        output[...,right:] = 0
   
    def transpose(self, input, output):
        right = -self.right if self.right != 0 else input.shape[-1]
        output[...] = input[...,self.left:right]

    def reshapein(self, shapein):
        return shapein[:-1] + (shapein[-1] + self.left + self.right,)
       
    def reshapeout(self, shapeout):
        return shapeout[:-1] + (shapeout[-1] - self.left - self.right,)


def ProjectionOperator(input, method=None, header=None, resolution=None,
                       npixels_per_sample=0, units=None, derived_units=None,
                       downsampling=False, packed=False, commin=MPI.COMM_WORLD,
                       commout=MPI.COMM_WORLD, onfly_func=None, onfly_ids=None,
                       onfly_shapeouts=None, **keywords):
    """
    Projection operator factory, to handle operations by one or more pointing
    matrices.

    It is assumed that each pointing matrix row has at most 'npixels_per_sample'
    non-null elements, or in other words, an output element can only intersect
    a fixed number of input elements.

    Given an input vector x and an output vector y, y = P(x) translates into:
        y[i] = sum(P.matrix[i,:].value * x[P.matrix[i,:].index])

    If the input is not MPI-distributed unlike the output, the projection
    operator is automatically multiplied by the operator MPIDistributionIdenti-
    tyOperator, to enable MPI reductions.

    If the input is MPI-distributed, this operator is automatically packed (see
    below) and multiplied by the operator MPIDistributionLocalOperator, which
    takes the local input as argument.

    Arguments
    ---------
    input : pointing matrix (or sequence of) or observation (deprecated)
    method : deprecated
    header : deprecated
    resolution : deprecated
    npixels_per_sample : deprecated
    downsampling : deprecated
    packed deprecated

    """
    # check if there is only one pointing matrix
    isonfly = input is None
    isobservation = hasattr(input, 'get_pointing_matrix')
    if isobservation:
        if hasattr(input, 'slice'):
            nmatrices = len(input.slice)
        else:
            nmatrices = 1
        commout = input.instrument.comm
    elif isonfly:
        nmatrices = len(onfly_ids)
    else:
        if isinstance(input, PointingMatrix):
            input = (input,)
        if any(not isinstance(i, PointingMatrix) for i in input):
            raise TypeError('The input is not a PointingMatrix, (nor a sequence'
                            ' of).')
        nmatrices = len(input)

    ismapdistributed = commin.size > 1
    istoddistributed = commout.size > 1

    # get the pointing matrix from the input observation
    if isobservation:
        if header is None:
            if not hasattr(input, 'get_map_header'):
                raise AttributeError("No map header has been specified and "
                    "the observation has no 'get_map_header' method.")
            header_global = input.get_map_header(resolution=resolution,
                                                 downsampling=downsampling)
        else:
            if isinstance(header, str):
                header = str2fitsheader(header)
            header_global = gather_fitsheader_if_needed(header, comm=commin)
        #XXX we should hand over the local header
        input = input.get_pointing_matrix(header_global, npixels_per_sample,
                    method=method, downsampling=downsampling, comm=commin)
        if isinstance(input, PointingMatrix):
            input = (input,)
    elif isonfly:
        header_global = header
    else:
        header_global = input[0].info['header']

    # check shapein
    if not isonfly:
        shapeins = [i.shape_input for i in input]
        if any(s != shapeins[0] for s in shapeins):
            raise ValueError('The pointing matrices do not have the same input '
                             "shape: {0}.".format(', '.join(str(shapeins))))
        shapein = shapeins[0]
    else:
        shapein = fitsheader2shape(header_global)

    # the output is simply a ProjectionOperator instance
    if nmatrices == 1 and not ismapdistributed and not istoddistributed \
       and not packed:
        return ProjectionInMemoryOperator(input[0], units=units, derived_units=
            derived_units, commin=commin, commout=commout, **keywords)

    if packed or ismapdistributed:
        if isonfly:
            raise NotImplementedError()
        # compute the map mask before this information is lost while packing
        mask_global = Map.ones(shapein, dtype=np.bool8, header=header_global)
        for i in input:
            i.get_mask(out=mask_global)
        for i in input:
            i.pack(mask_global)
        if ismapdistributed:
            shapes = distribute_shapes(mask_global.shape, comm=commin)
            shape = shapes[commin.rank]
            mask = np.empty(shape, bool)
            commin.Reduce_scatter([mask_global, MPI.BYTE], [mask, MPI.BYTE],
                                  [product(s) for s in shapes], op=MPI.BAND)
        else:
            mask = mask_global

    if isonfly:
        place_holder = {}
        operands = [ProjectionOnFlyOperator(place_holder, id, onfly_func,
                                            shapein=shapein, shapeout=shapeout,
                                            units=units,
                                            derived_units=derived_units)
                    for id, shapeout in zip(onfly_ids, onfly_shapeouts)]
    else:
        operands = [ProjectionInMemoryOperator(i, units=units, derived_units=
                    derived_units, commout=commout, **keywords) for i in input]
    result = BlockColumnOperator(operands, axisout=-1)

    if nmatrices > 1:
        def apply_mask(self, mask):
            mask = np.asarray(mask, np.bool8)
            dest = 0
            if mask.shape != self.shapeout:
                raise ValueError("The mask shape '{0}' is incompatible with tha"
                    "t of the projection operator '{1}'.".format(mask.shape,
                    self.shapeout))
            if any(isinstance(p, ProjectionOnFlyOperator)
                   for p in self.operands):
                blocks = self.copy()
                self.__class__ = CompositionOperator
                self.__init__([MaskOperator(mask), blocks])
                return
            for p in self.operands:
                n = p.matrix.shape[1]
                p.apply_mask(mask[...,dest:dest+n])
                dest += n
        def get_mask(self, out=None):
            for p in self.operands:
                out = p.get_mask(out=out)
            return out
        def get_pTp(self, out=None):
            for p in self.operands:
                out = p.get_pTp(out=out)
            return out
        def get_pTx_pT1(self, x, out=None, mask=None):
            dest = 0
            for p in self.operands:
                n = p.matrix.shape[1]
                out = p.get_pTx_pT1(x[...,dest:dest+n], out=out, mask=mask)
                dest += n
            return out
        def intersects(self, out=None):
            raise NotImplementedError('email-me')

        result.apply_mask = apply_mask.__get__(result)
        result.get_mask = get_mask.__get__(result)
        result.get_pTp = get_pTp.__get__(result)
        result.get_pTx_pT1 = get_pTx_pT1.__get__(result)
        result.intersects = intersects.__get__(result)

    if not istoddistributed and not ismapdistributed and not packed:
        return result

    if packed or ismapdistributed:
        def get_mask(self, out=None):
            if out is not None:
                out &= mask
            else:
                out = mask
            return out
        if ismapdistributed:
            header = scatter_fitsheader(header_global, comm=commin)
            result *= MPIDistributionLocalOperator(mask_global, commin=commin,
                                                   attrin={'header':header})
        elif packed:
            result *= PackOperator(mask)
        result.get_mask = get_mask.__get__(result)

    if istoddistributed and not ismapdistributed:
        def get_mask(self, out=None):
            out = self.operands[0].get_mask(out=out)
            commout.Allreduce(MPI.IN_PLACE, [out, MPI.BYTE], op=MPI.BAND)
            return out        
        result *= MPIDistributionIdentityOperator(commout=commout)
        result.get_mask = get_mask.__get__(result)

    def not_implemented(out=None):
        raise NotImplementedError('email-me')

    def apply_mask(self, mask):
        self.operands[0].apply_mask(mask)
    result.apply_mask = apply_mask.__get__(result)
    result.get_pTp = not_implemented
    result.get_pTx_pT1 = not_implemented
    result.intersects = not_implemented

    return result


@real
@linear
@square
@inplace
@separable
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
@contiguous
class DiscreteDifferenceOperator(Operator):
    """
    Discrete difference operator.

    Calculate the nth order discrete difference along given axis.

    """
    def __init__(self, axis=-1, **keywords):
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        self.axis = axis
        self.set_rule('.T.', self._rule_ddtdd, CompositionOperator)

    def direct(self, input, output):
        diff(input, output, self.axis, comm=self.commin)

    def transpose(self, input, output):
        diffT(input, output, self.axis, comm=self.commin)

    @staticmethod
    def _rule_ddtdd(dT, self):
        return DdTddOperator(self.axis, commin=self.commin)


@real
@symmetric
@inplace
@contiguous
class DdTddOperator(Operator):
    """
    Calculate operator dX.T dX along a given axis.

    """
    def __init__(self, axis=-1, scalar=1., **keywords):
        if scalar == 0:
            self.__class__ = ZeroOperator
            self.__init__(flags='square', **keywords)
            return
        Operator.__init__(self, dtype=var.FLOAT_DTYPE, **keywords)
        self.axis = axis
        self.scalar = scalar
        self.set_rule('{HomothetyOperator}.', lambda o,s: DdTddOperator(
                      self.axis, o.data * s.scalar), CompositionOperator)

    def direct(self, input, output):
        diffTdiff(input, output, self.axis, self.scalar, comm=self.commin)

    def __str__(self):
        s = (str(self.scalar) + ' ') if self.scalar != 1 else ''
        s += super(DdTddOperator, self).__str__()
        return s


@real
@linear
class MPIDistributionLocalOperator(Operator):
    """
    Scatter a distributed map to different MPI processes under the control of a
    local non-distributed mask.

    """
    def __init__(self, mask, operator=MPI.SUM, commin=MPI.COMM_WORLD,
                 **keywords):
        commout = commin
        shapeout = mask.size - int(np.sum(mask))
        shapein = distribute_shape(mask.shape, comm=commin)
        Operator.__init__(self, shapein=shapein, shapeout=shapeout,
                          commin=commin, commout=commout, **keywords)
        self.mask = mask
        self.chunk_size = product(mask.shape[1:])
        self.operator = operator

    def direct(self, input, output):
        ninput = input.size
        noutput = output.size
        input = input.ravel() if ninput > 0 else np.array(0.)
        output = output.ravel() if noutput > 0 else np.array(0.)
        status = tmf.operators_mpi.allscatterlocal(input, self.mask.view(
            np.int8).T, output, self.chunk_size, self.commin.py2f(),
            ninput=ninput, noutput=noutput)
        if status == 0: return
        if status < 0:
            rank = self.commin.rank
            if status == -1:
                raise ValueError('Invalid mask on node {0}.'.format(rank))
            if status == -2:
                raise ValueError('Invalid input on node {0}.'.format(rank))
            if status == -3:
                raise RuntimeError('Transmission failure {0}.'.format(rank))
            assert False
        raise MPI.Exception(status)

    def transpose(self, input, output):
        output[...] = 0
        status = tmf.operators_mpi.allreducelocal(input.T, self.mask.view(
            np.int8).T, output.T, self.chunk_size, self.operator.py2f(),
            self.commin.py2f())
        if status == 0: return
        if status < 0:
            rank = self.commin.rank
            if status == -1:
                raise ValueError('Invalid mask on node {0}.'.format(rank))
            assert False
        raise MPI.Exception(status)


@real
@linear
@inplace
@contiguous
class PackOperator(Operator):
    """
    Convert an nd array in a 1d map, under the control of a mask.
    The elements for which the mask is true are discarded.

    """
    def __init__(self, mask, **keywords):
        mask = np.array(mask, np.bool8, copy=False)
        Operator.__init__(self, shapein=mask.shape, shapeout=np.sum(mask == 0),
                          dtype=var.FLOAT_DTYPE, **keywords)
        self.mask = mask
        self.set_rule('.T', lambda s: UnpackOperator(s.mask))
        self.set_rule('.{UnpackOperator}', self._rule_left_unpack, CompositionOperator)
        self.set_rule('{UnpackOperator}.', self._rule_right_unpack, CompositionOperator)

    def direct(self, input, output):
        mask = self.mask.view(np.int8).ravel()
        if isalias(input, output):
            tmf.operators.pack_inplace(input.ravel(), mask, input.size)
        else:
            tmf.operators.pack_outplace(input.ravel(), mask, output)

    @staticmethod
    def _rule_left_unpack(self, op):
        if self.mask.shape != op.mask.shape:
            return
        if np.any(self.mask != op.mask):
            return
        return IdentityOperator()

    @staticmethod
    def _rule_right_unpack(op, self):
        if self.mask.shape != op.mask.shape:
            return
        if np.any(self.mask != op.mask):
            return
        return MaskOperator(self.mask)


@real
@linear
@inplace
@contiguous
class UnpackOperator(Operator):
    """
    Convert an nd array in a 1d map, under the control of a mask.
    The elements for which the mask is true are set to zero.

    """
    def __init__(self, mask, **keywords):
        mask = np.array(mask, np.bool8, copy=False)
        Operator.__init__(self, shapein=np.sum(mask == 0), shapeout=mask.shape,
                          dtype=var.FLOAT_DTYPE, **keywords)
        self.mask = mask
        self.set_rule('.T', lambda s: PackOperator(s.mask))

    def direct(self, input, output):
        mask = self.mask.view(np.int8).ravel()
        if isalias(input, output):
            tmf.operators.unpack_inplace(output.ravel(), mask, input.size)
        else:
            tmf.operators.unpack_outplace(input, mask, output.ravel())


@real
@linear
@square
@inplace
@contiguous
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


@real
@orthogonal
@inplace
@contiguous
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
        if isalias(input, output):
            tmf.fft_plan_inplace(input_, ishape[0], ishape[1], istride,
                                 self.ifplan._get_parameter())
            return
        output_, oshape, ostride = _ravel_strided(output)
        tmf.fft_plan_outplace(input_, ishape[0], ishape[1], istride, output_,
                     oshape[1], ostride, self.fplan._get_parameter())

    def transpose(self, input, output):
        input_, ishape, istride = _ravel_strided(input)
        if isalias(input, output):
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
    if type(cls) is not type:
        def func(*args, **kwargs):
            print("The class factory '{0}' has been renamed. Use '{1}' instead."
                  .format(cls.__name__, cls.__name__ + 'Operator'))
            return cls(*args, **kwargs)
        return func

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
def Projection(*args, **kwargs):
    return ProjectionOperator(*args, **kwargs)
@obsolete
class ResponseTruncatedExponential(ConvolutionTruncatedExponentialOperator):
    pass
@obsolete
class Unpacking(UnpackOperator):
    pass
