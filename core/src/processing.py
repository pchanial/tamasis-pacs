# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import numpy as np
import scipy
import time

from pyoperators.utils import strelapsed
from pysimulators.acquisitionmodels import ProjectionBaseOperator
from . import tamasisfortran as tmf

__all__ = [ 'deglitch_l2std',
            'deglitch_l2mad',
            'filter_median',
            'filter_polynomial',
            'filter_nonfinite',
            'interpolate_linear' ]

class ndarraywrap(np.ndarray):
    pass

def _deglitch(tod, projection, nsigma, method):
    time0 = time.time()
    if method not in ('std', 'mad'):
        raise ValueError('Invalid deglitching method.')
    func = {'std':tmf.deglitch_l2b_std, 'mad':tmf.deglitch_l2b_mad}[method]
    if isinstance(projection, ProjectionBaseOperator):
        matrix = projection.matrix
        npixels_per_sample = matrix.shape[-1]
    elif hasattr(projection, 'operands') and \
         all([isinstance(p, ProjectionBaseOperator) for p in projection.operands]):
        ops = projection.operands
        npixels_per_sample = max(p.matrix.shape[-1] for p in ops)
        if tod.shape[-1] != sum(p.matrix.shape[-2] for p in ops):
            raise ValueError('The tod has an incompatible shape.')
        matrix = np.zeros(tod.shape+(npixels_per_sample,), ops[0].matrix.dtype)
        matrix['index'] = -1
        dest = 0
        for p in ops:
            nsamples = p.matrix.shape[-2]
            matrix[...,dest:dest+nsamples,:p.matrix.shape[-1]] = p.matrix
            dest += nsamples
    else:
        raise TypeError('The operator is not a projection.')

    ny, nx = projection.shapein
    if tod.mask is None:
        mask = np.zeros(tod.shape, np.bool8)
    else:
        mask = tod.mask.copy()
    percent = func(matrix.view(int).ravel(), nx, ny, tod.T,
                   mask.view(np.int8).T, nsigma, npixels_per_sample)

    # IFORT stores .true. as -1, it conflicts with numpy algebra on booleans
    np.abs(mask, mask)

    print(strelapsed(time0, 'Deglitching ({0})'.format(method)) +
          ' (flagged samples: {0}%)'.format(percent))

    return mask
    
def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching (standard deviation).

    Each detector frame is back-projected onto the sky. The values associated
    with a sky map pixel are sigma clipped, using the standard deviation to the
    mean. In case of rejection (i. e. at a given time), each detector sample
    which contributes to the sky pixel value in the frame is masked.
    """
    return _deglitch(tod, projection, nsigma, 'std')


#-------------------------------------------------------------------------------


def deglitch_l2mad(tod, projection, nsigma=25.):
    """
    Second level deglitching (MAD).

    Each detector frame is back-projected onto the sky. The values associated
    with a sky map pixel are sigma clipped, using the MAD (median absolute
    deviation to the median). In case of rejection (i. e. at a given time),
    each detector sample which contributes to the sky pixel value in the frame
    is masked.
    """
    return _deglitch(tod, projection, nsigma, 'mad')


#-------------------------------------------------------------------------------


def filter_median(x, length=10, mask=None, partition=None, out=None):
    """
    Median filtering, O(1) in window length
    """
    from pyoperators.utils import strelapsed
    time0 = time.time()
    if out is None:
        out = x.copy()
    if mask is not None :
        mask = np.ascontiguousarray(mask, np.bool8).view(np.int8)
    elif hasattr(x, 'mask') and x.mask is not None and \
         x.mask is not np.ma.nomask:
        mask = x.mask.view(np.int8)
    else:
        mask = np.zeros(x.shape, np.int8)

    n = x.shape[-1]
    if partition is None:
        partition = (n,)
    status = tmf.filter_median(out.reshape((-1,n)).T,
        mask.reshape((-1,n)).T, length,
        np.array(partition, np.int32, ndmin=1))
    if status != 0:
        raise RuntimeError()
    print(strelapsed(time0, 'Median filtering') + ' (length={0})'.format(length))
    return out


#-------------------------------------------------------------------------------


def filter_polynomial(x, degree, mask=None, partition=None, out=None):
    """
    Filter by subtracting a fitted polynomial of arbitrary degree along
    the last dimension.

    Parameters
    ----------
    x : array_like
        The input array from which a polynomial fit will be removed.
        If it contains a mask attribute of the same shape, the masked values
        are not used to compute the polynomial fit. A True value in the mask
        means that the argument value is masked.
    degree : int
        The polynomial degree.
    mask : array_like, boolean
        Mask array compatible with argument x. If supplied, it overrides
        the first argument mask.
    partition : sequence of int
        Partition along the last dimension, to perform filtering along
        independent chunks.

    Example
    -------
    >>> y = np.array([np.arange(100.)*3, np.arange(100.)*2])
    >>> y += np.random.normal(size=y.shape)
    >>> y_filtered = filter_polynomial(y, 1)
    >>> std(y_filtered, axis=1)
    array([ 1.03496984,  1.0357345 ])

    """
    x = np.asanyarray(x)
    if out is None:
        out = x.copy()
    out_ = out.reshape((-1,x.shape[-1]))
    if mask is None:
        mask = getattr(x, 'mask', None)
    if mask is not None:
        mask = np.asarray(mask)
        mask = mask.reshape((-1,x.shape[-1]))
    if partition is None:
        partition = (x.shape[-1],)
    elif np.sum(partition) != x.shape[-1]:
        raise ValueError('The partition is not compatible with the input.')

    dest = 0
    for n in partition:
        arange = np.arange(n, dtype=x.dtype)
        for i in xrange(out_.shape[0]):
            if mask is None:
                x = arange
                y = out_[i,dest:dest+n]
            else:
                m = ~mask[i,dest:dest+n]
                x = arange[m]
                if x.size == 0:
                    continue
                y = out_[i,dest:dest+n][m]
            slope = scipy.polyfit(x, y, deg=degree)
            out_[i,dest:dest+n] -= scipy.polyval(slope, arange)
        dest += n

    return out


#-------------------------------------------------------------------------------


def filter_nonfinite(x, out=None):
    """
    Replace non-finite values with zeros in an array inplace and update
    the input mask accordingly, if present.

    Parameters
    ----------
    x : array_like
        The input array whose non-finite values will be set to zero.
        If the array contains a mask attribute, the mask entries for which the
        corresponding array entries are non-finite will be set to True.
    out : array_like, optional
        Array in which the result is stored.
    
    """
    x = np.asanyarray(x)
    
    if out is None:
        out = np.empty(x.shape)
        if type(x) is not np.ndarray:
            out = out.view(type(x))
            if out.__array_finalize__ is not None:
                out.__array_finalize__(x)
            if hasattr(x, 'mask') and isinstance(x.mask, np.ndarray):
                out.mask = x.mask.copy()
    else:
        if not isinstance(out, np.ndarray):
            raise TypeError('The output argument is not an ndarray.')

    mask = getattr(out, 'mask', None)    
    if mask is not None and mask.shape != x.shape:
        raise ValueError("The mask shape '{0}' is incompatible with that of the"
                         " input '{1}'.".format(mask.shape, x.shape))

    if x.__array_interface__['data'][0] == out.__array_interface__['data'][0]:
        if mask is None:
            tmf.processing.filter_nonfinite_inplace(x.ravel())
        else:
            tmf.processing.filter_nonfinite_mask_inplace(x.ravel(),
                mask.view(np.int8).ravel())
    else:
        if mask is None:
            tmf.processing.filter_nonfinite_outplace(x.ravel(), out.ravel())
        else:
            tmf.processing.filter_nonfinite_mask_outplace(x.ravel(),
                out.ravel(), mask.view(np.int8).ravel())

    return out


#-------------------------------------------------------------------------------


def interpolate_linear(tod, mask=None, partition=None, out=None):
    """
    Perform a linear interpolation along the last dimension under the control
    of a mask.

    Parameters
    ----------
    tod : array_like
        The input array whose masked values are to be interpolated.
        If the tod does not  contain a mask attribute, no interpolation is
        performed. A True value in the mask means that the argument value is
        masked.
    mask : array_like, boolean, optional
        Mask array compatible with argument tod. If supplied, it overrides
        the first argument mask.
    partition : sequence of int, optional
        Partition along the last dimension, to perform interpolation on
        independent chunks and avoid side-effects.
    out : array_like, optional
        Array in which the result is stored.
    
    """
    if mask is None:
        mask = getattr(tod, 'mask', None)
    if mask is None:
        return
    tod = np.asanyarray(tod)
    if out is None:
        out = tod.copy()
    else:
        out[...] = tod
    mask = np.asarray(mask)
    out_ = np.asarray(out)
    if partition is None:
        partition = (tod.shape[-1],)
    elif np.sum(partition) != tod.shape[-1]:
        raise ValueError('The partition is not compatible with the input.')

    mask = mask.reshape((-1,tod.shape[-1]))
    out_ = out_.reshape((-1,tod.shape[-1]))
    dest = 0
    for n in partition:
        x = np.arange(n)
        for i in range(out_.shape[0]):
            y = out_[i,dest:dest+n]
            invalid = mask[i,dest:dest+n]
            y[invalid] = np.interp(x[invalid], x[~invalid], y[~invalid])
        dest = dest + n

    return out
