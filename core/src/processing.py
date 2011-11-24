# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import numpy as np
import scipy
import tamasisfortran as tmf
from .datatypes import Tod
from acquisitionmodels import Projection

__all__ = [ 'deglitch_l2std',
            'deglitch_l2mad',
            'filter_median',
            'filter_polynomial',
            'interpolate_linear',
            'remove_nonfinite' ]

def _deglitch(tod, projection, nsigma, func):
    if isinstance(projection, Projection):
        pmatrix = projection._pmatrix
        npixels_per_sample = projection.npixels_per_sample
    elif hasattr(projection, 'operands') and \
         all([isinstance(p, Projection) for p in projection.operands]):
        ops = projection.operands
        npixels_per_sample = max([p.npixels_per_sample for p in ops])
        ndetectors = ops[0].ndetectors
        nsamples = sum([p.nsamples for p in ops])
        pmatrix = np.zeros((ndetectors, nsamples, npixels_per_sample),
                           ops[0].pmatrix.dtype)
        pmatrix['pixel'] = -1
        dest = 0
        for p in ops:
            pmatrix[:,dest:dest+p.nsamples,:p.npixels_per_sample] = p.pmatrix
            dest += p.nsamples
    else:
        raise TypeError('The operator is not a projection.')

    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    if tod.mask is None:
        mask = np.zeros(tod.shape, np.bool8)
    else:
        mask = tod.mask.copy()
    func(pmatrix.view(int).ravel(), nx, ny, tod.T, mask.view(np.int8).T,
         nsigma, npixels_per_sample)
    return mask
    
def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching (standard deviation).

    Each detector frame is back-projected onto the sky. The values associated
    with a sky map pixel are sigma clipped, using the standard deviation to the
    mean. In case of rejection (i. e. at a given time), each detector sample
    which contributes to the sky pixel value in the frame is masked.
    """
    return _deglitch(tod, projection, nsigma, tmf.deglitch_l2b_std)


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
    return _deglitch(tod, projection, nsigma, tmf.deglitch_l2b_mad)


#-------------------------------------------------------------------------------


def filter_median(tod, length=10, mask=None, nsamples=None):
    """
    Median filtering, O(1) in window length
    """
    filtered = Tod(tod)
    if mask is not None :
        mask = np.ascontiguousarray(mask, np.bool8).view(np.int8)
    elif hasattr(tod, 'mask') and tod.mask is not None and \
         tod.mask is not np.ma.nomask:
        mask = tod.mask.view(np.int8)
    else:
        mask = np.zeros(tod.shape, np.int8)

    if nsamples is None:
        nsamples = tod.shape[-1]
    nsamples_tot = tod.shape[-1]
    status = tmf.filter_median(filtered.reshape((-1,nsamples_tot)).T,
        mask.reshape((-1,nsamples_tot)).T, length,
        np.array(nsamples, np.int32, ndmin=1))
    if status != 0:
        raise RuntimeError()
    return filtered


#-------------------------------------------------------------------------------


def filter_polynomial(tod, degree, nsamples=None):
    """
    Filter by subtracting a fitted polynomial of arbitrary degree.
    """
    if nsamples is None:
        nsamples = (tod.shape[-1],)
    filtered = tod.copy()
    dest = 0
    for islice, n in enumerate(nsamples):
        x = np.arange(n)
        slope = scipy.polyfit(x, tod[:,dest:dest+n].T, deg=degree)

        for idetector in range(tod.shape[0]):
            filtered[idetector,dest:dest+n] -= \
                scipy.polyval(slope[:,idetector], x)
       
        dest = dest + n

    return filtered


#-------------------------------------------------------------------------------


def interpolate_linear(tod, nsamples=None):
    """
    In-place interpolation of masked values of a Tod
    """
    if tod.mask is None:
        return
    if nsamples is None:
        nsamples = (tod.shape[-1],)
    dest = 0
    for islice, n in enumerate(nsamples):
        x = np.arange(n)
        for idetector in range(tod.shape[0]):
            tod_ = tod[idetector,dest:dest+n]
            invalid = tod.mask[idetector,dest:dest+n]
            tod_[invalid] = np.interp(x[invalid], x[~invalid], tod_[~invalid])
        dest = dest + n
    
    return


#-------------------------------------------------------------------------------


def remove_nonfinite(tod):
    """
    Replace NaN values in a Tod with zeros and update the mask.
    """
    if tod.mask is None:
        tod.mask = np.zeros(tod.shape, np.bool8)
    tmf.remove_nonfinite_mask(tod.T, tod.mask.view(np.int8).T)
