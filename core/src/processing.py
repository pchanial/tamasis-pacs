# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import numpy as np
import pyfits
import scipy
import tamasisfortran as tmf
from .datatypes import Tod

__all__ = [ 'deglitch_l2std',
            'deglitch_l2mad',
            'filter_median',
            'filter_polynomial',
            'interpolate_linear',
            'remove_nonfinite' ]

def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching (standard deviation).

    Each detector frame is back-projected onto the sky. The values associated
    with a sky map pixel are sigma clipped, using the standard deviation to the
    mean. In case of rejection (i. e. at a given time), each detector sample
    which contributes to the sky pixel value in the frame is masked.
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    if tod.mask is None:
        mask = np.zeros(tod.shape, np.bool8)
    else:
        mask = tod.mask.copy()
    tmf.deglitch_l2b_std(projection._pmatrix, nx, ny, tod.T,
                         mask.view(np.int8).T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def deglitch_l2mad(tod, projection, nsigma=5.):
    """
    Second level deglitching (MAD).

    Each detector frame is back-projected onto the sky. The values associated
    with a sky map pixel are sigma clipped, using the MAD (median absolute
    deviation to the median). In case of rejection (i. e. at a given time),
    each detector sample which contributes to the sky pixel value in the frame
    is masked.
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    if tod.mask is None:
        mask = np.zeros(tod.shape, np.bool8)
    else:
        mask = tod.mask.copy()
    tmf.deglitch_l2b_mad(projection._pmatrix, nx, ny, tod.T,
                         mask.view(np.int8).T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def filter_median(tod, length=10, mask=None):
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

    nsamples = getattr(tod, 'nsamples', (tod.shape[-1],))
    nsamples_tot = tod.shape[-1]
    status = tmf.filter_median(filtered.reshape((-1,nsamples_tot)).T,
        mask.reshape((-1,nsamples_tot)).T, length, np.array(nsamples, np.int32))
    if status != 0:
        raise RuntimeError()
    return filtered


#-------------------------------------------------------------------------------


def filter_polynomial(tod, degree):
    """
    Filter by subtracting a fitted polynomial of arbitrary degree.
    """
    filtered = tod.copy()
    dest = 0
    for islice in range(len(tod.nsamples)):
        nsamples = tod.nsamples[islice]
        x = np.arange(nsamples)
        slope = scipy.polyfit(x, tod[:,dest:dest+nsamples].T, deg=degree)

        for idetector in range(tod.shape[0]):
            filtered[idetector,dest:dest+nsamples] -= \
                scipy.polyval(slope[:,idetector], x)
       
        dest = dest + nsamples

    return filtered


#-------------------------------------------------------------------------------


def interpolate_linear(tod):
    """
    In-place interpolation of masked values of a Tod
    """
    if tod.mask is None:
        return

    dest = 0
    for islice in range(len(tod.nsamples)):
        nsamples = tod.nsamples[islice]
        x = np.arange(nsamples)
        for idetector in range(tod.shape[0]):
            tod_ = tod[idetector,dest:dest+nsamples]
            invalid = tod.mask[idetector,dest:dest+nsamples]
            tod_[invalid] = np.interp(x[invalid], x[~invalid], tod_[~invalid])
        dest = dest + nsamples
    
    return


#-------------------------------------------------------------------------------


def remove_nonfinite(tod):
    """
    Replace NaN values in a Tod with zeros and update the mask.
    """
    if tod.mask is None:
        tod.mask = np.zeros(tod.shape, np.bool8)
    tmf.remove_nonfinite_mask(tod.T, tod.mask.view(np.int8).T)
