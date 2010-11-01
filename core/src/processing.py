import numpy
import pyfits
import scipy
import tamasisfortran as tmf

__all__ = [ 'deglitch_l2std', 'deglitch_l2mad', 'filter_median', 'filter_polynomial', 'interpolate_linear', 'remove_nan' ]

def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the standard deviation to the mean.
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    if tod.mask is None:
        mask = numpy.zeros(tod.shape, numpy.bool8)
    else:
        mask = tod.mask.copy()
    tmf.deglitch_l2b_std(projection._pmatrix, nx, ny, tod.T, mask.view(numpy.int8).T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def deglitch_l2mad(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the MAD (median absolute deviation to the median)
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    if tod.mask is None:
        mask = numpy.zeros(tod.shape, numpy.bool8)
    else:
        mask = tod.mask.copy()
    tmf.deglitch_l2b_mad(projection._pmatrix, nx, ny, tod.T, mask.view(numpy.int8).T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def filter_median(tod, length=10, mask=None):
    """
    Median filtering, O(1) in window length
    """
    filtered = tod.copy()
    if mask is None:
        mask = tod.mask
        if mask is None:
            mask = numpy.zeros(tod.shape, numpy.bool8)
    else:
        mask = numpy.ascontiguousarray(mask, numpy.bool8)
    status = tmf.filter_median(filtered.T, mask.view(numpy.int8).T, length, numpy.array(tod.nsamples, dtype='int32'))
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
        x = numpy.arange(nsamples)
        slope = scipy.polyfit(x, tod[:,dest:dest+nsamples].T, deg=degree)

        for idetector in range(tod.shape[0]):
            filtered[idetector,dest:dest+nsamples] -= scipy.polyval(slope[:,idetector], x)
       
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
        x = numpy.arange(nsamples)
        for idetector in range(tod.shape[0]):
            tod_ = tod[idetector,dest:dest+nsamples]
            invalid = tod.mask[idetector,dest:dest+nsamples]
            tod_[invalid] = numpy.interp(x[invalid], x[~invalid], tod_[~invalid])
        dest = dest + nsamples
    
    return


#-------------------------------------------------------------------------------


def remove_nan(tod):
    """
    Replace NaN values in a Tod with zeros and update the mask.
    """
    if tod.mask is None:
        tod.mask = numpy.zeros(tod.shape, numpy.bool8)
    tmf.remove_nan(tod.T, tod.mask.view(numpy.int8).T)
