import numpy
import pyfits
import scipy
import tamasisfortran as tmf

__all__ = [ 'deglitch_l2std', 'deglitch_l2mad', 'filter_median', 'filter_polynomial', 'interpolate_linear' ]

def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the standard deviation to the mean.
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    mask = tod.mask.copy()
    tmf.deglitch_l2b_std(projection._pmatrix, nx, ny, tod.T, mask.T, nsigma, npixels_per_sample)
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
    mask = tod.mask.copy()
    tmf.deglitch_l2b_mad(projection._pmatrix, nx, ny, tod.T, mask.T, nsigma, npixels_per_sample)
    return mask


#-------------------------------------------------------------------------------


def filter_median(tod, length=10):
    """
    Median filtering, O(1) in window length
    """
    filtered = tod.copy()
    status = tmf.filter_median(filtered.T, length, numpy.array(tod.nsamples))
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
    Interpolate masked values of a Tod
    """
    tod = tod.copy()
    if tod.mask is None:
        return tod

    dest = 0
    for islice in range(len(tod.nsamples)):
        nsamples = tod.nsamples[islice]
        x = numpy.arange(nsamples)
        for idetector in range(tod.shape[0]):
            tod_ = tod[idetector,dest:dest+nsamples]
            valid = tod.mask[idetector,dest:dest+nsamples] == 0
            tod_[~valid] = numpy.interp(x[~valid], x[valid], tod_[valid])
        dest = dest + nsamples
    
    return tod
