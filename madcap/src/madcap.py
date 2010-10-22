import numpy
import pyfits
import re
from tamasis.core import *
from tamasis.observations import *

__all__ = [ 'MadMap1Observation' ]

class MadMap1Observation(Observation):
    """Class for the handling of an observation in the MADMAP1 format"""
    def __init__(self, todfile, invnttfile, mapmaskfile, convert, ndetectors, missing_value=None):

        # Get information from files
        nslices, status = tmf.madmap1_nslices(invnttfile, ndetectors)
        if status != 0: raise RuntimeError()
        npixels_per_sample, nsamples, ncorrelations, status = tmf.madmap1_info(todfile, invnttfile, convert, ndetectors, nslices)
        if status != 0: raise RuntimeError()

        m=re.search(r'(?P<filename>.*)\[(?P<extname>\w+)\]$', mapmaskfile)
        if m is None:
            mask = pyfits.fitsopen(mapmaskfile)[0].data
        else:
            filename = m.group('filename')
            extname  = m.group('extname')
            mask = pyfits.fitsopen(filename)[extname].data
        if mask is None:
            raise IOError('HDU '+mapmaskfile+' has no data.')
        mapmask = numpy.zeros(mask.shape, dtype='int8')
        if missing_value is None:
            mapmask[mask != 0] = 1
        elif numpy.isnan(missing_value):
            mapmask[numpy.isnan(mask)] = 1
        elif numpy.isinf(missing_value):
            mapmask[numpy.isinf(mask)] = 1
        else:
            mapmask[mask == missing_value] = 1

        # Store instrument information
        self.instrument = Instrument('Unknown', numpy.zeros((ndetectors,)))

        # Store observation information
        class MadMap1ObservationInfo(object):
            pass
        self.info = MadMap1ObservationInfo()
        self.info.todfile = todfile
        self.info.invnttfile = invnttfile
        self.info.ncorrelations = ncorrelations
        self.info.npixels_per_sample = npixels_per_sample
        self.info.mapmaskfile = mapmaskfile
        self.info.convert = convert
        self.info.missing_value = missing_value
        self.info.mapmask = mapmask

        # Store slice information
        self.slice = numpy.recarray(nslices, dtype=[('nsamples', int), ('invnttfile', 'S256')])
        self.slice.nsamples = nsamples
        self.slice.nfinesamples = nsamples
        self.slice.invnttfile = [invnttfile+'.'+str(i) for i in range(nslices)]        

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, method=None, oversampling=False):
        """
        Method to get the pointing matrix.
        """
        if npixels_per_sample is not None and npixels_per_sample != self.info.npixels_per_sample:
            raise ValueError('The npixels_per_sample value is incompatible with the MADMAP1 file.')
        if header is not None:
            raise ValueError('The map header cannot be specified for MADmap1 observations.')
        if resolution is not None:
            raise ValueError('The map resolution cannot be specified for MADmap1 observations.')
        header = pyfits.Header()
        header.update('simple', True)
        header.update('bitpix', -64)
        header.update('extend', True)
        header.update('naxis', 1)
        header.update('naxis1', numpy.sum(self.info.mapmask == 0))

        tod = Tod.zeros((self.instrument.ndetectors,self.slice.nsamples))
        sizeofpmatrix = self.info.npixels_per_sample * numpy.sum(self.slice.nsamples) * self.instrument.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.madmap1_read_tod(self.info.todfile, self.info.invnttfile, self.info.convert, self.info.npixels_per_sample, tod.T, pmatrix)
        if status != 0: raise RuntimeError()
        return pmatrix, header, self.slice.nsamples, self.info.npixels_per_sample

    def get_tod(self):
        """
        Method to get the Tod from this observation
        """
        tod = Tod.zeros((self.instrument.ndetectors, self.slice.nsamples))
        sizeofpmatrix = self.info.npixels_per_sample * numpy.sum(self.slice.nsamples) * self.instrument.ndetectors
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=int)
        status = tmf.madmap1_read_tod(self.info.todfile, self.info.invnttfile, self.info.convert, self.info.npixels_per_sample, tod.T, pmatrix)
        if status != 0: raise RuntimeError()
        return tod

    def get_filter_uncorrelated(self):
        """
        Method to get the invNtt for uncorrelated detectors.
        """
        data, nsamples, status = tmf.madmap1_read_filter(self.info.invnttfile, self.info.convert, self.info.ncorrelations, self.instrument.ndetectors, self.slice.size)
        if status != 0: raise RuntimeError()
        return data.T
