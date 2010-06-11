import numpy
import pyfits
import re
import tamasisfortran as tmf
from datatypes import *

__all__ = [ 'MadMap1Observation' ]

class MadMap1Observation(object):
    """Class for the handling of an observation in the MADMAP1 format"""
    def __init__(self, todfile, invnttfile, mapmaskfile, convert, ndetectors, missing_value=None):
        nslices, status = tmf.read_madmap1_nslices(invnttfile, ndetectors)
        if (status != 0): raise RuntimeError()
        self.npixels_per_sample, nsamples, status = tmf.read_madmap1_info(todfile, invnttfile, convert, ndetectors, nslices)
        if (status != 0): raise RuntimeError()

        self.todfile = todfile
        self.invnttfile = invnttfile
        m=re.search(r'(?P<filename>.*)\[(?P<extname>\w+)\]$', mapmaskfile)
        if m is None:
            mask = pyfits.fitsopen(mapmaskfile)[0].data
        else:
            filename = m.group('filename')
            extname  = m.group('extname')
            mask = pyfits.fitsopen(filename)[extname].data
        if mask is None:
            raise IOError('HDU '+mapmaskfile+' has no data.')
        self.mapmask = numpy.zeros(mask.shape, dtype='int8')
        if missing_value is None:
            self.mapmask[mask != 0] = 1
        elif numpy.isnan(missing_value):
            self.mapmask[numpy.isnan(mask)] = 1
        elif numpy.isinf(missing_value):
            self.mapmask[numpy.isinf(mask)] = 1
        else:
            self.mapmask[mask == missing_value] = 1
        self.convert = convert
        self.ndetectors = ndetectors
        self.nsamples = tuple(nsamples)
        self.nfinesamples = tuple(nsamples)

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, oversampling=False):
        if npixels_per_sample is not None and npixels_per_sample != self.npixels_per_sample:
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
        header.update('naxis1', numpy.sum(self.mapmask == 0))

        tod = Tod.zeros((self.ndetectors,self.nsamples))
        sizeofpmatrix = self.npixels_per_sample * numpy.sum(self.nsamples) * self.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.read_madmap1(self.todfile, self.invnttfile, self.convert, self.npixels_per_sample, tod.T, pmatrix)
        if (status != 0): raise RuntimeError()
        return pmatrix, header, self.ndetectors, self.nsamples, self.npixels_per_sample

    def get_tod(self):
        tod = Tod.zeros((self.ndetectors,self.nsamples))
        sizeofpmatrix = self.npixels_per_sample * numpy.sum(self.nsamples) * self.ndetectors
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.read_madmap1(self.todfile, self.invnttfile, self.convert, self.npixels_per_sample, tod.T, pmatrix)
        if (status != 0): raise RuntimeError()
        return tod
