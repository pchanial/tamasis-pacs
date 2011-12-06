# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import gc
import numpy as np
import pyfits
import re

from . import var
from . import tmf
from tamasis.datatypes import Tod
from tamasis.instruments import Instrument
from tamasis.observations import Observation

__all__ = [ 'MadMap1Observation' ]

class MadMap1Observation(Observation):
    """Class for the handling of an observation in the MADMAP1 format"""
    def __init__(self, todfile, invnttfile, mapmaskfile, convert, ndetectors,
                 missing_value=None, comm_tod=None):

        # Get information from files
        nslices, status = tmf.madmap1_nslices(invnttfile, ndetectors)
        if status != 0: raise RuntimeError()
        npixels_per_sample, nsamples, ncorrelations, status = tmf.madmap1_info(\
            todfile, invnttfile, convert, ndetectors, nslices)
        if status != 0: raise RuntimeError()

        m=re.search(r'(?P<filename>.*)\[(?P<extname>\w+)\]$', mapmaskfile)
        if m is None:
            mask = pyfits.open(mapmaskfile)[0].data
        else:
            filename = m.group('filename')
            extname  = m.group('extname')
            mask = pyfits.open(filename)[str(extname)].data #XXX Python3
        if mask is None:
            raise IOError('HDU '+mapmaskfile+' has no data.')
        mapmask = np.zeros(mask.shape, dtype=bool)
        if missing_value is None:
            mapmask[mask != 0] = True
        elif np.isnan(missing_value):
            mapmask[np.isnan(mask)] = True
        elif np.isinf(missing_value):
            mapmask[np.isinf(mask)] = True
        else:
            mapmask[mask == missing_value] = True

        # Store instrument information
        self.instrument = Instrument('Unknown', (ndetectors,))

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
        self.slice = np.recarray(nslices, dtype=[('nsamples_all', int),
                                                 ('invnttfile', 'S256')])
        self.slice.nsamples_all = nsamples
        self.slice.nfinesamples = nsamples
        self.slice.invnttfile = [invnttfile+'.'+str(i) for i in range(nslices)]

        # Store pointing information
        self.pointing = np.recarray(np.sum(nsamples), [('removed', np.bool_)])
        self.pointing.removed = False

        # Store communicator information
        self.comm_tod = comm_tod or var.comm_tod
        if self.comm_tod.Get_size() > 1:
            raise NotImplementedError('The parallelisation of the TOD is not ' \
                                      'implemented')

    def get_map_header(self, resolution=None, **keywords):
        header = pyfits.Header()
        header.update('simple', True)
        header.update('bitpix', -64)
        header.update('extend', True)
        header.update('naxis', 1)
        header.update('naxis1', np.sum(~self.info.mapmask))
        return header

    def get_pointing_matrix(self, header, npixels_per_sample, method=None,
                            oversampling=False, islice=None):
        """
        Method to get the pointing matrix.
        """
        if npixels_per_sample not in (0, self.info.npixels_per_sample):
            raise ValueError('The npixels_per_sample value is incompatible wi' \
                             'th the MADMAP1 file.')
        if method is None:
            method = 'default'
        method = method.lower()
        if method != 'default':
            raise ValueError("Invalid pointing matrix method '" + method + "'.")

        ndetectors = self.get_ndetectors()
        if islice is None:
            nsamples = int(np.sum(self.get_nsamples()))
            islice = -1
        else:
            nsamples = self.get_nsamples()[islice]
        tod = np.empty((ndetectors, nsamples))

        shape = (ndetectors, nsamples, self.info.npixels_per_sample)
        dtype = [('weight', 'f4'), ('pixel', 'i4')]
        if npixels_per_sample != 0:
            print('Info: Allocating ' + str(np.product(shape) / 2**17) + ' MiB '
                  'for the pointing matrix.')
        try:
            pmatrix = np.empty(shape, dtype).view(np.recarray)
        except MemoryError:
            gc.collect()
            pmatrix = np.empty(shape, dtype).view(np.recarray)
        status = tmf.madmap1_read_tod(self.info.todfile, self.info.invnttfile,
            self.info.convert, self.info.npixels_per_sample, islice + 1, tod.T,
            pmatrix.ravel().view(np.int64))
        if status != 0: raise RuntimeError()
        return pmatrix, method, None, None

    def get_tod(self, unit=None):
        """
        Method to get the Tod from this observation
        """
        tod = Tod.empty((self.get_ndetectors(), np.sum(self.get_nsamples())))
        sizeofpmatrix = self.info.npixels_per_sample * tod.size
        pmatrix = np.zeros(sizeofpmatrix, dtype=int)
        status = tmf.madmap1_read_tod(self.info.todfile, self.info.invnttfile,
            self.info.convert, self.info.npixels_per_sample, 0, tod.T, pmatrix)
        if status != 0: raise RuntimeError()
        if unit is not None:
            tod.unit = unit
        return tod

    def get_filter_uncorrelated(self, **keywords):
        """
        Method to get the invNtt for uncorrelated detectors.
        """
        data, nsamples, status = tmf.madmap1_read_filter(self.info.invnttfile,
            self.info.convert, self.info.ncorrelations, self.get_ndetectors(),
            self.slice.size)
        if status != 0: raise RuntimeError()
        return data.T

    def get_nsamples(self):
        nsamples = np.sum(~self.pointing.removed)
        if nsamples != np.sum(self.pointing.size):
            raise ValueError("The pointing attribute 'removed' cannot be set " \
                             'to True.')
        return self.slice.nsamples_all
