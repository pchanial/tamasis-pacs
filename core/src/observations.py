# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import numpy as np
from mpi4py import MPI
from pyoperators.utils import isscalar, strenum
from .instruments import Instrument

__all__ = ['Observation', 'MaskPolicy']

class Observation(object):

    def __init__(self):
        self.instrument = None
        self.pointing = None
        self.policy = None
        self.slice = None

    comm_map = MPI.COMM_WORLD
    comm_tod = MPI.COMM_WORLD

    def get_ndetectors(self):
        return self.instrument.get_ndetectors()
    get_ndetectors.__doc__ = Instrument.get_ndetectors.__doc__

    def get_filter_uncorrelated(self):
        """
        Return the invNtt for uncorrelated detectors.
        """
        raise NotImplementedError()

    def get_map_header(self, resolution=None, downsampling=False):
        """
        Return the FITS header of the smallest map that encompasses
        the observation, by taking into account the instrument geometry.

        Parameters
        ----------
        resolution : float
            Sky pixel increment, in arc seconds. Default is .default_resolution.

        Returns
        -------
        header : pyfits.Header
            The resulting FITS header.
        """
        return self.instrument.get_map_header(self.pointing,
                                              resolution=resolution)

    def get_pointing_matrix(self, header, npixels_per_sample=0, method=None,
                            downsampling=False):
        """
        Return the pointing matrix.
        """
        return self.instrument.get_pointing_matrix(self.pointing, header,
            npixels_per_sample, method)

    def get_nsamples(self):
        """
        Return the number of valid pointings for each slice
        They are those for which self.pointing.removed is False
        """
        return tuple([int(np.sum(~self.pointing[s.start:s.stop].removed)) \
                      for s in self.slice])

    def get_nfinesamples(self):
        """
        Return the number of valid samplings for each slice, by taking
        into account compression factor and fine sampling.
        They are those whose self.pointing is not removed.
        """
        return tuple(np.asarray(self.get_nsamples()) * \
                     self.slice.compression_factor * \
                     self.instrument.fine_sampling_factor)

    def get_tod(self, unit=None, flatfielding=True, subtraction_mean=True):
        """
        Return the Tod from this observation
        """
        raise NotImplementedError()

    def save(self, filename, tod):
        """
        Save this observation and tod as a FITS file.
        """
        raise NotImplementedError()

    def pack(self, input, masked=False):
        return self.instrument.pack(input, masked=masked)
    pack.__doc__ = Instrument.pack.__doc__

    def unpack(self, input, masked=False):
        return self.instrument.unpack(input, masked=masked)
    unpack.__doc__ = Instrument.unpack.__doc__


#-------------------------------------------------------------------------------


class MaskPolicy(object):
    KEEP   = 0
    MASK   = 1
    REMOVE = 2
    def __init__(self, flags, values, description=None):
        self._description = description
        if isscalar(flags):
            if isinstance(flags, str):
                flags = flags.split(',')
            else:
                flags = (flags,)
        if isscalar(values):
            if isinstance(values, str):
                values = values.split(',')
            else:
                values = (values,)
        if len(flags) != len(values):
            raise ValueError('The number of policy flags is different from th' \
                             'e number of policies.')

        self._policy = []
        for flag, value in zip(flags, values):
            if flag[0] == '_':
                raise ValueError('A policy flag should not start with an unde' \
                                 'rscore.')
            value = value.strip().lower()
            choices = ('keep', 'mask', 'remove')
            if value not in choices:
                raise KeyError('Invalid policy ' + flag + "='" + value + \
                    "'. Expected policies are " + strenum(choices) + '.')
            self._policy.append({flag:value})
            setattr(self, flag, value)
        self._policy = tuple(self._policy)

    def __array__(self, dtype=int):
        conversion = { 'keep' : self.KEEP,
                       'mask' : self.MASK,
                       'remove' : self.REMOVE }
        return np.array([conversion[policy.values()[0]] \
                         for policy in self._policy], dtype=dtype)

    def __str__(self):
        str = self._description + ': ' if self._description is not None else ''
        str_ = []
        for policy in self._policy:
            str_.append(policy.values()[0] + " '" + policy.keys()[0] + "'")
        str += ', '.join(str_)
        return str
