import numpy

from tamasis.config import get_default_dtype_float
from tamasis.numpyutils import _my_isscalar
from tamasis.unit import Quantity
from tamasis.datatypes import Map, Tod

__all__ = ['Observation', 'Instrument', 'FlatField', 'MaskPolicy', 'Pointing']

class Observation(object):

    def __init__(self):
        self.instrument = None
        self.pointing = None
        self.policy = None
        self.slice = None

    def get_ndetectors(self):
        """
        Method to get the number of detectors, accordingly to the detector mask.
        """
        return int(numpy.sum(self.instrument.detector_mask == 0))

    def get_filter_uncorrelated(self):
        """
        Method to get the invNtt for uncorrelated detectors.
        Not implemented.
        """
        raise NotImplementedError()

    def get_map_header(self, resolution=None, oversampling=True):
        """
        Method to get the map FITS header which encompasses this observation.
        Not implemented.
        """
        raise NotImplementedError()

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, method=None, oversampling=True):
        """
        Method to get the pointing matrix.
        Not implemented.
        """
        raise NotImplementedError()

    def get_nsamples(self):
        """
        Method to get the number of valid pointings for each slice
        They are those for which self.pointing.removed is False
        """
        result = []
        dest = 0
        for slice in self.slice:
            result.append(numpy.sum(~self.pointing[dest:dest+slice.nsamples_all].removed))
            dest += slice.nsamples_all
        return tuple(result)

    def get_nfinesamples(self):
        """
        Method to get the number of valid samplings for each slice, by taking
        into account compression factor and fine sampling.
        They are those whose self.pointing is not removed.
        """
        return tuple(numpy.asarray(self.get_nsamples()) * self.slice.compression_factor * self.instrument.fine_sampling_factor)

    def get_tod(self, unit=None, flatfielding=True, subtraction_mean=True):
        """
        Method to get the Tod from this observation
        Not implemented.
        """
        raise NotImplementedError()

    def save(self, filename, tod):
        """
        Method to save this observation.
        Not implemented.
        """
        raise NotImplementedError()

    def unpack(self, tod):
        """
        Convert a Tod which only includes the valid detectors into 
        another Tod which contains all the detectors under the control
        of the detector mask

        Parameters
        ----------
        tod : Tod
              rank-2 Tod to be unpacked.

        Returns
        -------
        unpacked_tod : Tod
              This will be a new view object if possible; otherwise, it will
              be a copy.

        See Also
        --------
        pack: inverse method

        """
        if numpy.rank(tod) != 2:
            raise ValueError('The input Tod is not packed.')
        ndetectors = self.instrument.detector_mask.size
        nvalids  = self.get_ndetectors()
        nsamples = tod.shape[-1]

        newshape = numpy.concatenate((self.instrument.shape, (nsamples,)))
        if nvalids != tod.shape[0]:
            raise ValueError("The detector mask has a number of valid detectors '" + str(numpy.sum(mask == 0)) + 
                             "' incompatible with the packed input '" + str(tod.shape[0]) + "'.")

        # return a view if all detectors are valid
        if nvalids == ndetectors:
            return tod.reshape(newshape)

        # otherwise copy the detector timelines one by one
        mask = self.instrument.detector_mask.ravel()
        utod = Tod.zeros(newshape, nsamples=tod.nsamples, unit=tod.unit, dtype=tod.dtype, 
                         mask=numpy.ones(newshape, dtype='int8'))
        rtod = utod.reshape((ndetectors, nsamples))
        
        i = 0
        for iall in range(mask.size):
            if mask[iall] != 0:
                rtod[iall,:] = 0
                rtod.mask[iall,:] = 1
                continue
            rtod[iall,:] = tod[i,:]
            if tod.mask is not None:
                rtod.mask[iall,:] = tod.mask[i,:]
            else:
                rtod.mask[iall,:] = 0
            i += 1
        return utod

    def pack(self, tod):
        """
        Convert a Tod which only includes all the detectors into 
        another Tod which contains only the valid ones under the control
        of the detector mask

        Parameters
        ----------
        tod : Tod
              rank-2 Tod to be packed.

        Returns
        -------
        packed_tod : Tod
              This will be a new view object if possible; otherwise, it will
              be a copy.

        See Also
        --------
        unpack: inverse method

        """
        if numpy.rank(tod) != numpy.rank(self.instrument.detector_mask)+1:
            raise ValueError('The input Tod is not unpacked.')
        ndetectors = self.instrument.detector_mask.size
        nvalids = self.get_ndetectors()
        nsamples = tod.shape[-1]

        newshape = (nvalids, nsamples)
        if tod.shape[0:-1] != self.instrument.shape:
            raise ValueError("The detector mask has a shape '" + str(self.instrument.shape) + 
                             "' incompatible with the unpacked input '" + str(tod.shape[0:-1]) + "'.")

        # return a view if all detectors are valid
        if nvalids == ndetectors:
            return tod.reshape((ndetectors, nsamples))

        # otherwise copy the detector timelines one by one
        mask = self.instrument.detector_mask.ravel()
        ptod = Tod.empty(newshape, nsamples=tod.nsamples, unit=tod.unit, dtype=tod.dtype, 
                         mask=numpy.empty(newshape, dtype='int8'))
        rtod = tod.reshape((ndetectors, nsamples))
        
        i = 0
        for iall in range(mask.size):
            if mask[iall] != 0:
                 continue
            ptod[i,:] = rtod[iall,:]
            if tod.mask is not None:
                ptod.mask[i,:] = rtod.mask[iall,:]
            else:
                ptod.mask[i,:] = 0
            i += 1
        return ptod


#-------------------------------------------------------------------------------


class Instrument(object):
    """
    Class storing information about the instrument.
    Attributes include: name, ndetectors, shape, detector_mask
    """
    def __init__(self, name, detector_mask):
        self.name = name
        self.detector_mask = Map(detector_mask, origin='upper', dtype='uint8', copy=True)
        self.shape = detector_mask.shape


#-------------------------------------------------------------------------------


class FlatField(object):
    def __init__(self, optical, detector):
        self.optical = Map(optical, origin='upper')
        self.detector = Map(detector, origin='upper')

    @property
    def total(self):
        return Map(self.optical * self.detector, origin='upper')


#-------------------------------------------------------------------------------


class MaskPolicy(object):
    KEEP   = 0
    MASK   = 1
    REMOVE = 2
    def __init__(self, flags, values, description=None):
        self._description = description
        if _my_isscalar(flags):
            if isinstance(flags, str):
                flags = flags.split(',')
            else:
                flags = (flags,)
        if _my_isscalar(values):
            if isinstance(values, str):
                values = values.split(',')
            else:
                values = (values,)
        if len(flags) != len(values):
            raise ValueError('The number of policy flags is different from the number of policies.')

        self._policy = []
        for flag, value in zip(flags, values):
            if flag[0] == '_':
                raise ValueError('A policy flag should not start with an underscore.')
            value = value.strip().lower()
            if value not in ('keep', 'mask', 'remove'):
                raise KeyError("Invalid policy "+flag+"='" + value + "'. Valid policies are 'keep', 'mask' or 'remove'.")
            self._policy.append({flag:value})
            setattr(self, flag, value)
        self._policy = tuple(self._policy)

    def __array__(self, dtype=int):
        conversion = { 'keep':self.KEEP, 'mask':self.MASK, 'remove':self.REMOVE }
        return numpy.array([conversion[policy.values()[0]] for policy in self._policy], dtype=dtype)

    def __str__(self):
        str = self._description + ': ' if self._description is not None else ''
        str_ = []
        for policy in self._policy:
            str_.append(policy.values()[0] + " '" + policy.keys()[0] + "'")
        str += ', '.join(str_)
        return str


#-------------------------------------------------------------------------------


class Pointing(numpy.recarray):
    INSCAN     = 1
    TURNAROUND = 2
    OTHER      = 3
    def __new__(cls, time, ra, dec, pa, info=None, masked=None, removed=None, nsamples=None, dtype=None):

        if nsamples is None:
            nsamples = (time.size,)
        else:
            if _my_isscalar(nsamples):
                nsamples = (nsamples,)
            else:
                nsamples = tuple(nsamples)
        nsamples_tot = numpy.sum(nsamples)

        if info is None:
            info = Pointing.INSCAN

        if masked is None:
            masked = False

        if removed is None:
            removed = False

        if dtype is None:
            ftype = get_default_dtype_float()
            dtype = numpy.dtype([('time', ftype), ('ra', ftype), ('dec', ftype), ('pa', ftype), ('info', numpy.int64),
                                 ('masked', numpy.bool8), ('removed', numpy.bool8)])


        time    = numpy.asarray(time)
        ra      = numpy.asarray(ra)
        dec     = numpy.asarray(dec)
        pa      = numpy.asarray(pa)
        info    = numpy.asarray(info)
        masked  = numpy.asarray(masked)
        removed = numpy.asarray(removed)
        if numpy.any(map(lambda x: x.size not in (1, nsamples_tot), [ra,dec,pa,info,masked,removed])):
            raise ValueError('The pointing inputs do not have the same size.')

        result = numpy.recarray(nsamples_tot, dtype=dtype)
        result.time    = time
        result.ra      = ra
        result.dec     = dec
        result.pa      = pa
        result.info    = info
        result.masked  = masked
        result.removed = removed
        result = result.view(cls)
        result.nsamples = nsamples
        return result

    @property
    def velocity(self):
        ra = Quantity(self.ra, 'deg')
        dec = Quantity(self.dec, 'deg')
        dra  = numpy.diff(ra)
        ddec = numpy.diff(dec)
        dtime = Quantity(numpy.diff(self.time), 's')
        vel = numpy.sqrt((dra*numpy.cos(dec[0:-1].inunit('rad')))**2 + ddec**2) / dtime
        vel.unit = 'arcsec/s'
        u = vel._unit
        vel = numpy.append(vel, vel[-1])
        # BUG: append eats the unit...
        vel._unit = u
        return vel
