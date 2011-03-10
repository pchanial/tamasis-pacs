import numpy as np

from . import var
from .datatypes import FitsArray, Map, Tod, create_fitsheader
from .numpyutils import _my_isscalar
from .quantity import Quantity
from .stringutils import strenum

__all__ = ['Observation', 'Instrument', 'MaskPolicy', 'Pointing']

class Observation(object):

    def __init__(self):
        self.instrument = None
        self.pointing = None
        self.policy = None
        self.slice = None

    def get_ndetectors(self):
        """
        Return the number of detectors, according to the detector mask.
        """
        return int(np.sum(~self.instrument.detector.removed))

    def get_filter_uncorrelated(self):
        """
        Return the invNtt for uncorrelated detectors.
        """
        raise NotImplementedError()

    def get_map_header(self, resolution=None, oversampling=True):
        """
        Return the map FITS header which encompasses this observation.
        """
        raise NotImplementedError()

    def get_pointing_matrix(self, header, resolution, npixels_per_sample,
                            method=None, oversampling=True):
        """
        Return the pointing matrix.
        """
        raise NotImplementedError()

    def get_nsamples(self):
        """
        Return the number of valid pointings for each slice
        They are those for which self.pointing.removed is False
        """
        result = []
        dest = 0
        for slice in self.slice:
            result.append(
                int(np.sum(~self.pointing[dest:dest+slice.nsamples_all].removed)))
            dest += slice.nsamples_all
        return tuple(result)

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
        if np.rank(tod) != 2:
            raise ValueError('The input Tod is not packed.')
        ndetectors = self.instrument.detector.size
        nvalids  = self.get_ndetectors()
        nsamples = tod.shape[-1]

        newshape = np.concatenate((self.instrument.detector.shape, (nsamples,)))
        if nvalids != tod.shape[0]:
            raise ValueError("The detector mask has a number of valid detecto" \
                "rs '" + str(np.sum(mask == 0)) + "' incompatible with the pac" \
                "ked input '" + str(tod.shape[0]) + "'.")

        # return a view if all detectors are valid
        if nvalids == ndetectors:
            return tod.reshape(newshape)

        # otherwise copy the detector timelines one by one
        mask = self.instrument.detector.removed.ravel()
        utod = Tod.empty(newshape, nsamples=tod.nsamples, unit=tod.unit,
                         derived_units=tod.derived_units, dtype=tod.dtype, 
                         mask=np.empty(newshape, dtype='int8'))
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
        if np.rank(tod) != np.rank(self.instrument.detector.removed)+1:
            raise ValueError('The input Tod is not unpacked.')
        ndetectors = self.instrument.detector.size
        nvalids = self.get_ndetectors()
        nsamples = tod.shape[-1]

        newshape = (nvalids, nsamples)
        if tod.shape[0:-1] != self.instrument.shape:
            raise ValueError("The detector mask has a shape '" + \
                str(self.instrument.shape) + "' incompatible with the unpacke" \
                "d input '" + str(tod.shape[0:-1]) + "'.")

        # return a view if all detectors are valid
        if nvalids == ndetectors:
            return tod.reshape((ndetectors, nsamples))

        # otherwise copy the detector timelines one by one
        mask = self.instrument.detector.removed.ravel()
        ptod = Tod.empty(newshape, nsamples=tod.nsamples, unit=tod.unit,
                         derived_units=tod.derived_units, dtype=tod.dtype, 
                         mask=np.empty(newshape, dtype='int8'))
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
    Attributes include: name, detector
    """
    def __init__(self, name, shape, dtype=None):

        if dtype is None:
            dtype = [('masked', np.bool8), ('removed', np.bool8)]
        self.name = name
        self.detector = Map.zeros(shape, dtype=dtype, origin='upper')

    def __getattribute__(self, name):
        table = {'detector_time_constant':'time_constant',
                 'detector_mask':'removed', 'detector_bad':'masked'}
        if name in table:
            print("\nXXXXXX\nPlease update your script: 'obs.instrument." + name + "' is DEPRECATED and will be removed in the next version of Tamasis. Use 'obs.instrument.detector." + table[name] + "' instead.\nXXXXXX\n")
            return self.detector[table[name]]
        return object.__getattribute__(self, name)


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
                    "'. Expected policies are " + strenum(choices, 'or') + '.')
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


#-------------------------------------------------------------------------------

POINTING_DTYPE = [('time', var.FLOAT_DTYPE), ('ra', var.FLOAT_DTYPE),
                  ('dec', var.FLOAT_DTYPE), ('pa', var.FLOAT_DTYPE),
                  ('info', np.int64), ('masked', np.bool8),
                  ('removed', np.bool8)]

class Pointing(FitsArray):
    INSCAN     = 1
    TURNAROUND = 2
    OTHER      = 3
    def __new__(cls, time, ra, dec, pa, info=None, masked=None, removed=None,
                nsamples=None, dtype=None):

        if info is None:
            info = Pointing.INSCAN

        if masked is None:
            masked = False

        if removed is None:
            removed = False

        if dtype is None:
            dtype = POINTING_DTYPE

        time    = np.asarray(time)
        ra      = np.asarray(ra)
        dec     = np.asarray(dec)
        pa      = np.asarray(pa)
        info    = np.asarray(info)
        masked  = np.asarray(masked)
        removed = np.asarray(removed)

        if nsamples is None:
            nsamples = (time.size,)
        else:
            if _my_isscalar(nsamples):
                nsamples = (nsamples,)
            else:
                nsamples = tuple(nsamples)
        nsamples_tot = np.sum(nsamples)

        if np.any([x.size not in (1, nsamples_tot) \
                   for x in [ra, dec, pa, info, masked, removed]]):
            raise ValueError('The pointing inputs do not have the same size.')

        result = FitsArray.zeros(nsamples_tot, dtype=dtype)
        result.time    = time
        result.ra      = ra
        result.dec     = dec
        result.pa      = pa
        result.info    = info
        result.masked  = masked
        result.removed = removed
        result.header = create_fitsheader(None, naxis=result.size)
        result._unit = {}
        result._derived_units = {}
        result = result.view(cls)
        result.nsamples = nsamples

        return result.view(cls)

    def __array_finalize__(self, array):
        FitsArray.__array_finalize__(self, array)
        nsamples = getattr(array, 'nsamples', None)
        self.nsamples = (self.size,) if nsamples is None else nsamples

    def __getattr__(self, name):
        if self.dtype.names is None or name not in self.dtype.names:
            raise AttributeError("'" + self.__class__.__name__ + "' object ha" \
                                 "s no attribute '" + name + "'")
        return self[name].magnitude

    @property
    def velocity(self):
        if self.size <= 1:
            return Quantity(np.nan, 'arcsec/s')
        ra = Quantity(self.ra, 'deg')
        dec = Quantity(self.dec, 'deg')
        dra  = np.diff(ra)
        ddec = np.diff(dec)
        dtime = Quantity(np.diff(self.time), 's')
        vel = np.sqrt((dra * np.cos(dec[0:-1].tounit('rad')))**2 + ddec**2) \
              / dtime
        vel.inunit('arcsec/s')
        u = vel._unit
        vel = np.append(vel, vel[-1])
        # BUG: append eats the unit...
        vel._unit = u
        return vel


#-------------------------------------------------------------------------------


def create_scan(ra0, dec0, scan_acceleration, sampling_period, scan_angle=0.,
                scan_length=30., scan_nlegs=3, scan_step=20., scan_speed=10.,
                dtype=None):
    """
    compute the pointing timeline of the instrument reference point
    from the description of a scan map
    Authors: R. Gastaud
    """
    
    scan_angle -= 90.

    scan_acceleration = float(scan_acceleration)
    if scan_acceleration <= 0:
        raise ValueError('Input scan_acceleration must be strictly positive.')

    sampling_period = float(sampling_period)
    if sampling_period <= 0:
        raise ValueError('Input sampling_period must be strictly positive.')
    
    scan_length = float(scan_length)
    if scan_length <= 0:
        raise ValueError('Input scan_length must be strictly positive.')

    scan_nlegs = int(scan_nlegs)
    if scan_nlegs <= 0:
        raise ValueError('Input scan_nlegs must be strictly positive.')
    
    scan_step = float(scan_step)
    if scan_step <= 0: 
        raise ValueError('Input scan_step must be strictly positive.')

    scan_speed = float(scan_speed)
    if scan_speed <= 0:
        raise ValueError('Input scan_speed must be strictly positive.')

    # compute the different times and the total number of points
    # acceleration time at the beginning of a leg, and deceleration time
    # at the end
    extra_time1 = scan_speed / scan_acceleration
    # corresponding length 
    extralength = 0.5 * scan_acceleration * extra_time1 * extra_time1
    # Time needed to go from a scan line to the next 
    extra_time2 = np.sqrt(scan_step / scan_acceleration)
    # Time needed to turn around (not used)
    turnaround_time = 2 * (extra_time1 + extra_time2)

    # Time needed to go along the scanline at constant speed
    line_time = scan_length / scan_speed
    # Total time for a scanline
    full_line_time = extra_time1 + line_time + extra_time1 + extra_time2 + \
                     extra_time2
    # Total duration of the observation
    total_time = full_line_time * scan_nlegs - 2 * extra_time2

    # Number of samples
    nsamples = int(np.ceil(total_time / sampling_period))

    # initialization
    time          = np.zeros(nsamples)
    latitude      = np.zeros(nsamples)
    longitude     = np.zeros(nsamples)
    infos         = np.zeros(nsamples, dtype=int)
    line_counters = np.zeros(nsamples, dtype=int)

    # Start of computations, alpha and delta are the longitude and
    # latitide in arc seconds in the referential of the map.
    signe = 1
    delta = -extralength - scan_length/2.
    alpha = -scan_step * (scan_nlegs-1)/2.
    alpha0 = alpha
    line_counter = 0
    working_time = 0.

    for i in range(nsamples):
        info = 0
   
        # check if new line
        if working_time > full_line_time:
            working_time = working_time - full_line_time
            signe = -signe
            line_counter = line_counter + 1
            alpha = -scan_step * (scan_nlegs-1) / 2. + line_counter * scan_step
            alpha0 = alpha
   
        # acceleration at the beginning of a scan line to go from 0 to the
        # scan_speed. 
        if working_time < extra_time1:
            delta = -signe * (extralength + scan_length/2) + signe * 0.5 * \
                    scan_acceleration * working_time * working_time
            info  = Pointing.TURNAROUND
   
        # constant speed
        if working_time >=  extra_time1 and \
           working_time < extra_time1 + line_time:
            delta = signe * (-scan_length / 2. + (working_time - extra_time1) *\
                    scan_speed)
            info  = Pointing.INSCAN
 
        # Deceleration at then end of the scanline to stop
        if working_time >= extra_time1 + line_time and \
           working_time < extra_time1 + line_time + extra_time1:
            dt = working_time - extra_time1 - line_time
            delta = signe * (scan_length / 2. + scan_speed * dt - \
                    0.5 * scan_acceleration * dt * dt)
            info  = Pointing.TURNAROUND
  
        # Acceleration to go toward the next scan line
        if working_time >= 2 * extra_time1 + line_time and \
           working_time < 2 * extra_time1 + line_time + extra_time2:
            dt = working_time - 2 * extra_time1 - line_time
            alpha = alpha0 + 0.5 * scan_acceleration * dt * dt
            info  = Pointing.TURNAROUND
   
        # Deceleration to stop at the next scan line
        if working_time >= 2 * extra_time1 + line_time + extra_time2 and \
           working_time < full_line_time:
            dt = working_time - 2 * extra_time1 - line_time - extra_time2
            speed = scan_acceleration * extra_time2
            alpha = (alpha0 + scan_step / 2.) + speed * dt - \
                    0.5 * scan_acceleration * dt * dt
            info  = Pointing.TURNAROUND

        time[i] = i * sampling_period
        infos[i] = info
        latitude[i] = delta
        longitude[i] = alpha
        line_counters[i] = line_counter
        working_time = working_time + sampling_period

    # Convert the longitude and latitude *expressed in degrees) to ra and dec
    ra, dec = _change_coord(ra0, dec0, scan_angle, longitude / 3600.,
                            latitude / 3600.)

    scan = Pointing(time, ra, dec, 0., infos, dtype=dtype)
    header = create_fitsheader(scan)
    header.update('ra', ra0)
    header.update('dec', dec0)
    header.update('HIERARCH scan_angle', scan_angle)
    header.update('HIERARCH scan_length', scan_length)
    header.update('HIERARCH scan_nlegs', scan_nlegs)
    header.update('HIERARCH scan_step', scan_step)
    header.update('HIERARCH scan_speed', scan_speed)
    header.update('HIERARCH scan_acceleration', scan_acceleration)
    scan.header = header
    return scan


#-------------------------------------------------------------------------------


def _change_coord(ra0, dec0, pa0, lon, lat):
    """
    Transforms the longitude and latitude coordinates expressed in the
    native coordinate system attached to a map into right ascension and
    declination.
    Author: R. Gastaud
    """
    from numpy import arcsin, arctan2, cos, deg2rad, mod, sin, rad2deg

    # Arguments in radian
    lambd = deg2rad(lon)
    beta  = deg2rad(lat) 
    alpha = deg2rad(ra0)
    delta = deg2rad(dec0)
    eps   = deg2rad(pa0)

    # Cartesian coordinates from longitude and latitude
    z4 = sin(beta)
    y4 = sin(lambd)*cos(beta)
    x4 = cos(lambd)*cos(beta)

    # rotation about the x4 axe of -eps 
    x3 =  x4
    y3 =  y4*cos(eps) + z4*sin(eps)
    z3 = -y4*sin(eps) + z4*cos(eps)
    
    # rotation about the axis Oy2, angle delta
    x2 = x3*cos(delta) - z3*sin(delta)
    y2 = y3
    z2 = x3*sin(delta) + z3*cos(delta)

    # rotation about the axis Oz1, angle alpha
    x1 = x2*cos(alpha) - y2*sin(alpha)
    y1 = x2*sin(alpha) + y2*cos(alpha)
    z1 = z2

    # Compute angles from cartesian coordinates
    # it is the only place where we can get nan with arcsinus
    dec = rad2deg(arcsin(np.clip(z1, -1., 1.)))
    ra  = mod(rad2deg(arctan2(y1, x1)), 360.)

    return ra, dec

