from __future__ import division

import kapteyn
import numpy as np
from matplotlib import pyplot
from pyoperators.utils import isscalar

from . import var
from .datatypes import FitsArray, Map, create_fitsheader
from .quantity import Quantity
from .wcsutils import angle_lonlat, barycenter_lonlat

__all__ = ['Pointing', 'plot_scan']
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
            if isscalar(nsamples):
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
        result.header = create_fitsheader(result.size)
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
        if self.size == 0:
            return Quantity([], 'arcsec/s')
        elif self.size == 1:
            return Quantity([np.nan], 'arcsec/s')
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


def create_scan(ra0, dec0, length, step, sampling_period, speed, acceleration,
                nlegs=3, angle=0., instrument_angle=45, cross_scan=True,
                dtype=POINTING_DTYPE):
    """
    Return a sky scan.

    The output is a Pointing instance that can be handed to PacsSimulation to
    create a simulation.

    Parameters
    ----------
    ra0 : float
        Right Ascension of the scan center.
    dec0 : float
        Declination of the scan center.
    length : float
        Length of the scan lines, in arcseconds.
    step : float
        Separation between scan legs, in arcseconds.
    sampling_period : float
        Duration between two pointings.
    speed : float
        Scan speed, in arcsec/s.
    acceleration : float
        Acceleration, in arcsec/s^2.
    nlegs : integer
        Number of scan legs.
    angle : float
        Angle between the scan line direction and the North minus 90 degrees,
        in degrees.
    instrument_angle : float
        Angle between the scan line direction and the instrument second axis,
        in degrees.
    cross_scan : boolean
        If true, a cross-scan is appended to the pointings.
    dtype : flexible dtype
        Pointing data type.
    """

    scan = _create_scan(ra0, dec0, length, step, sampling_period, speed,
                        acceleration, nlegs, angle, dtype)
    if cross_scan:
        cross = _create_scan(ra0, dec0, length, step, sampling_period, speed,
                             acceleration, nlegs, angle + 90, dtype)
        cross.time += scan.time[-1] + sampling_period
        scan, scan.header = Pointing(np.hstack([scan.time, cross.time]),
                                     np.hstack([scan.ra, cross.ra]),
                                     np.hstack([scan.dec, cross.dec]),
                                     0.,
                                     info=np.hstack([scan.info, cross.info]),
                                     dtype=dtype), scan.header

    scan.pa = angle + instrument_angle
    scan.header.update('HIERARCH instrument_angle', instrument_angle)
    return scan


#-------------------------------------------------------------------------------


def _create_scan(ra0, dec0, length, step, sampling_period, speed, acceleration,
                 nlegs, angle, dtype):
    """
    compute the pointing timeline of the instrument reference point
    from the description of a scan map
    Authors: R. Gastaud
    """
    
    length = float(length)
    if length <= 0:
        raise ValueError('Input length must be strictly positive.')

    step = float(step)
    if step <= 0: 
        raise ValueError('Input step must be strictly positive.')

    sampling_period = float(sampling_period)
    if sampling_period <= 0:
        raise ValueError('Input sampling_period must be strictly positive.')
    
    speed = float(speed)
    if speed <= 0:
        raise ValueError('Input speed must be strictly positive.')

    acceleration = float(acceleration)
    if acceleration <= 0:
        raise ValueError('Input acceleration must be strictly positive.')

    nlegs = int(nlegs)
    if nlegs <= 0:
        raise ValueError('Input nlegs must be strictly positive.')
    
    # compute the different times and the total number of points
    # acceleration time at the beginning of a leg, and deceleration time
    # at the end
    # The time needed to turn around is 2 * (extra_time1 + extra_time2)
    extra_time1 = speed / acceleration
    # corresponding length 
    extralength = 0.5 * acceleration * extra_time1 * extra_time1
    # Time needed to go from a scan line to the next 
    extra_time2 = np.sqrt(step / acceleration)

    # Time needed to go along the scanline at constant speed
    line_time = length / speed
    # Total time for a scanline
    full_line_time = extra_time1 + line_time + extra_time1 + extra_time2 + \
                     extra_time2
    # Total duration of the observation
    total_time = full_line_time * nlegs - 2 * extra_time2

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
    delta = -extralength - length/2
    alpha = -step * (nlegs-1)/2
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
            alpha = -step * (nlegs-1) / 2 + line_counter * step
            alpha0 = alpha
   
        # acceleration at the beginning of a scan line to go from 0 to the
        # speed. 
        if working_time < extra_time1:
            delta = -signe * (extralength + length / 2) + signe * 0.5 * \
                    acceleration * working_time * working_time
            info  = Pointing.TURNAROUND
   
        # constant speed
        if working_time >=  extra_time1 and \
           working_time < extra_time1 + line_time:
            delta = signe * (-length / 2 + (working_time - extra_time1) * speed)
            info  = Pointing.INSCAN
 
        # Deceleration at then end of the scanline to stop
        if working_time >= extra_time1 + line_time and \
           working_time < extra_time1 + line_time + extra_time1:
            dt = working_time - extra_time1 - line_time
            delta = signe * (length / 2 + speed * dt - \
                    0.5 * acceleration * dt * dt)
            info  = Pointing.TURNAROUND
  
        # Acceleration to go toward the next scan line
        if working_time >= 2 * extra_time1 + line_time and \
           working_time < 2 * extra_time1 + line_time + extra_time2:
            dt = working_time - 2 * extra_time1 - line_time
            alpha = alpha0 + 0.5 * acceleration * dt * dt
            info  = Pointing.TURNAROUND
   
        # Deceleration to stop at the next scan line
        if working_time >= 2 * extra_time1 + line_time + extra_time2 and \
           working_time < full_line_time:
            dt = working_time - 2 * extra_time1 - line_time - extra_time2
            alpha = (alpha0 + step / 2) + acceleration * extra_time2 * dt - \
                    0.5 * acceleration * dt * dt
            info  = Pointing.TURNAROUND

        time[i] = i * sampling_period
        infos[i] = info
        latitude[i] = delta
        longitude[i] = alpha
        line_counters[i] = line_counter
        working_time = working_time + sampling_period

    # Convert the longitude and latitude *expressed in degrees) to ra and dec
    ra, dec = _change_coord(ra0, dec0, angle - 90, longitude / 3600,
                            latitude / 3600)

    scan = Pointing(time, ra, dec, 0., infos, dtype=dtype)
    header = create_fitsheader(nsamples)
    header.update('ra', ra0)
    header.update('dec', dec0)
    header.update('HIERARCH scan_angle', angle)
    header.update('HIERARCH scan_length', length)
    header.update('HIERARCH scan_nlegs', nlegs)
    header.update('HIERARCH scan_step', step)
    header.update('HIERARCH scan_speed', speed)
    header.update('HIERARCH scan_acceleration', acceleration)
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


#-------------------------------------------------------------------------------


def plot_scan(input, map=None, title=None, new_figure=True, linewidth=2, **kw):

    def plot_scan_(ra, dec, nsamples, image, linewidth=linewidth, **kw):
        x, y = image.topixel(ra, dec)
        dest = 0
        for n in nsamples:
            p = pyplot.plot(x[dest:dest+n], y[dest:dest+n],
                            linewidth=linewidth, **kw)
            for i in xrange(dest, dest+n):
                if np.isfinite(x[i]) and np.isfinite(y[i]):
                    pyplot.plot(x[i], y[i], 'o', color=p[0]._color)
                    break
            dest += n

    if type(input) not in (list, tuple):
        input = [input]
    else:
        input = list(input)

    # format the input as a list of (ra, dec, nsamples)
    for islice, slice in enumerate(input):
        if hasattr(slice, 'pointing'):
            slice = slice.pointing
        if hasattr(slice, 'ra') and hasattr(slice, 'dec'):
            ra = slice.ra
            dec = slice.dec
            if hasattr(slice, 'removed'):
                ra = ra.copy()
                dec = dec.copy()
                ra [slice.removed] = np.nan
                dec[slice.removed] = np.nan
            nsamples = getattr(slice, 'nsamples', (ra.size,))
            input[islice] = (ra, dec, nsamples)
        elif type(slice) not in (list,tuple) or len(slice) != 2:
            raise TypeError('Invalid input.')
        else:
            ra = np.array(slice[0], ndmin=1, copy=False)
            dec = np.array(slice[1], ndmin=1, copy=False)
            if ra.shape != dec.shape:
                raise ValueError('Input R.A. and Dec do not have the same len' \
                                 'gth.')
            input[islice] = (ra, dec, (ra.size,))

    npointings = sum([s[0].size for s in input])
    ra = np.empty(npointings)
    dec = np.empty(npointings)

    dest = 0
    nsamples = []
    for slice in input:
        n = slice[0].size
        ra [dest:dest+n] = slice[0]
        dec[dest:dest+n] = slice[1]
        nsamples += slice[2]
        dest += n

    nvalids = np.sum(np.isfinite(ra*dec))
    if nvalids == 0:
        raise ValueError('There is no valid coordinates.')

    if isinstance(map, Map) and map.has_wcs():
        image = map.imshow(title=title)
        plot_scan_(ra, dec, nsamples, image, linewidth=linewidth, **kw)
        return image

    crval = barycenter_lonlat(ra, dec)
    angles = angle_lonlat((ra,dec), crval)
    angle_max = np.nanmax(angles)
    if angle_max >= 90.:
        print('Warning: some coordinates have an angular distance to the proj' \
              'ection point greater than 90 degrees.')
        mask = angles >= 90
        ra[mask] = np.nan
        dec[mask] = np.nan
        angles[mask] = np.nan
        angle_max = np.nanmax(angles)
    cdelt = np.rad2deg(np.tan(np.deg2rad(angle_max))) / 100
    header = create_fitsheader(naxis=[201,201], cdelt=cdelt, crval=crval,
                               crpix=[101,101])

    fitsobj = kapteyn.maputils.FITSimage(externalheader=header)
    if new_figure:
        fig = pyplot.figure()
        frame = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    else:
        frame = pyplot.gca()
    if title is not None:
        frame.set_title(title)
    image = fitsobj.Annotatedimage(frame, blankcolor='w')
    image.Graticule()
    image.plot()
    image.interact_toolbarinfo()
    image.interact_writepos()
    pyplot.show()
    plot_scan_(ra, dec, nsamples, image, linewidth=linewidth, **kw)
    return image
