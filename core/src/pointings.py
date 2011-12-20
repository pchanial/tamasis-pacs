from __future__ import division

import kapteyn
import numpy as np
from matplotlib import pyplot
from pyoperators.utils import isscalar

from . import var
from .datatypes import FitsArray, Map
from .quantity import Quantity
from .wcsutils import angle_lonlat, barycenter_lonlat, create_fitsheader

__all__ = ['Pointing']
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

    def get_map_header(self, resolution=None, naxis=None):
        """
        Returns a FITS header that encompasses all non-removed and
        non-masked pointings.
        """
        mask = ~self.masked & ~self.removed
        ra = self.ra[mask]
        dec = self.dec[mask]
        crval = barycenter_lonlat(ra, dec)
        angles = angle_lonlat((ra, dec), crval)
        angle_max = np.nanmax(angles)
        if angle_max >= 90.:
            print('Warning: some coordinates have an angular distance to the pr'
                  'ojection point greater than 90 degrees.')
            angles[angles >= 90] = np.nan
            angle_max = np.nanmax(angles)
        angle_max = np.deg2rad(angle_max)
        if resolution is None:
            cdelt = 2 * np.rad2deg(np.tan(angle_max)) / naxis
        else:
            cdelt = resolution / 3600
            naxis = 2 * np.ceil(np.rad2deg(np.tan(angle_max)) / cdelt)
        return create_fitsheader(naxis=2*[naxis], cdelt=cdelt, crval=crval,
                                 crpix=2*[naxis//2+1])

    def plot(self, map=None, header=None, title=None, new_figure=True,
             linewidth=2, **kw):
        """
        Plot the pointings' celestial coordinates.

        Parameters
        ----------
        map : ndarray of dim 2
            An optional map to be displayed as background.
        header : pyfits.Header
            The optional map's FITS header.
        title : str
            The figure's title.
        new_figure : boolean
            If true, the display will be done in a new window.
        linewidth : float
            The scan line width.
        """

        mask = ~self.masked & ~self.removed
        ra = self.ra[mask]
        dec = self.dec[mask]
        if ra.size == 0 and new_figure:
            raise ValueError('There is no valid coordinates.')

        if header is None:
            header = getattr(map, 'header', None)
        elif map is not None:
            map = Map(map, header=header, copy=False)

        if isinstance(map, Map) and map.has_wcs():
            image = map.imshow(title=title, new_figure=new_figure)
            _plot_scan(image, ra, dec, linewidth=linewidth, **kw)
            return image

        if header is None:
            header = self.get_map_header(naxis=1)
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
        _plot_scan(image, ra, dec, linewidth=linewidth, **kw)
        return image

def _plot_scan(image, ra, dec, linewidth=None, **kw):
    x, y = image.topixel(ra, dec)
    p = pyplot.plot(x, y, linewidth=linewidth, **kw)
    for i in xrange(ra.size):
        if np.isfinite(x[i]) and np.isfinite(y[i]):
            pyplot.plot(x[i], y[i], 'o', color=p[0]._color)
            break
