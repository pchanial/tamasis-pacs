# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
from __future__ import division

import kapteyn.maputils
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import pickle
import pyfits
import StringIO

try:
    import ds9
    _imported_ds9 = True
except:
    _imported_ds9 = False

from functools import reduce

from . import MPI
from .wcsutils import create_fitsheader
from .mpiutils import read_fits, write_fits
from .quantity import Quantity

__all__ = [ 'FitsArray', 'Map', 'Tod' ]


class FitsArray(Quantity):

    _header = None

    def __new__(cls, data, header=None, unit=None, derived_units=None,
                dtype=None, copy=True, order='C', subok=False, ndmin=0,
                comm=None):

        if isinstance(data, (str, unicode)):

            if comm is None:
                comm = MPI.COMM_WORLD

            try:
                derived_units = pyfits.open(data)['derived_units'].data
                derived_units = pickle.loads(str(derived_units.data))
            except KeyError:
                pass

            data, header = read_fits(data, None, comm)

            copy = False
            if unit is None:
                if 'BUNIT' in header:
                    unit = header['BUNIT']
                elif 'QTTY____' in header:
                    unit = header['QTTY____'] # HCSS crap
                    if unit == '1':
                        unit = ''

        elif comm is not None:
            raise ValueError('The MPI communicator can only be set for input FI'
                             'TS files.')

        # get a new FitsArray instance (or a subclass if subok is True)
        result = Quantity.__new__(cls, data, unit, derived_units, dtype, copy,
                                  order, True, ndmin)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)

        # copy header attribute
        if header is not None:
            result.header = header
        elif hasattr(data, '_header') and \
             data._header.__class__ is pyfits.Header:
            if copy:
                result._header = data._header.copy()
            else:
                result._header = data._header
        else:
            result._header = None

        return result

    def __array_finalize__(self, array):
        Quantity.__array_finalize__(self, array)
        self._header = getattr(array, '_header', None)

    def __getattr__(self, name):
        if self.dtype.names and name in self.dtype.names:
            return self[name]
        return super(FitsArray, self).__getattribute__(name)

    def __setattr__(self, name, value):
        if self.dtype.names and name in self.dtype.names:
            self[name] = value
        else:
            super(FitsArray,self).__setattr__(name, value)

    @staticmethod
    def empty(shape, header=None, unit=None, derived_units=None, dtype=None,
              order=None):
        return FitsArray(np.empty(shape, dtype, order), header, unit,
            derived_units, dtype, copy=False)

    @staticmethod
    def ones(shape, header=None, unit=None, derived_units=None, dtype=None,
             order=None):
        return FitsArray(np.ones(shape, dtype, order), header, unit,
            derived_units, dtype, copy=False)

    @staticmethod
    def zeros(shape, header=None, unit=None, derived_units=None, dtype=None,
              order=None):
        return FitsArray(np.zeros(shape, dtype, order), header, unit,
            derived_units, dtype, copy=False)

    def has_wcs(self):
        """
        Returns True is the array has a FITS header with a defined World
        Coordinate System.
        """
        if self.header is None:
            return False

        required = 'CRPIX,CRVAL,CTYPE'.split(',')
        keywords = np.concatenate(
            [(lambda i: [r+str(i+1) for r in required])(i) 
             for i in range(self.header['NAXIS'])])

        return all([k in self.header for k in keywords])

    @property
    def header(self):
        if self._header is None and not np.iscomplexobj(self):
            self._header = create_fitsheader(self.shape[::-1], fromdata=self)
        return self._header

    @header.setter
    def header(self, header):
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header (' + \
                            str(type(header))+').')
        self._header = header

    def tofile(self, fid, sep='', format='%s'):
        super(FitsArray,self).tofile(fid, sep, format)

    def save(self, filename, header=None, comm=None):
        """
        Write a FitsArray instance as a FITS file.
        """

        header = (header or self.header).copy()
        if comm is None:
            comm = MPI.COMM_WORLD

        if len(self._unit) != 0:
            header.update('BUNIT', self.unit)

        write_fits(filename, self, header, False, None, comm)
        if not self.derived_units:
            return
        if comm is None or comm.rank == 0:
            buf = StringIO.StringIO()
            pickle.dump(self.derived_units, buf, pickle.HIGHEST_PROTOCOL)
            data = np.frombuffer(buf.getvalue(), np.uint8)
            write_fits(filename, data, None, True, 'derived_units',
                       MPI.COMM_SELF)
        comm.Barrier()

    def imsave(self, filename, colorbar=True, **kw):
        is_interactive = matplotlib.is_interactive()
        matplotlib.interactive(False)
        dpi = 80.
        figsize = np.clip(np.max(np.array(self.shape[::-1])/dpi), 8, 50)
        figsize = (figsize + (2 if colorbar else 0), figsize)
        fig = pyplot.figure(figsize=figsize, dpi=dpi)
        self.imshow(colorbar=colorbar, new_figure=False, **kw)
        fig = pyplot.gcf()
        fig.savefig(filename)
        matplotlib.interactive(is_interactive)

    def imshow(self, mask=None, title=None, colorbar=True, aspect=None,
               interpolation='nearest', origin=None, xlabel='', ylabel='',
               new_figure=True, **kw):
        """
        A simple graphical display function for the Tod class

        mask: array-like
            True means masked.
        """

        unfinite = ~np.isfinite(np.asarray(self))
        if mask is None:
            mask = unfinite
        else:
            mask = np.logical_or(mask, unfinite)

        data = np.ma.MaskedArray(np.asarray(self), mask=mask, copy=False)
        if np.iscomplexobj(data):
            data = abs(data)
        mean   = np.mean(data)
        stddev = np.std(data)
        # casting to float because of a bug numpy1.4 + matplotlib
        minval = float(max(mean - 2*stddev, np.min(data)))
        maxval = float(min(mean + 5*stddev, np.max(data)))

        if new_figure:
            fig = pyplot.figure()
        else:
            fig = pyplot.gcf()
        fontsize = 12. * fig.get_figheight() / 6.125

        image = pyplot.imshow(data, aspect=aspect, interpolation=interpolation,
                              origin=origin, **kw)
        image.set_clim(minval, maxval)

        ax = pyplot.gca()
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)

        if title is not None:
            pyplot.title(title, fontsize=fontsize)
        if colorbar:
            colorbar = pyplot.colorbar()
            for tick in colorbar.ax.get_yticklabels():
                tick.set_fontsize(fontsize)

        pyplot.draw()
        return image

    def ds9(self, xpamsg=None, origin=None, new=True, **keywords):
        """
        Display the array using ds9.

        The ds9 process will be given an random id. By default, the
        following access point are set before the array is loaded:
            -cmap heat
            -scale scope local
            -scale mode 99.5
        Other access points can be set before the data is loaded though
        the keywords (see examples below).
        After the array is loaded, the map's header is set and the user
        may add other XPA messages through the xpamsg argument or by
        setting them through the returned ds9 instance.

        Parameters
        ----------
        xpamsg : string or tuple of string
            XPA access point message to be set after the array is loaded.
            (see http://hea-www.harvard.edu/RD/ds9/ref/xpa.html).
        origin: string
            Set origin to 'upper' for Y increasing downwards
        new: boolean
            If true, open the array in a new ds9 instance.
        **keywords : string or tuple of string
            Specify more access points to be set before array loading.
            a keyword such as 'height=400' will be appended to the command
            that launches ds9 in the form 'ds9 [...] -height 400'

        Returns
        -------
        The returned object is a ds9 instance. It can be manipulated using
        XPA access points.

        Examples
        --------
        >>> m = Map('myfits.fits')
        >>> d=m.ds9(('zoom to fit','saveimage png myfits.png'),scale='histequ', 
                    cmap='invert yes', height=400)
        >>> d.set('exit')
        """
        if not _imported_ds9:
            raise RuntimeError('The library pyds9 has not been installed.')
        import os, time, sys, uuid, xpa

        id = None
        if not new:
            list = ds9.ds9_targets()
            if list is not None:
                id = list[-1]

        if id is None:
            if 'cmap' not in keywords:
                keywords['cmap'] = 'heat'

            if 'scale' not in keywords:
                keywords['scale'] = ('scope local', 'mode 99.5')

            if origin == 'upper' or \
               'orient' not in keywords and self.origin == 'upper':
                keywords['orient'] = 'y'

            wait = 10

            id = 'ds9_' + str(uuid.uuid1())[4:8]

            command = 'ds9 -title ' + id

            for k, v in keywords.items():
                k = str(k)
                if type(v) is not tuple:
                    v = (v,)
                command += reduce(lambda x,y: \
                                  str(x) + ' -' + k + ' ' + str(y),v,'')

            os.system(command + ' &')

            # start the xpans name server
            if xpa.xpaaccess("xpans", None, 1) == None:
                _cmd = None
                # look in install directories
                for _dir in sys.path:
                    _fname = os.path.join(_dir, "xpans")
                    if os.path.exists(_fname):
                        _cmd = _fname + " -e &"
                if _cmd:
                    os.system(_cmd)

            for i in range(wait):
                list = xpa.xpaaccess(id, None, 1024)
                if list: break
                time.sleep(1)
            if not list:
                raise ValueError('No active ds9 running for target: %s' % list)

        # get ds9 instance with given id
        d = ds9.ds9(id)

        # load array
        input = self.view(np.ndarray)
        if input.dtype.kind in ('b', 'i'):
            input = np.array(input, np.int32, copy=False)
        d.set_np2arr(input.T)

        # load header
        if self.has_wcs():
            d.set('wcs append', str(self.header))

        if xpamsg is not None:
            if isinstance(xpamsg, str):
                xpamsg = (xpamsg,)
            for v in xpamsg:
                d.set(v)

        return d


class Map(FitsArray):

    """
    Represent a map, complemented with unit and FITS header.
    """

    coverage = None
    error = None
    origin = None

    def __new__(cls, data,  header=None, unit=None, derived_units=None,
                coverage=None, error=None, origin=None, dtype=None, copy=True,
                order='C', subok=False, ndmin=0, comm=None):

        # get a new Map instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units,
            dtype, copy, order, True, ndmin, comm)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)

        if isinstance(data, (str, unicode)):

            if comm is None:
                comm = MPI.COMM_WORLD

            if 'DISPORIG' in result.header:
                if origin is None:
                    origin = result.header['DISPORIG']
                del result.header['DISPORIG']
            if coverage is None:
                try:
                    coverage, junk = read_fits(data, 'Coverage', comm)
                except KeyError:
                    pass
            if error is None:
                try:
                    error, junk = read_fits(data, 'Error', comm)
                except KeyError:
                    pass

        if origin is not None:
            origin = origin.strip().lower()
            if origin not in ('upper', 'lower', 'none'):
                raise ValueError("Invalid origin '" + origin + "'. Expected v" \
                                 "alues are None, 'upper' or 'lower'.")
            if origin != 'none':
                result.origin = origin

        if coverage is not None:
            result.coverage = coverage
        elif copy and result.coverage is not None:
            result.coverage = result.coverage.copy()

        if error is not None:
            result.error = error
        elif copy and result.error is not None:
            result.error = result.error.copy()

        return result

    def __array_finalize__(self, array):
        FitsArray.__array_finalize__(self, array)
        self.coverage = getattr(array, 'coverage', None)
        self.error = getattr(array, 'error', None)
        self.origin = getattr(array, 'origin', 'lower')

    def __getitem__(self, key):
        item = super(Map, self).__getitem__(key)
        if not isinstance(item, Map):
            return item
        if item.coverage is not None:
            item.coverage = item.coverage[key]
        if item.error is not None:
            item.error = item.error[key]
        return item

    @staticmethod
    def empty(shape, coverage=None, error=None, origin='lower', header=None,
              unit=None, derived_units=None, dtype=None, order=None):
        return Map(np.empty(shape, dtype, order), header, unit, derived_units,
                   coverage, error, origin, dtype, copy=False)

    @staticmethod
    def ones(shape, coverage=None, error=None, origin='lower', header=None,
             unit=None, derived_units=None, dtype=None, order=None):
        return Map(np.ones(shape, dtype, order), header, unit, derived_units,
                   coverage, error, origin, dtype, copy=False)

    @staticmethod
    def zeros(shape, coverage=None, error=None, origin='lower', header=None,
              unit=None, derived_units=None, dtype=None, order=None):
        return Map(np.zeros(shape, dtype, order), header, unit, derived_units,
                   coverage, error, origin, dtype, copy=False)

    def imshow(self, mask=None, title=None, new_figure=True, origin=None, **kw):
        """A simple graphical display function for the Map class"""

        if mask is None and self.coverage is not None:
            mask = self.coverage <= 0
        if mask is not None:
            data = np.array(self, copy=True)
            data[mask] = np.nan
        else:
            data = np.asarray(self)

        if origin is None:
            origin = self.origin

        # check if the map has no astrometry information
        if not self.has_wcs():
            if 'xlabel' not in kw:
                kw['xlabel'] = 'X'
            if 'ylabel' not in kw:
                kw['ylabel'] = 'Y'
            image = super(Map, self).imshow(title=title, new_figure=new_figure,
                                            origin=origin, **kw)
            return image

        fitsobj = kapteyn.maputils.FITSimage(externaldata=data,
                                             externalheader=self.header)
        if new_figure:
            fig = pyplot.figure()
            frame = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        else:
            frame = pyplot.gca()
        if title is not None:
            frame.set_title(title)
        annim = fitsobj.Annotatedimage(frame, blankcolor='w')
        annim.Image(interpolation='nearest')
        grat = annim.Graticule()
        grat.setp_gratline(visible=False)
        annim.plot()
        annim.interact_imagecolors()
        annim.interact_toolbarinfo()
        annim.interact_writepos()
        pyplot.show()
        return annim

    def save(self, filename, header=None, comm=None):

        header = (header or self.header).copy()
        if comm is None:
            comm = MPI.COMM_WORLD

        header.update('DISPORIG', self.origin, 'Map display convention')
        FitsArray.save(self, filename, header, comm=comm)
        if self.coverage is not None:
            write_fits(filename, self.coverage, None, True, 'Coverage', comm)
        if self.error is not None:
            write_fits(filename, self.error, None, True, 'Error', comm)


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    _mask = None

    def __new__(cls, data, mask=None, header=None, unit=None,
                derived_units=None, dtype=None, copy=True, order='C',
                subok=False, ndmin=0, comm=None):

        # get a new Tod instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units,
            dtype, copy, order, True, ndmin, comm)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)
        
        # mask attribute
        if mask is np.ma.nomask:
            mask = None

        if mask is None and isinstance(data, (str, unicode)):

            if comm is None:
                comm = MPI.COMM_WORLD

            try:
                mask, junk = read_fits(data, 'mask', comm)
                mask = mask.view(np.bool8)
                copy = False
            except:
                pass

        if mask is None and hasattr(data, 'mask') and \
           data.mask is not np.ma.nomask:
            mask = data.mask

        if mask is not None:
            result._mask = np.array(mask, np.bool8, copy=copy)

        return result

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        if mask is None or mask is np.ma.nomask:
            self._mask = None
            return

        # enforce bool8 dtype
        if not isinstance(mask, np.ndarray):
            mask = np.array(mask, np.bool8)
        elif mask.dtype.type != np.bool8:
            if mask.dtype.itemsize == 1:
                mask = mask.view(np.bool8)
            else:
                mask = np.asarray(mask, np.bool8)

        # handle the scalar case
        if np.rank(mask) == 0:
            if self._mask is None:
                func = np.zeros if mask == 0 else np.ones
                self._mask = func(self.shape, dtype=np.bool8)
            else:
                self._mask[:] = mask
            return

        # check shape compatibility
        if self.shape != mask.shape:
            raise ValueError("The input mask has a shape '" + str(mask.shape) +\
                "' incompatible with that of the Tod '" + str(self.shape) +"'.")
        
        self._mask = mask

    def __array_finalize__(self, array):
        FitsArray.__array_finalize__(self, array)
        self._mask = getattr(array, 'mask', None)

    def __getitem__(self, key):
        item = super(Tod, self).__getitem__(key)
        if not isinstance(item, Tod):
            return item
        if item.mask is not None:
            item.mask = item.mask[key]
        return item

    def reshape(self, newdims, order='C'):
        result = np.ndarray.reshape(self, newdims, order=order)
        if self.mask is not None:
            result.mask = self.mask.reshape(newdims, order=order)
        result.derived_units = self.derived_units.copy()
        return result

    @staticmethod
    def empty(shape, mask=None, header=None, unit=None, derived_units=None,
              dtype=None, order=None):
        return Tod(np.empty(shape, dtype, order), mask, header, unit,
                   derived_units, dtype, copy=False)

    @staticmethod
    def ones(shape, mask=None, header=None, unit=None, derived_units=None,
             dtype=None, order=None):
        return Tod(np.ones(shape, dtype, order), mask, header, unit,
                   derived_units, dtype, copy=False)

    @staticmethod
    def zeros(shape, mask=None, header=None, unit=None, derived_units=None,
              dtype=None, order=None):
        return Tod(np.zeros(shape, dtype, order), mask, header, unit,
                   derived_units, dtype, copy=False)

    def imshow(self, title=None, aspect='auto', **kw):
        """
        A simple graphical display function for the Map class
        """

        xlabel = 'Sample'
        ylabel = 'Detector number'
        image = super(Tod, self).imshow(mask=self.mask, title=title,
                                        origin='upper', xlabel=xlabel,
                                        ylabel=ylabel, aspect=aspect, **kw)
        return image

    def __str__(self):
        if np.rank(self) == 0:
            if np.iscomplexobj(self):
                return 'Tod ' + str(complex(self))
            else:
                return 'Tod ' + str(float(self))
        output = 'Tod ['
        if np.rank(self) > 1:
            output += str(self.shape[-2])+' detector'
            if self.shape[-2] > 1:
                output += 's'
            output += ', '
        output += str(self.shape[-1]) + ' sample'
        if self.shape[-1] > 1:
            output += 's'
        return output + ']'

    def save(self, filename, header=None, comm=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        FitsArray.save(self, filename, header, comm=comm)
        if self.mask is not None:
            mask = self.mask.view('uint8')
            write_fits(filename, mask, None, True, 'Mask', comm)
