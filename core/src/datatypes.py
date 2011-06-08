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

import tamasisfortran as tmf

from functools import reduce
from mpi4py import MPI

from . import var
from .numpyutils import _my_isscalar
from .wcsutils import create_fitsheader
from .mpiutils import read_fits, write_fits, split_shape, split_work
from .quantity import Quantity, UnitError, _extract_unit, _strunit

__all__ = [ 'FitsArray', 'Map', 'Tod' ]


class DistributedArray(np.ndarray):
    """ndarray subclass, to handle MPI-distributed arrays.
    """
    def __new__(cls, data, shape_global=None, comm=MPI.COMM_SELF):
        if shape_global is None and comm.Get_size() > 1:
            raise ValueError('The global shape of the local array is not speci'\
                             'fied.')
        if not isinstance(data, DistributedArray):
            data = data.view(cls)
        data.shape_global = shape_global or data.shape
        data.comm = comm
        return data

    def __array_finalize__(self, array):
        if array is None:
            return
        self.shape_global = getattr(array, 'shape_global', array.shape)
        self.comm = getattr(array, 'comm', MPI.COMM_SELF)

    def __reduce__(self):
        state = list(np.ndarray.__reduce__(self))
        subclass_state = self.__dict__.copy()
        subclass_state['comm'] = self.comm.py2f()
        state[2] = (state[2],subclass_state)
        return tuple(state)
    
    def __setstate__(self,state):
        ndarray_state, subclass_state = state
        np.ndarray.__setstate__(self,ndarray_state)        
        for k, v in subclass_state.items():
            if k == 'comm':
                v = MPI.Comm.f2py(v)
            setattr(self, k, v)

    def toglobal(self, attr=()):
        """Gather the local images into a global image."""
        if self.comm is None:
            raise RuntimeError('This array is not a local image.')
        counts = []
        offsets = [0]
        for rank in range(self.comm.Get_size()):
            s = split_work(self.shape_global[0], rank=rank, comm=self.comm)
            n = (s.stop - s.start) * np.product(self.shape_global[1:])
            counts.append(n)
            offsets.append(offsets[-1] + n)
        offsets.pop()
        s = split_work(self.shape_global[0], comm=self.comm)
        n = s.stop - s.start
        output = self.empty(self.shape_global, dtype=self.dtype,
                            comm=MPI.COMM_SELF)
        output.__array_finalize__(self)
        t = MPI.BYTE.Create_contiguous(self.dtype.itemsize)
        t.Commit()
        self.comm.Allgatherv([self[0:n], t], [output.view(np.byte), (counts, offsets), t])

        for a in attr:
            i = getattr(self, a, None)
            if i is None:
                continue
            o = np.empty(self.shape_global, dtype=i.dtype)
            t = MPI.BYTE.Create_contiguous(i.dtype.itemsize)
            t.Commit()
            self.comm.Allgatherv([i[0:n], t], [o, (counts, offsets), t])
            setattr(output, a, o)
            
        output.comm = MPI.COMM_SELF

        return output

    def tolocal(self, comm=MPI.COMM_WORLD):
        """Scatter a global image into local ones."""
        if hasattr(self, 'shape_global') and self.shape != self.shape_global or\
           hasattr(self, 'comm') and self.comm.Get_size() > 1:
            raise ValueError('This array is not a global image.')

        order = 'f' if np.isfortran(self) else 'c'
        if hasattr(self, 'empty'):
            output = self.empty(self.shape, dtype=self.dtype, order=order,
                                comm=comm)
        else:
            shape = split_shape(self.shape, comm)
            output = np.empty(shape, dtype=self.dtype, order=order)
            output = DistributedArray.__new__(DistributedArray, output,
                                              self.shape, comm)
        if hasattr(self, '__dict__'):
            for k,v in self.__dict__.items():
                setattr(output, k, v)
        output.comm = comm

        s = split_work(self.shape[0], comm=comm)
        n = s.stop - s.start
        output[0:n] = self[s.start:s.stop]
        if n < output.shape[0]:
            output[n:] = 0
        return output


class FitsArray(DistributedArray, Quantity):

    def __new__(cls, data, header=None, unit=None, derived_units=None,
                dtype=None, copy=True, order='C', subok=False, ndmin=0,
                shape_global=None, comm=None):

        comm = comm or getattr(data, 'comm', MPI.COMM_SELF)
        shape_global = shape_global or getattr(data, 'shape_global', None)

        if type(data) is str:
            ihdu = 0
            fits = pyfits.open(data)
            while True:
                try:
                    hdu = fits[ihdu]
                except IndexError:
                    raise IOError('The FITS file has no data.')
                if hdu.header['NAXIS'] == 0:
                    ihdu += 1
                    continue
                if hdu.data is not None:
                    break
            if comm.Get_size() > 1:
                files = comm.allgather(data)
                # before splitting the array into local image, check that the
                # file name is the same for all MPI process
                if any([f != files[0] for f in files]):
                    raise ValueError('The file name is not the same for all MP'\
                                     'I processses.')
            data, header, shape_global = read_fits(hdu, comm)

            copy = False
            if unit is None:
                if 'BUNIT' in header:
                    unit = header['BUNIT']
                elif 'QTTY____' in header:
                    unit = header['QTTY____'] # HCSS crap
            try:
                derived_units = fits['derived_units'].data
                derived_units = pickle.loads(str(derived_units.data))
                comm.Barrier()
            except KeyError:
                pass

        # get a new FitsArray instance (or a subclass if subok is True)
        result = Quantity.__new__(cls, data, unit, derived_units, dtype, copy,
                                  order, True, ndmin)
        result = DistributedArray.__new__(cls, result, shape_global, comm)
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
        DistributedArray.__array_finalize__(self, array)
        self._header = getattr(array, '_header', None)

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if not isinstance(item, self.__class__):
            return item
        item.shape_global = item.shape
        item.comm = MPI.COMM_SELF
        return item

    def __getattr__(self, name):
        if self.dtype.names is None or name not in self.dtype.names:
            raise AttributeError("'" + self.__class__.__name__ + "' object ha" \
                "s no attribute '" + name + "'")
        return self[name]

    def __setattr__(self, name, value):
        if self.dtype.names and name in self.dtype.names:
            self[name] = value
        else:
            super(FitsArray,self).__setattr__(name, value)

    @staticmethod
    def empty(shape, header=None, unit=None, derived_units=None, dtype=None,
              order=None, comm=MPI.COMM_SELF):
        shape_local = split_shape(shape, comm)
        return FitsArray(np.empty(shape_local, dtype, order), header, unit,
            derived_units, dtype, copy=False, shape_global=shape, comm=comm)

    @staticmethod
    def ones(shape, header=None, unit=None, derived_units=None, dtype=None,
             order=None, comm=MPI.COMM_SELF):
        shape_local = split_shape(shape, comm)
        return FitsArray(np.ones(shape_local, dtype, order), header, unit,
            derived_units, dtype, copy=False, shape_global=shape, comm=comm)

    @staticmethod
    def zeros(shape, header=None, unit=None, derived_units=None, dtype=None,
              order=None, comm=MPI.COMM_SELF):
        shape_local = split_shape(shape, comm)
        return FitsArray(np.zeros(shape_local, dtype, order), header, unit,
            derived_units, dtype, copy=False, shape_global=shape, comm=comm)

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
            self._header = create_fitsheader(fromdata=self)
        return self._header

    @header.setter
    def header(self, header):
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header (' + \
                            str(type(header))+').')
        self._header = header

    def tofile(self, fid, sep='', format='%s'):
        super(FitsArray,self).tofile(fid, sep, format)

    def save(self, filename, fitskw=None):
        """Save a FitsArray instance to a fits file given a filename
       
        If the same file already exist it overwrites it.
        """

        if self.header is not None:
            header = self.header.copy()
        else:
            header = create_fitsheader(fromdata=self)
       
        if len(self._unit) != 0:
            header.update('BUNIT', self.unit)

        if fitskw is not None:
            for k, v in fitskw.items():
                if len(k) > 8:
                    k = 'HIERARCH ' + k
                if isinstance(v, tuple):
                    if len(v) != 2:
                        raise ValueError('A tuple input for fitskw must have t'\
                                         'wo elements (value, comment).')
                    v, c = v
                else:
                    c = None
                v = v if isinstance(v, (str, int, long, float, complex, bool,
                    np.floating, np.integer, np.complexfloating)) else repr(v)
                header.update(k, v, c)

        if np.rank(self) == 0:
            value = self.reshape((1,))
        else:
            value = self.T if np.isfortran(self) else self
        write_fits(filename, self, header, self.shape_global, False, self.comm)
        if not isinstance(self, Map) and not isinstance(self, Tod):
            _save_derived_units(filename, self.derived_units, self.comm)

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
        import ds9, os, time, sys, uuid, xpa

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

            for k,v in keywords.items():
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
    def __new__(cls, data,  header=None, unit=None, derived_units=None,
                coverage=None, error=None, origin=None, dtype=None, copy=True,
                order='C', subok=False, ndmin=0, shape_global=None, comm=None):

        # get a new Map instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units,
            dtype, copy, order, True, ndmin, shape_global, comm)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)

        if type(data) is str:
            if 'DISPORIG' in result.header:
                if origin is None:
                    origin = result.header['DISPORIG']
                del result.header['DISPORIG']
            fits = pyfits.open(data)
            try:
                if coverage is None:
                    coverage, j1, j2 = read_fits(fits['Coverage'], comm)
            except:
                pass
            try:
                if error is None:
                    error, j1, j2 = read_fits(fits['Error'], comm)
            except:
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
              unit=None, derived_units=None, dtype=None, order=None,
              comm=MPI.COMM_SELF):
        shape_local = split_shape(shape, comm)
        return Map(np.empty(shape_local, dtype, order), header, unit,
                   derived_units, coverage, error, origin, dtype, copy=False,
                   shape_global=shape, comm=comm)

    @staticmethod
    def ones(shape, coverage=None, error=None, origin='lower', header=None,
             unit=None, derived_units=None, dtype=None, order=None,
             comm=MPI.COMM_SELF):
        shape_local = split_shape(shape, comm)
        return Map(np.ones(shape_local, dtype, order), header, unit,
                   derived_units, coverage, error, origin, dtype, copy=False,
                   shape_global=shape, comm=comm)

    @staticmethod
    def zeros(shape, coverage=None, error=None, origin='lower', header=None,
              unit=None, derived_units=None, dtype=None, order=None,
              comm=MPI.COMM_SELF):
        shape_local = split_shape(shape, comm)
        return Map(np.zeros(shape_local, dtype, order), header, unit,
                   derived_units, coverage, error, origin, dtype, copy=False,
                   shape_global=shape, comm=comm)

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

    def save(self, filename, fitskw=None):
        fitskw = fitskw or {}
        if 'DISPORIG' not in [k.upper() for k in fitskw.keys()]:
            fitskw['DISPORIG'] = (self.origin, 'Map display convention')
        FitsArray.save(self, filename, fitskw=fitskw)
        if self.coverage is not None:
            write_fits(filename, self.coverage, None, self.shape_global, True,
                       self.comm, extname='Coverage')
        if self.error is not None:
            write_fits(filename, self.error, None, self.shape_global, True,
                       self.comm, extname='Error')
        _save_derived_units(filename, self.derived_units, self.comm)

    def toglobal(self):
        return DistributedArray.toglobal(self, attr=('coverage', 'error'))

    def tolocal(self, comm=MPI.COMM_WORLD):
        output = DistributedArray.tolocal(self, comm)
        tolocal = DistributedArray.__dict__['tolocal']
        if self.coverage is not None:
            output.coverage = tolocal(self.coverage, comm)
        if self.error is not None:
            output.error = tolocal(self.error, comm)
        return output


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None,
                derived_units=None, dtype=None, copy=True, order='C',
                subok=False, ndmin=0, shape_global=None, comm=None):

        # get a new Tod instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units,
            dtype, copy, order, True, ndmin, shape_global, comm)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)
        
        # mask attribute
        if mask is np.ma.nomask:
            mask = None

        if mask is None and isinstance(data, str):
            try:
                mask, j, k = read_fits(pyfits.open(data)['Mask'], comm)
                mask = mask.view(np.bool8)
                copy = False
            except:
                pass

        if mask is None and hasattr(data, 'mask') and \
           data.mask is not np.ma.nomask:
            mask = data.mask

        if mask is not None:
            result._mask = np.array(mask, np.bool8, copy=copy)
        
        # nsamples attribute
        if type(data) is str and nsamples is None:
            if 'nsamples' in result.header:
                nsamples = result.header['nsamples'][1:-1].replace(' ', '')
                if len(nsamples) > 0:
                    nsamples = [int(float(x))
                                for x in nsamples.split(',') if x.strip() != '']
                del result.header['NSAMPLES']
        if nsamples is None:
            return result
        shape = validate_sliced_shape(result.shape, nsamples)
        result.nsamples = shape[-1] if len(shape) > 0 else shape

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
        self.nsamples = getattr(array, 'nsamples', () if np.rank(self) == 0 \
                        else (self.shape[-1],))

    def __getitem__(self, key):
        item = super(Tod, self).__getitem__(key)
        if not isinstance(item, Tod):
            return item
        if item.mask is not None:
            item.mask = item.mask[key]
        if not isinstance(key, tuple):
            return item
        if len(key) > 1:
            if not isinstance(key[-1], slice):
                return item
            else:
                if key[-1].start is not None or key[-1].stop is not None or \
                   key[-1].step is not None:
                    item.nsamples = (item.shape[-1],)
        return item

    def reshape(self, newdims, order='C'):
        result = np.ndarray.reshape(self, newdims, order=order)
        if self.mask is not None:
            result.mask = self.mask.reshape(newdims, order=order)
        result.derived_units = self.derived_units.copy()
        return result

    @staticmethod
    def empty(shape, mask=None, nsamples=None, header=None, unit=None,
              derived_units=None, dtype=None, order=None, comm=MPI.COMM_SELF):
        nsamples = nsamples or shape[-1]
        shape = flatten_sliced_shape(shape)
        if len(shape) <= 1 and comm.Get_size() > 1:
            raise ValueError('Scalar or vector Tod can not be distributed.')
        shape_local = split_shape(shape, comm)
        return Tod(np.empty(shape_local, dtype, order), mask, nsamples, header,
                   unit, derived_units, dtype, copy=False, shape_global=shape,
                   comm=comm)

    @staticmethod
    def ones(shape, mask=None, nsamples=None, header=None, unit=None,
             derived_units=None, dtype=None, order=None, comm=MPI.COMM_SELF):
        nsamples = nsamples or shape[-1]
        shape = flatten_sliced_shape(shape)
        if len(shape) <= 1 and comm.Get_size() > 1:
            raise ValueError('Scalar or vector Tod can not be distributed.')
        shape_local = split_shape(shape, comm)
        return Tod(np.ones(shape_local, dtype, order), mask, nsamples, header,
                   unit, derived_units, dtype, copy=False, shape_global=shape,
                   comm=comm)

    @staticmethod
    def zeros(shape, mask=None, nsamples=None, header=None, unit=None,
              derived_units=None, dtype=None, order=None, comm=MPI.COMM_SELF):
        nsamples = nsamples or shape[-1]
        shape = flatten_sliced_shape(shape)
        if len(shape) <= 1 and comm.Get_size() > 1:
            raise ValueError('Scalar or vector Tod can not be distributed.')
        shape_local = split_shape(shape, comm)
        return Tod(np.zeros(shape_local, dtype, order), mask, nsamples, header,
                   unit, derived_units, dtype, copy=False, shape_global=shape,
                   comm=comm)

    def flatten(self, order='C'):
        """
        Return a copy of the array collapsed into one dimension.
        """
        result = super(self.__class__, self).flatten(order)
        result.nsamples = None
        return result

    def ravel(self, order='C'):
        """
        Return a flattened view of the array
        """
        result = super(self.__class__, self).ravel(order)
        result.nsamples = None
        return result
        
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
        nslices = len(self.nsamples)
        if nslices > 1:
            strsamples = ','.join((str(self.nsamples[i])
                                   for i in range(nslices)))
            output += ' in ' + str(nslices) + ' slices ('+strsamples+')'
        return output + ']'

    def save(self, filename, fitskw=None):
        fitskw = fitskw or {}
        if 'NSAMPLES' not in [k.upper() for k in fitskw.keys()]:
            fitskw['NSAMPLES'] = (self.nsamples, 'Number of samples per slice')
        FitsArray.save(self, filename, fitskw=fitskw)
        if self.mask is not None:
            write_fits(filename, self.mask.view('uint8'), None,
                       self.shape_global, True, self.comm, extname='Mask')
        _save_derived_units(filename, self.derived_units, self.comm)

    def toglobal(self):
        return DistributedArray.toglobal(self, attr=('mask',))

    def tolocal(self, comm=MPI.COMM_WORLD):
        output = DistributedArray.tolocal(self, comm)
        if self.mask is not None:
            tolocal = DistributedArray.__dict__['tolocal']
            output.mask = tolocal(self.mask, comm)
        return output


#-------------------------------------------------------------------------------


def flatten_sliced_shape(shape):
    if shape is None:
        return shape
    if not shape:
        return ()
    if not isinstance(shape, (list, tuple, np.ndarray)):
        return (int(shape),)
    shape_1 = shape[-1]
    if not isinstance(shape_1, (list, tuple, np.ndarray)):
        return tuple(shape)
    shape = list(shape)
    shape[-1] = sum(tuple(shape_1))
    return tuple(shape)

   
#-------------------------------------------------------------------------------


def combine_sliced_shape(shape, nsamples):
    if not isinstance(shape, (list, tuple, np.ndarray)):
        shape = [ shape ]
    else:
        shape = list(shape) # list makes a shallow copy
    if not isinstance(nsamples, (list, tuple, np.ndarray)):
        nsamples = (int(nsamples),)
    else:
        nsamples = tuple(nsamples)
    shape.append(nsamples)
    return tuple(shape)

   
#-------------------------------------------------------------------------------


def validate_sliced_shape(shape, nsamples=None):

    if isinstance(shape, np.ndarray):
        shape = tuple(shape)

    if not shape:
        if nsamples:
            raise ValueError("The input is scalar or None, but nsamples is eq" \
                             "ual to '" + str(nsamples) + "'.")
        return None if shape is None else ()

    if not isinstance(shape, (list, tuple, np.ndarray)):
        shape = (int(shape),)

    shape_1 = shape[-1]
    if nsamples is None:
        if isinstance(shape_1, (list, np.ndarray)):
            shape = list(shape)
            shape[-1] = tuple(shape_1)
        return tuple(shape)
    
    if not isinstance(nsamples, (list, tuple, np.ndarray)):
        nsamples = (nsamples,)
    else:
        nsamples = tuple(nsamples)

    if len(nsamples) == 0:
        raise ValueError('The input is not scalar and is incompatible with ns' \
            'amples.')

    if isinstance(shape_1, (list, tuple, np.ndarray)):
        shape_1 = tuple(shape_1)
        if shape_1 != nsamples:
            raise ValueError('The input has an incompatible slicing of samples'\
                " '" + str(nsamples) + "' instead of '" + str(shape_1) + "'.")

    sum_shape_1 = sum(shape_1) if isinstance(shape_1, tuple) else shape_1
    if sum_shape_1 != sum(nsamples):
        raise ValueError('The sliced input has an incompatible number of samp' \
            "les '" + str(nsamples) + "' instead of '" + str(shape_1) + "'.")
    
    if isinstance(shape_1, tuple):
        return tuple(shape)

    shape = list(shape)
    shape[-1] = nsamples

    return tuple(shape)


#-------------------------------------------------------------------------------


def _save_derived_units(filename, du, comm):
    if not du:
        return
    if comm is None:
        return
    if comm.Get_rank() == 0:
        buffer = StringIO.StringIO()
        pickle.dump(du, buffer, pickle.HIGHEST_PROTOCOL)
        data = np.frombuffer(buffer.getvalue(), np.uint8)
        header = create_fitsheader(fromdata=data, extname='derived_units')
        pyfits.append(filename, data, header)
    comm.Barrier()
