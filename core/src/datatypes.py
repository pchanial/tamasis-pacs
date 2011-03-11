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
from .numpyutils import _my_isscalar
from .wcsutils import create_fitsheader
from .quantity import Quantity, UnitError, _extract_unit, _strunit

__all__ = [ 'FitsArray', 'Map', 'Tod' ]


class FitsArray(Quantity):

    __slots__ = ('_header',)
    def __new__(cls, data, header=None, unit=None, derived_units=None,
                dtype=None, copy=True, order='C', subok=False, ndmin=0):

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
            data = hdu.data
            header = hdu.header
            copy = False
            if unit is None:
                if 'BUNIT' in header:
                    unit = header['BUNIT']
                elif 'QTTY____' in header:
                    unit = header['QTTY____'] # HCSS crap
            del header['BUNIT']
            try:
                derived_units = fits['derived_units'].data
                derived_units = pickle.loads(str(derived_units.data))
            except KeyError:
                pass

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
            self._header = create_fitsheader(self)
        return self._header

    @header.setter
    def header(self, header):
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header (' + \
                            str(type(header))+').')
        self._header = header

    def tofile(self, fid, sep='', format='%s'):
        super(FitsArray,self).tofile(fid, sep, format)

    def save(self, filename, fitskw={}):
        """Save a FitsArray instance to a fits file given a filename
       
        If the same file already exist it overwrites it.
        """

        if self.header is not None:
            header = self.header.copy()
        else:
            header = create_fitsheader(self)
       
        if len(self._unit) != 0:
            header.update('BUNIT', self.unit)

        for k,v in fitskw.items():
            if hasattr(self, k):
                header.update(v, str(getattr(self, k)))

        if np.rank(self) == 0:
            value = self.reshape((1,))
        else:
            value = self.T if np.isfortran(self) else self
        hdu = pyfits.PrimaryHDU(value, header)
        hdu.writeto(filename, clobber=True)
        if not isinstance(self, Map) and not isinstance(self, Tod):
            _save_derived_units(filename, self.derived_units)

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

    __slots__ = ('coverage', 'error', 'origin')

    """
    Represent a map, complemented with unit and FITS header.
    """
    def __new__(cls, data,  header=None, unit=None, derived_units=None,
                coverage=None, error=None, origin=None, dtype=None, copy=True,
                order='C', subok=False, ndmin=0):

        # get a new Map instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units,
                                   dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)

        if type(data) is str:
            if 'DISPORIG' in result.header:
                if origin is None:
                    origin = result.header['DISPORIG']
                del result.header['DISPORIG']
            try:
                if error is None: error = pyfits.open(data)['Error'].data
            except:
                pass
            try:
                if coverage is None: coverage = pyfits.open(data)['Coverage'].data
            except:
                pass

        if origin is not None:
            origin = origin.strip().lower()
            if origin not in ('upper', 'lower', 'none'):
                raise ValueError("Invalid origin '" + origin + "'. Expected v" \
                                 "alues are None, 'upper' or 'lower'.")
            if origin != 'none':
                result.origin = origin

        if error is not None:
            result.error = error
        elif copy and result.error is not None:
            result.error = result.error.copy()

        if coverage is not None:
            result.coverage = coverage
        elif copy and result.coverage is not None:
            result.coverage = result.coverage.copy()

        return result

    def __array_finalize__(self, array):
        FitsArray.__array_finalize__(self, array)
        self.coverage = getattr(array, 'coverage', None)
        self.error = getattr(array, 'error', None)
        self.origin = getattr(array, 'origin', 'lower')

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if not isinstance(item, Map):
            return item
        if item.coverage is not None:
            item.coverage = item.coverage[key]
        return item

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header (' + \
                            str(type(header))+').')
        self._header = header

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

    def save(self, filename):
        FitsArray.save(self, filename, fitskw={'origin':'DISPORIG'})
        if self.error is not None:
            header = create_fitsheader(self.error, extname='Error')
            pyfits.append(filename, self.error, header)
        if self.coverage is not None:
            header = create_fitsheader(self.coverage, extname='Coverage')
            pyfits.append(filename, self.coverage, header)
        _save_derived_units(filename, self.derived_units)


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    __slots__ = ('_mask', 'nsamples')

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None,
                derived_units=None, dtype=None, copy=True, order='C',
                subok=False, ndmin=0):

        # get a new Tod instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units,
                                   dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)
        
        # mask attribute
        if mask is np.ma.nomask:
            mask = None

        if mask is None and isinstance(data, str):
            try:
                mask = pyfits.open(data)['Mask'].data.view(np.bool8)
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
        if not nsamples:
            return result
        shape = validate_sliced_shape(result.shape, nsamples)
        result.nsamples = shape[-1]
        if type(result.nsamples) is not tuple:
            result.nsamples = (result.nsamples,)

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
        item = super(Quantity, self).__getitem__(key)
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
              derived_units=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(np.empty(shape_flat, dtype, order), mask, shape[-1], header,
                   unit, derived_units, dtype, copy=False)

    @staticmethod
    def ones(shape, mask=None, nsamples=None, header=None, unit=None,
             derived_units=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(np.ones(shape_flat, dtype, order), mask, shape[-1], header,
                   unit, derived_units, dtype, copy=False)

    @staticmethod
    def zeros(shape, mask=None, nsamples=None, header=None, unit=None,
              derived_units=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(np.zeros(shape_flat, dtype, order), mask, shape[-1], header,
                   unit, derived_units, dtype, copy=False)
   
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

    def save(self, filename):
        FitsArray.save(self, filename, fitskw={'nsamples':'NSAMPLES'})
        if self.mask is None:
            return
        header = create_fitsheader(self.mask, extname='Mask')
        pyfits.append(filename, self.mask.view(np.uint8), header)
        _save_derived_units(filename, self.derived_units)


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


def _save_derived_units(filename, du):
    if not du:
        return
    buffer = StringIO.StringIO()
    pickle.dump(du, buffer, pickle.HIGHEST_PROTOCOL)
    data = np.frombuffer(buffer.getvalue(), np.uint8)
    header = create_fitsheader(data, extname='derived_units')
    pyfits.append(filename, data, header)
