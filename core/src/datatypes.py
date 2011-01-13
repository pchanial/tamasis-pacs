import kapteyn.maputils
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pyfits
try:
    import ds9
    _imported_ds9 = True
except:
    _imported_ds9 = False

import tamasisfortran as tmf

from functools import reduce
from .numpyutils import _my_isscalar
from .quantity import Quantity, UnitError, _extract_unit, _strunit

__all__ = [ 'FitsArray', 'Map', 'Tod', 'create_fitsheader' ]


class FitsArray(Quantity):

    __slots__ = ('_header',)
    def __new__(cls, data, header=None, unit=None, derived_units=None, dtype=None, copy=True, order='C', subok=False, ndmin=0):

        if type(data) is str:
            ihdu = 0
            while True:
                try:
                    hdu = pyfits.open(data)[ihdu]
                except IndexError:
                    raise IOError('The FITS file has no data.')
                if hdu.data is not None:
                    break
                ihdu += 1
            data = hdu.data
            header = hdu.header
            copy = False
            if unit is None:
                if 'bunit' in header:
                    unit = header['bunit']
                elif 'QTTY____' in header:
                    unit = header['QTTY____'] # HCSS crap

        # get a new FitsArray instance (or a subclass if subok is True)
        result = Quantity.__new__(cls, data, unit, derived_units, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)

        # copy header attribute
        if header is not None:
            result.header = header
        elif hasattr(data, '_header') and data._header.__class__ is pyfits.Header:
            if copy:
                result._header = data._header.copy()
            else:
                result._header = data._header
        elif not numpy.iscomplexobj(result):
            result._header = create_fitsheader(result)
        else:
            result._header = None

        return result

    def __array_finalize__(self, obj):
        Quantity.__array_finalize__(self, obj)
        # obj might be None without __new__ being called... (ex: append)
        #if obj is None: return
        self._header = getattr(obj, '_header', None)

    def __array_wrap__(self, obj, context=None):
        result = Quantity.__array_wrap__(self, obj, context).view(type(self))
        result._header = self._header
        return result

    @staticmethod
    def empty(shape, header=None, unit=None, derived_units=None, dtype=None, order=None):
        return FitsArray(numpy.empty(shape, dtype, order), header, unit, derived_units, dtype, copy=False)

    @staticmethod
    def ones(shape, header=None, unit=None, derived_units=None, dtype=None, order=None):
        return FitsArray(numpy.ones(shape, dtype, order), header, unit, derived_units, dtype, copy=False)

    @staticmethod
    def zeros(shape, header=None, unit=None, derived_units=None, dtype=None, order=None):
        return FitsArray(numpy.zeros(shape, dtype, order), header, unit, derived_units, dtype, copy=False)

    def copy(self, order='C'):
        return FitsArray(self, copy=True, order=order)

    def has_wcs(self):
        """
        Returns True is the array has a FITS header with a defined World
        Coordinate System.
        """
        if self.header is None:
            return False

        required = 'CRPIX,CRVAL,CTYPE'.split(',')
        keywords = numpy.concatenate(
            [(lambda i: [r+str(i+1) for r in required])(i) 
             for i in range(self.header['NAXIS'])])

        return all([k in self.header for k in keywords])

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header ('+str(type(header))+').')
        self._header = header

    def tofile(self, fid, sep='', format='%s'):
        super(FitsArray,self).tofile(fid, sep, format)

    def save(self, filename):
        """Save a FitsArray instance to a fits file given a filename
       
        If the same file already exist it overwrites it.
        """

        if self.header is not None:
            header = self.header.copy()
        else:
            header = create_fitsheader(self)
       
        unit = self.unit
        if unit != '':
            header.update('BUNIT', unit)

        if hasattr(self, 'nsamples'):
            header.update('NSAMPLES', str(self.nsamples))

        if numpy.rank(self) == 0:
            value = self.reshape((1,))
        else:
            value = self.T if numpy.isfortran(self) else self
        hdu = pyfits.PrimaryHDU(value, header)
        hdu.writeto(filename, clobber=True)

    def imsave(self, filename, colorbar=True, **kw):
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        is_interactive = matplotlib.is_interactive()
        matplotlib.interactive(False)
        dpi = 80.
        figsize = numpy.clip(numpy.max(numpy.array(self.shape[::-1])/dpi), 8, 50)
        figsize = (figsize + (2 if colorbar else 0), figsize)
        fig = pyplot.figure(figsize=figsize, dpi=dpi)
        self.imshow(colorbar=colorbar, new_figure=False, **kw)
        fig = pyplot.gcf()
        fig.savefig(filename)
        matplotlib.interactive(is_interactive)

    def imshow(self, mask=None, title=None, colorbar=True, aspect=None, interpolation='nearest', origin=None, xlabel='', ylabel='', new_figure=True, **kw):
        """
        A simple graphical display function for the Tod class

        mask: array-like
            True means masked.
        """

        unfinite = ~numpy.isfinite(numpy.asarray(self))
        if mask is None:
            mask = unfinite
        else:
            mask = numpy.logical_or(mask, unfinite)

        data = numpy.ma.MaskedArray(numpy.asarray(self), mask=mask, copy=False)
        if numpy.iscomplexobj(data):
            data = abs(data)
        mean   = numpy.mean(data)
        stddev = numpy.std(data)
        # casting to float because of a bug numpy1.4 + matplotlib
        minval = float(max(mean - 2*stddev, numpy.min(data)))
        maxval = float(min(mean + 5*stddev, numpy.max(data)))

        if new_figure:
            fig = pyplot.figure()
        else:
            fig = pyplot.gcf()
        fontsize = 12. * fig.get_figheight() / 6.125

        image = pyplot.imshow(data, aspect=aspect, interpolation=interpolation, origin=origin, **kw)
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
            -zoom to fit
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
        >>> d=m.ds9('saveimage png myfits.png', scale='histequ', 
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

            if 'zoom' not in keywords:
                keywords['zoom'] = 'to fit'

            if 'scale' not in keywords:
                keywords['scale'] = ('scope local', 'mode 99.5')

            if origin == 'upper' or 'orient' not in keywords and self.origin == 'upper':
                keywords['orient'] = 'y'

            wait = 10

            id = 'ds9_' + str(uuid.uuid1())[4:8]

            command = 'ds9 -title ' + id

            for k,v in keywords.items():
                k = str(k)
                if type(v) is not tuple:
                    v = (v,)
                command += reduce(lambda x,y: str(x) + ' -' + k + ' ' + str(y),v,'')

            os.system(command + '&')

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
        d.set_np2arr(self.view(numpy.ndarray).T)

        if self.has_wcs():
            d.set('wcs append', str(self.header))
    
        if xpamsg is not None:
            if isinstance(xpamsg, str):
                xpamsg = (xpamsg,)
            for v in xpamsg:
                d.set(v)

        return d

    @property
    def unit(self):
        """
        Return the Quantity unit as a string
        
        Example
        -------
        >>> Quantity(32., 'm/s').unit
        'm / s'
        """
        return super(FitsArray, self).unit

    @unit.setter
    def unit(self, unit):
        """
        Convert a Quantity into a new unit

        Parameters
        ----------
        unit : string or None
             A string representing the unit into which the Quantity
             should be converted

        Example
        -------
        >>> q = Quantity(1., 'km')
        >>> q.unit = 'm'
        >>> print(q)
        1000.0 m
        """
        newunit = _extract_unit(unit)
        if self._unit is None or newunit is None:
            self._unit = newunit
            return

        q1 = FitsArray(1., self.header, self._unit, self.derived_units).SI
        q2 = FitsArray(1., self.header, newunit, self.derived_units).SI
        if q1._unit != q2._unit:
            raise UnitError("Units '" + self.unit + "' and '" + \
                            _strunit(newunit) + "' are incompatible.")
        factor = q1.magnitude / q2.magnitude
        if numpy.rank(self) == 0:
            self.magnitude = self.magnitude * factor
        else:
            self.magnitude.T[:] *= factor
        self._unit = newunit


#-------------------------------------------------------------------------------


class Map(FitsArray):

    __slots__ = ('coverage', 'error', 'origin')

    """
    Represent a map, complemented with unit and FITS header.
    """
    def __new__(cls, data,  header=None, unit=None, derived_units=None, coverage=None, error=None, origin='lower', dtype=None, copy=True, order='C', subok=False, ndmin=0):

        # get a new Map instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)

        if type(data) is str:
            try:
                if error is None: error = pyfits.open(data)['Error'].data
            except:
                pass
            try:
                if coverage is None: coverage = pyfits.open(data)['Coverage'].data
            except:
                pass

        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        if origin is not None:
            origin = origin.strip().lower()
            if origin not in ('upper', 'lower'):
                raise ValueError("Invalid origin '"+origin+"'. Valid values are None, 'upper' or 'lower'.")
        
        result.coverage = coverage
        result.error = error
        result.origin = origin

        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None:
            self.coverage = None
            self.error = None
            self.origin = 'lower'
            return

        self.coverage = getattr(obj, 'coverage', None)
        self.error = getattr(obj, 'error', None)
        self.origin = getattr(obj, 'origin', 'lower')

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
            raise TypeError('Incorrect type for the input header ('+str(type(header))+').')
        self._header = header
        if self.has_wcs():
            area = 0.
            if all([key in header for key in ('CD1_1', 'CD2_1', 'CD1_2', 'CD2_2')]):
                cd = numpy.array([ [header['cd1_1'],header['cd1_2']], [header['cd2_1'],header['cd2_2']] ])
                area = abs(numpy.linalg.det(cd))
            if all([key in header for key in ('CDELT1', 'CDELT2')]):
                area = abs(header['CDELT1'] * header['CDELT2'])
            if area != 0:
                self.derived_units.update({'pixel': (projection_scale, Quantity(1., 'pixel_reference')), 'pixel_reference' : Quantity(area, 'deg^2').tounit('arcsec^2')})

    @staticmethod
    def empty(shape, coverage=None, error=None, origin='lower', header=None, unit=None, derived_units=None, dtype=None, order=None):
        return Map(numpy.empty(shape, dtype, order), header, unit, derived_units, coverage, error, origin, dtype, copy=False)

    @staticmethod
    def ones(shape, coverage=None, error=None, origin='lower', header=None, unit=None, derived_units=None, dtype=None, order=None):
        return Map(numpy.ones(shape, dtype, order), header, unit, derived_units, coverage, error, origin, dtype, copy=False)

    @staticmethod
    def zeros(shape, coverage=None, error=None, origin='lower', header=None, unit=None, derived_units=None, dtype=None, order=None):
        return Map(numpy.zeros(shape, dtype, order), header, unit, derived_units, coverage, error, origin, dtype, copy=False)

    def copy(self, order='C'):
        return Map(self, copy=True, order=order)

    def imshow(self, mask=None, title=None, new_figure=True, origin=None, **kw):
        """A simple graphical display function for the Map class"""

        if mask is None and self.coverage is not None:
            mask = self.coverage <= 0
        if mask is not None:
            data = numpy.array(self, copy=True)
            data[mask] = numpy.nan
        else:
            data = numpy.asarray(self)

        if origin is None:
            origin = self.origin

        # check if the map has no astrometry information
        if not self.has_wcs():
            if 'xlabel' not in kw:
                kw['xlabel'] = 'X'
            if 'ylabel' not in kw:
                kw['ylabel'] = 'Y'
            image = super(Map, self).imshow(title=title, new_figure=new_figure, origin=origin, **kw)
            return image

        fitsobj = kapteyn.maputils.FITSimage(externaldata=data, externalheader=self.header)
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
        super(Map, self).save(filename)
        if self.error is not None:
            header = create_fitsheader(self.error, extname='Error')
            pyfits.append(filename, self.error, header)
        if self.coverage is not None:
            header = create_fitsheader(self.coverage, extname='Coverage')
            pyfits.append(filename, self.coverage, header)


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    __slots__ = ('_mask', 'nsamples')

    def __new__(cls, data, mask=None, nsamples=None, header=None, unit=None, derived_units=None, dtype=None, copy=True, order='C', subok=False, ndmin=0):

        # get a new Tod instance (or a subclass if subok is True)
        result = FitsArray.__new__(cls, data, header, unit, derived_units, dtype, copy, order, True, ndmin)
        if not subok and result.__class__ is not cls:
            result = result.view(cls)
        
        # mask attribute
        if mask is numpy.ma.nomask:
            mask = None

        if mask is None and isinstance(data, str):
            try:
                mask = pyfits.open(data)['Mask'].data.view(numpy.bool8)
                copy = False
            except:
                pass

        if mask is None and hasattr(data, 'mask') and data.mask is not numpy.ma.nomask:
            mask = data.mask

        if mask is not None:
            result._mask = numpy.array(mask, numpy.bool8, copy=copy)
        else:
            result._mask = None
        
        # nsamples attribute
        if type(data) is str and nsamples is None:
            if 'nsamples' in result.header:
                nsamples = result.header['nsamples'][1:-1].replace(' ', '')
                if len(nsamples) > 0:
                    nsamples = [int(float(x)) for x in nsamples.split(',') if x.strip() != '']
        if nsamples is None:
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
        if mask is None or mask is numpy.ma.nomask:
            self._mask = None
            return

        # enforce bool8 dtype
        if not isinstance(mask, numpy.ndarray):
            mask = numpy.array(mask, numpy.bool8)
        elif mask.dtype.type != numpy.bool8:
            if mask.dtype.itemsize == 1:
                mask = mask.view(numpy.bool8)
            else:
                mask = numpy.asarray(mask, numpy.bool8)

        # handle the scalar case
        if numpy.rank(mask) == 0:
            if self._mask is None:
                func = numpy.zeros if mask == 0 else numpy.ones
                self._mask = func(self.shape, dtype=numpy.bool8)
            else:
                self._mask[:] = mask
            return

        # check shape compatibility
        if self.shape != mask.shape:
            raise ValueError("The input mask has a shape '" + str(mask.shape) + \
                             "' incompatible with that of the Tod '" + str(self.shape) + "'.")
        
        self._mask = mask

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        # obj might be None without __new__ being called (ex: append, concatenate)
        if obj is None:
            self.nsamples = () if numpy.rank(self) == 0 else (self.shape[-1],)
            self.mask = None
            return
        if hasattr(obj, 'mask') and obj.mask is not numpy.ma.nomask:
            self._mask = None if obj.mask is numpy.ma.nomask else obj.mask
        self.nsamples = getattr(obj, 'nsamples', () if numpy.rank(self) == 0 else (self.shape[-1],))

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.mask = self.mask
        result.nsamples = self.nsamples
        return result

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if not isinstance(item, Tod):
            return item
        if item.mask is not None:
            item.mask = self.mask[key]
        if not isinstance(key, tuple):
            return item
        if len(key) > 1:
            if not isinstance(key[-1], slice):
                return item
            else:
                if key[-1].start is not None or key[-1].stop is not None or key[-1].step is not None:
                    item.nsamples = (item.shape[-1],)
        return item

    def reshape(self, newdims, order='C'):
        result = super(Tod, self).reshape(newdims, order=order)
        if self.mask is not None:
            result.mask = self.mask.reshape(newdims, order=order)
        return result

    @staticmethod
    def empty(shape, mask=None, nsamples=None, header=None, unit=None, derived_units=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(numpy.empty(shape_flat, dtype, order), mask, shape[-1], header, unit, derived_units, dtype, copy=False)

    @staticmethod
    def ones(shape, mask=None, nsamples=None, header=None, unit=None, derived_units=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(numpy.ones(shape_flat, dtype, order), mask, shape[-1], header, unit, derived_units, dtype, copy=False)

    @staticmethod
    def zeros(shape, mask=None, nsamples=None, header=None, unit=None, derived_units=None, dtype=None, order=None):
        shape = validate_sliced_shape(shape, nsamples)
        shape_flat = flatten_sliced_shape(shape)
        return Tod(numpy.zeros(shape_flat, dtype, order), mask, shape[-1], header, unit, derived_units, dtype, copy=False)
   
    def copy(self, order='C'):
        return Tod(self, copy=True, order=order)

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
        image = super(Tod, self).imshow(mask=self.mask, title=title, origin='upper', xlabel=xlabel, ylabel=ylabel, aspect=aspect, **kw)
        return image

    def __str__(self):
        if numpy.rank(self) == 0:
            if numpy.iscomplexobj(self):
                return 'Tod ' + str(complex(self))
            else:
                return 'Tod ' + str(float(self))
        output = 'Tod ['
        if numpy.rank(self) > 1:
            output += str(self.shape[-2])+' detector'
            if self.shape[-2] > 1:
                output += 's'
            output += ', '
        output += str(self.shape[-1]) + ' sample'
        if self.shape[-1] > 1:
            output += 's'
        nslices = len(self.nsamples)
        if nslices > 1:
            strsamples = ','.join((str(self.nsamples[i]) for i in range(nslices)))
            output += ' in ' + str(nslices) + ' slices ('+strsamples+')'
        return output + ']'

    def save(self, filename):
        super(Tod, self).save(filename)
        if self.mask is None:
            return
        header = create_fitsheader(self.mask, extname='Mask')
        pyfits.append(filename, self.mask.view(numpy.uint8), header)


#-------------------------------------------------------------------------------


def create_fitsheader(array, extname=None, crval=(0.,0.), crpix=None, ctype=('RA---TAN','DEC--TAN'), cunit='deg', cd=None, cdelt=None, naxis=None):
    """
    Return a FITS header

    Parameters
    ----------
    array : array_like
        An array from which the dimensions will be extracted. Note that
        by FITS convention, the dimension along X is the second value 
        of the array shape and that the dimension along the Y axis is 
        the first one. If None is specified, naxis keyword must be set
    extname : None or string
        if a string is specified ('' can be used), the returned header
        type will be an Image HDU (otherwise a Primary HDU)
    crval : 2 element array, optional
        Reference pixel values (FITS convention)
    crpix : 2 element array, optional
        Reference pixel (FITS convention)
    ctype : 2 element string array, optional
        Projection types
    cunit : string or 2 element string array
        Units of the CD matrix (default is degrees/pixel)
    cd : 2 x 2 array
        Astrometry parameters
            CD1_1 CD1_2
            CD2_1 CD2_2
    cdelt : 2 element array
        Physical increment at the reference pixel
    naxis : 2 element array
        (NAXIS1,NAXIS2) tuple, to be specified only if array argument is None

    Examples
    --------
    >>> map = Map.ones((10,100), unit='Jy/pixel')
    >>> map.header = create_fitsheader(map, cd=[[-1,0],[0,1]])
    """

    if array is None:
        if naxis is None:
            raise ValueError('An array argument or naxis keyword should be specified.')
        typename = 'float64'
    else:
        if not isinstance(array, numpy.ndarray):
            raise TypeError('The input is not an ndarray.')
        naxis = tuple(reversed(array.shape))
        if array.dtype.itemsize == 1:
            typename = 'uint8'
        elif array.dtype.names is not None:
            typename = None
        else:
            typename = array.dtype.name

    if type(naxis) not in (list, tuple):
        naxis = (naxis,)
    numaxis = len(naxis)

    if extname is None:
        card = pyfits.createCard('simple', True)
    else:
        card = pyfits.createCard('xtension', 'IMAGE', 'Image extension')
    header = pyfits.Header([card])
    if typename is not None:
        header.update('bitpix', pyfits.PrimaryHDU.ImgCode[typename], 'array data type')
    header.update('naxis', numaxis, 'number of array dimensions')
    for dim in range(numaxis):
        header.update('naxis'+str(dim+1), naxis[dim])
    if extname is None:
        header.update('extend', True)
    else:
        header.update('pcount', 0, 'number of parameters')
        header.update('gcount', 1, 'number of groups')
        header.update('extname', extname)

    if cd is not None:
        cd = numpy.asarray(cd, dtype=numpy.float64)
        if cd.shape != (2,2):
            raise ValueError('The CD matrix is not a 2x2 matrix.')
    else:
        if cdelt is None:
            return header
        if _my_isscalar(cdelt):
            cdelt = (-cdelt, cdelt)
        cd = numpy.array(((cdelt[0], 0), (0,cdelt[1])))

    crval = numpy.asarray(crval, dtype=numpy.float64)
    if crval.size != 2:
        raise ValueError('CRVAL does not have two elements.')

    if crpix is None:
        crpix = (numpy.array(naxis) + 1) / 2.
    else:
        crpix = numpy.asarray(crpix, dtype=numpy.float64)
    if crpix.size != 2:
        raise ValueError('CRPIX does not have two elements.')

    ctype = numpy.asarray(ctype, dtype=numpy.string_)
    if ctype.size != 2:
        raise ValueError('CTYPE does not have two elements.')

    if _my_isscalar(cunit):
        cunit = (cunit, cunit)
    cunit = numpy.asarray(cunit, dtype=numpy.string_)
    if cunit.size != 2:
        raise ValueError('CUNIT does not have two elements.')

    header.update('crval1', crval[0])
    header.update('crval2', crval[1])
    header.update('crpix1', crpix[0])
    header.update('crpix2', crpix[1])
    header.update('cd1_1' , cd[0,0])
    header.update('cd2_1' , cd[1,0])
    header.update('cd1_2' , cd[0,1])
    header.update('cd2_2' , cd[1,1])
    header.update('ctype1', ctype[0])
    header.update('ctype2', ctype[1])
    header.update('cunit1', cunit[0])
    header.update('cunit2', cunit[1])

    return header


#-------------------------------------------------------------------------------


def projection_scale(input):
    """
    Returns the pixel area in units of reference pixel.
    """
    if not isinstance(input, FitsArray):
        raise TypeError('The input is a ' + type(input).__name__ + ' instead of a FitsArray.')
    if not input.has_wcs():
        return 1.
    scale, status = tmf.projection_scale(str(input.header).replace('\n',''), input.header['NAXIS1'], input.header['NAXIS2'])
    if status != 0:
        raise RuntimeError()

    return scale
    

#-------------------------------------------------------------------------------


def flatten_sliced_shape(shape):
    if shape is None: return shape
    if _my_isscalar(shape):
        return (int(shape),)
    return tuple(map(numpy.sum, shape))

   
#-------------------------------------------------------------------------------


def combine_sliced_shape(shape, nsamples):
    if _my_isscalar(shape):
        shape = [ shape ]
    else:
        shape = list(shape) # list makes a shallow copy
    if _my_isscalar(nsamples):
        nsamples = (int(nsamples),)
    else:
        nsamples = tuple(nsamples)
    shape.append(nsamples)
    return tuple(shape)

   
#-------------------------------------------------------------------------------

    
def validate_sliced_shape(shape, nsamples=None):
    # convert shape and nsamples to tuple
    if shape is None:
        if nsamples is None:
            return None
        shape = ()
    elif _my_isscalar(shape):
        shape = (int(shape),)
    else:
        shape = tuple(shape)
    if nsamples is not None:
        if _my_isscalar(nsamples):
            nsamples = (int(nsamples),)
        else:
            nsamples = tuple(nsamples)
    
    if len(shape) == 0:
        if nsamples is not None and len(nsamples) != 0:
            raise ValueError("The input is scalar, but nsamples is equal to '"+str(nsamples)+"'.")
        return (shape,)
    
    if nsamples is None:
        if _my_isscalar(shape[-1]):
            nsamples = int(shape[-1])
        else:
            nsamples = tuple(map(int, shape[-1]))
    else:
        if len(nsamples) == 0:
            raise ValueError("The input is not scalar and is incompatible with nsamples.")
        if numpy.any(numpy.array(nsamples) < 0):
            raise ValueError('The input nsamples has negative values.')
        elif numpy.sum(nsamples) != numpy.sum(shape[-1]):
            raise ValueError("The sliced input has an incompatible number of samples '" + str(nsamples) + "' instead of '" + str(shape[-1]) + "'.")
        if len(nsamples) == 1:
            nsamples = nsamples[0]
    
    l = list(shape[0:-1])
    l.append(nsamples)
    return tuple(l)
