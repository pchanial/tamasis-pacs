#-------------------------------------------------------------------------------
#
# Map class
#
#-------------------------------------------------------------------------------



class Map(numpy.ndarray):
    """
    A class containing  data and metadata representing a  sky map or a
    projection  It is  an numpy  array supplemented  with localization
    metadata and pixel physical size.
    Author: N. Barbey
    """
    def __new__(subtype, shape, dtype='float64', buffer=None, offset=0,
                strides=None, order='F',
                position=numpy.array((0, 0)), rotation=0,
                step=numpy.array((1, 1)), chop=0, crpix=numpy.array((1, 1))):
        # Create the ndarray instance of our type, given the usual
        # input arguments.  This will call the standard ndarray
        # constructor, but return an object of our type
        obj = numpy.ndarray.__new__(
            subtype, shape, dtype, buffer, offset, strides, order)
        # add the new attribute to the created instance
        obj.position = numpy.array(position)
        obj.rotation = rotation
        obj.step = numpy.array(step)
        obj.chop = chop
        obj.crpix = numpy.array(crpix)
        return obj

    def __array_finalize__(self, obj):
        # reset the attribute from passed original object
        self.position = getattr(obj, 'position', numpy.array((0, 0)))
        self.rotation = getattr(obj, 'rotation', 0)
        self.step = getattr(obj, 'step', numpy.array((1, 1)))
        self.chop = getattr(obj, 'chop', 0)
        self.crpix = getattr(obj, 'crpix', numpy.array((1, 1)))
        # We do not need to return anything

    def __repr__(self):
        str_rep = "values:\n %s" % self
        str_rep += "\n position:\n %s" % self.position
        str_rep += "\n rotation:\n %s" % self.rotation
        str_rep += "\n step:\n %s" % self.step
        str_rep += "\n chop:\n %s" % self.chop
        str_rep += "\n crpix:\n %s" % self.crpix
        return str_rep

    def __copy__(self, order='C'):
        array_copy = numpy.ndarray.__copy__(self[:], order)
        new_map = Map(array_copy.shape, position=self.position,
                      rotation=self.rotation, step=self.step,
                      chop=self.chop, crpix=self.crpix)
        new_map[:] = array_copy
        return new_map

    def __reduce__(self):
        object_state = list(numpy.ndarray.__reduce__(self))
        subclass_state = (self.position, self.rotation, self.step, 
                          self.chop, self.crpix)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        nd_state, own_state = state
        numpy.ndarray.__setstate__(self, nd_state)
        position, rotation, step, chop = own_state
        self.position = position
        self.rotation = rotation
        self.step = step
        self.chop = chop
        self.crpix = crpix

    def copy(self):
        return self.__copy__()

    def zeros(self):
        """Copy a Map instance replacing data values with zeros"""
        new_map = self.copy()
        new_map[:] = numpy.zeros(new_map.shape)
        return new_map

    def ones(self):
        """Copy a Map instance replacing data values with zeros"""
        new_map = self.copy()
        new_map[:] = numpy.ones(new_map.shape)
        return new_map

    def imshow(self, n_figure=None, axis=True, title_str=None):
        """A simple graphical display function for the Map class"""
        from matplotlib.pyplot import gray, figure, imshow, colorbar, \
            draw, show, xlabel, ylabel, title
        gray()
        if n_figure==None:
            figure()
        else:
            figure(n_figure)
        no_nan = self.copy()
        no_nan[numpy.where(numpy.isnan(no_nan))] = 0
        if axis:
            (x, y) = self.axis()
            imshow(no_nan, extent=(numpy.min(y), numpy.max(y),
                                   numpy.min(x), numpy.max(x)),
                   interpolation='nearest')
            xlabel("right ascension (degrees)")
            ylabel("declination")
        else:
            imshow(no_nan, interpolation='nearest')
        if title_str is not None:
            title(title_str)
        if n_figure==None:
            colorbar()
        draw()

    def axis(self):
        """Returns axis vectors of a Map"""
        return [self.position[i] + self.step[i] *
                numpy.array(range(self.shape[i])) for i in xrange(2)]

    @property
    def header(self):
        return self.header

    @header.setter
    def header(self, header):
        import pyfits
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header.') 
        self.header = header

#    def header(self):
#        """Convert a Map instance parameters to a fits-like header"""
#        header = {
#            'NAXIS' : len(self.shape),
#            'NAXIS1': self.shape[0],
#            'NAXIS2': self.shape[1],
#            'CTYPE1': 'RA---TAN',
#            'CUNIT1': 'deg',
#            'CRVAL1': self.position[0],
#            'CRPIX1': self.crpix[0],
#            'CDELT1': self.step[0],
#            'CTYPE2': 'DEC--TAN',
#            'CUNIT2': 'deg',
#            'CRVAL2': self.position[1],
#            'CRPIX2': self.crpix[1],
#            'CDELT2': self.step[1],
#            'CROTA2': self.rotation
#            }
#        return header


#    def pyfits_header(self):
#        h_dict = self.header()
#        cards = list()
#        for key in h_dict.keys():
#            cards.append(pyfits.Card(key=key, value=h_dict[key]))
#        cardlist = pyfits.CardList(cards)
#        header = pyfits.Header(cardlist)
#        return header

#    def header_str(self):
#        hs = self.pyfits_header().__str__()
#        hs = hs.replace('\n', '')
#        return hs

#    def str2fitsheader(string):
#        """
#        Convert a string into a pyfits.Header object
#        All cards are extracted from the input string until the END keyword is reached
#        """
#        import pyfits
#        header = pyfits.Header()
#        cards = header.ascardlist()
#        iline = 0
#        while (iline*80 < len(string)):
#            line = string[iline*80:(iline+1)*80]
#            cards.append(pyfits.Card().fromstring(line))
#            if line[0:3] == 'END': break
#            iline += 1
#        return header

    def writefits(self, filename):
        """Save map instance to a fits file given a filename
        
        If the same file already exist it overwrites it.
        It should correspond to the input required by the PACS simulator
        """
        from pyfits import PrimaryHDU
        from os.path import isfile
        from os import remove
        hdu = PrimaryHDU(self, self.header)
        try:
            remove(filename)
        except OSError:
            pass
        try:
            hdu.writeto(filename)
        except IOError:
            pass

#    def to_fits(self, filename):
#        """Save map instance to a fits file given a filename
#        
#        If the same file already exist it overwrites it.
#        It should correspond to the input required by the PACS simulator
#        """
#        from pyfits import PrimaryHDU
#        from os.path import isfile
#        from os import remove
#        hdu = PrimaryHDU(self)
#        h = self.header()
#        for key in h.keys():
#            hdu.header.update(key, h[key])
#        try:
#            remove(filename)
#        except OSError:
#            pass
#        try:
#            hdu.writeto(filename)
#        except IOError:
#            pass

    def bin(self, factor):
        """Bin the Map image and change step accordingly"""
        from scipy.signal import convolve2d
        from csh import map_data
        out = self[:]
        mask = numpy.ones((factor, factor))
        out = convolve2d(out, mask, 'same')[1::factor, 1::factor]
        out = map_data(out,
                       position=self.position,
                       rotation=self.rotation,
                       step=self.step * factor,
                       chop=self.chop, crpix=self.crpix)
        return out

    def crop(self, rect):
        from csh import map_data
        """Crop the Map image according to a rectangle in pixel
        coordinates"""
        out = self[:]
        out = out[rect[0]:rect[2], rect[1]:rect[3]]
        new_position = self.position + rect[0:2] * self.step
        out = map_data(out,
                       position=new_position,
                       rotation=self.rotation,
                       step=self.step,
                       chop=self.chop,
                       crpix=self.crpix)
        return out

    


