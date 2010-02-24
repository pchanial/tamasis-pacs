import tamasisfortran as tmmf
import numpy
    



#-------------------------------------------------------------------------------
#
# Acquisition models
#
#-------------------------------------------------------------------------------



class ValidationError(Exception): pass


#-------------------------------------------------------------------------------


class AcquisitionModel(object):
    """
    Class representing the Acquisition model M as in the equation 
         y = M.x,
    where x is the map and y is the instrumental response. An acquisition model can be the 
    combination of several submodels describing various parts of the instrumental chain:
         M = M3 * M2 * M1 ...
    Two methods are provided to the user:
         y = M.direct(x)s
         x = M.transpose(y)
    they return M.x (submodels are applied in right-to-left order) and transpose(M).y.
    Author: P. Chanial
    """
    def __init__(self, description):
        self.description = description

    #abstractmethod
    def direct(self, data):
        for model in self:
            data = model.direct(data)
        return data

    #abstractmethod
    def transpose(self, data):
        for model in reversed(self):
            data = model.transpose(data)
        return data

    def __validate_chain_direct(self, shapein):
        self.shapeout = shapein
        for model in self:
            self.shapeout = model._validate_shape_direct(self.shapeout)

    def __validate_chain_transpose(self, shapeout):
        self.shapein = shapeout
        for model in reversed(self):
            self.shapein = model._validate_shape_transpose(self.shapein)

    def _validate_chain(self):
        self.__validate_chain_direct(None)
        self.__validate_chain_transpose(None)

    def _validate_shape_direct(self, shapein):
        if shapein is None or self.shapein is None:
            if self.shapeout is not None:
               return self.shapeout
            else:
               return shapein
        if shapein != self.shapein:
            raise ValidationError('The input of '+self.__class__.__name__+' has incompatible shape '+str(shapein)+' instead of '+str(self.shapein)+'.')
        return self.shapeout

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None or self.shapeout is None:
            if self.shapein is not None:
                return self.shapein
            else:
                return shapeout
        if shapeout != self.shapeout:
            raise ValidationError('The input of '+self.__class__.__name__+' transpose has incompatible shape '+str(shapeout)+' instead of '+str(self.shapeout)+'.')
        return self.shapein

    def _validate_input_direct(self, cls, data):
        if not isinstance(data, cls):
            raise TypeError("The input of '+self.__class__.__name+' has an invalid type '"+cls.__name__+"'.")
        return self._validate_shape_direct(data.shape)

    def _validate_input_transpose(self, cls, data):
        if not isinstance(data, cls):
            raise TypeError("The input of '+self.__class__.__name+' transpose has an invalid type '"+cls.__name__+"'.")
        return self._validate_shape_transpose(data.shape)

    def _validate_output_direct(self, cls, shapeout, **options):
        if shapeout is None:
            raise ValueError('The shape of the output of '+self.__class__.__name__+' is not known.')
        if self._output_direct is not None and shapeout == self._output_direct.shape:
            return
        print 'Info: allocating '+str(numpy.product(shapeout)/2.**17)+' MiB for the output of '+self.__class__.__name__
        self._output_direct = cls.empty(shapeout, dtype=numpy.float64, order='fortran', **options)

    def _validate_output_transpose(self, cls, shapein, **options):
        if shapein is None:
            raise ValueError('The shape of the input of '+self.__class__.__name__+' is not known.')
        if self._output_transpose is not None and shapein == self._output_transpose.shape:
            return
        print 'Info: allocating '+str(numpy.product(shapein)/2.**17)+' MiB for the output of the transpose of '+self.__class__.__name__
        self._output_transpose = cls.empty(shapein, dtype=numpy.float64, order='fortran', **options)

    shapein           = None   # input is unconstrained
    shapeout          = None   # output is unconstrained
    _output_direct    = None   # stores the input of the transpose model. Its memory allocation is re-used as the output of the direct model
    _output_transpose = None   # stores the input of the direct model. Its memory allocation is re-used as the output of the transpose model
    blocks            = []     # components are ordered as in the direct order

    def __mul__(self, other):
        newmodel = AcquisitionModel(None)
        newmodel.blocks = []
        for block in other:
            newmodel.blocks.append(block)
        for block in self:
            newmodel.blocks.append(block)
        newmodel._validate_chain()
        return newmodel

    def __getitem__(self, index):
        if len(self) == 0:
            raise IndexError('This acquisition model has no component.')
        return self.blocks[index]

    def __setitem__(self, index, value):
        if len(self) == 0:
            raise IndexError('This acquisition model has no component.')

        try:
            oldvalue = self[index]
            oldshapein = self.shapein
            oldshapeout = self.shapeout
        except:
            raise IndexError('Only substitutions are allowed in an acquisition model.')

        self.blocks[index] = value
        try:
            self.validate()
        except ValidationError as inst:
            self.blocks[index] = oldvalue
            self.shapein = oldshapein
            self.shapeout = oldshapeout
            raise inst

    def __iter__(self):
        if len(self) == 0:
            return iter((self,))
        return self.blocks.__iter__()

    def __reversed__(self):
        if len(self) == 0:
            return iter((self,))
        return self.blocks.__reversed__()

    def __len__(self):
        return self.blocks.__len__()

    def __str__(self):
        result = self.__class__.__name__
        if self.description is not None:
            result = self.description+' ('+self.__class__.__name__+')'
        else:
            result = self.__class__.__name__
        if self.shapein is not None or self.shapeout is not None:
            result += ' [input:'
            result += 'unconstrained' if self.shapein is None else str(self.shapein).replace(' ','')
            result += ', output:'
            result += 'unconstrained' if self.shapeout is None else str(self.shapeout).replace(' ','')
            result += ']'
        if len(self) == 0:
            return result
        result += '\n  '+'\n  '.join((str(block) for block in reversed(self)))
        return result


#-------------------------------------------------------------------------------


class PacsProjectionSharpEdges(AcquisitionModel):
    """
    This class handles the integration of the sky inside the PACS detectors.
    It is assumed that the detectors have sharp edges, so that the integration is simply equal 
    to the sum of the sky pixels values weighted by their intersections with the detector surface.
    Author: P. Chanial
    """
    def __init__(self, pacs, description=None):
        AcquisitionModel.__init__(self, description)
        self.npixels_per_sample = pacs.npixels_per_sample
        sizeofpmatrix = pacs.npixels_per_sample * pacs.nfinesamples * pacs.ndetectors
        print 'Info: allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        self.pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        self.header = pacs.header
        self.shapein = (pacs.header["naxis1"], pacs.header["naxis2"])
        self.shapeout = (pacs.nfinesamples, pacs.ndetectors) 
        tmmf.pacs_pointing_matrix(pacs.array, pacs.pointing_time, pacs.ra, pacs.dec, pacs.pa, pacs.chop, pacs.fine_time, pacs.npixels_per_sample, pacs.ndetectors, pacs.transparent_mode, pacs.bad_pixel_mask.astype('int8'), pacs.keep_bad_detectors, str(pacs.header).replace('\n', ''), self.pmatrix)

    def direct(self, map2d):
        self._validate_output_direct(Tod, self._validate_input_direct(Map, map2d))
        map2d.header = self.header
        self._output_transpose = map2d
        tmmf.pacs_projection_sharp_edges_direct(self.pmatrix, map2d, self._output_direct, self.npixels_per_sample)
        return self._output_direct

    def transpose(self, signal):
        self._validate_output_transpose(Map, self._validate_input_transpose(Tod, signal), header=self.header)
        self._output_direct = signal
        tmmf.pacs_projection_sharp_edges_transpose(self.pmatrix, signal, self._output_transpose,  self.npixels_per_sample)
        return self._output_transpose

    def __str__(self):
        return super(PacsProjectionSharpEdges, self).__str__()#+' => '+self.filename+' ['+str(self.first)+','+str(self.last)+']'


#-------------------------------------------------------------------------------


class PacsMultiplexing(AcquisitionModel):
    """
    Performs the multiplexing of the PACS subarrays. The subarray columns are read one after the
    other, in a 0.025s cycle (40Hz).
    Author: P. Chanial
    """
    def __init__(self, pacs, description=None):
        AcquisitionModel.__init__(self, description)
        self.fine_sampling_factor = pacs.fine_sampling_factor
        self.ij = pacs.ij

    def direct(self, signal):
        self._validate_output_direct(Tod, self._validate_input_direct(Tod, signal))
        self._output_transpose = signal
        tmmf.pacs_multiplexing_direct(signal, self._output_direct, self.fine_sampling_factor, self.ij)
        return self._output_direct

    def transpose(self, signal):
        self._validate_output_transpose(Tod, self._validate_input_transpose(Tod, signal))
        self._output_direct = signal
        tmmf.pacs_multiplexing_transpose(signal, self._output_transpose, self.fine_sampling_factor, self.ij)
        return self._output_transpose

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        super(PacsMultiplexing, self)._validate_shape_direct(shapein)
        if shapein[0] % self.fine_sampling_factor != 0:
            raise ValidationError('The input timeline size ('+str(shapein[0])+') is not an integer times the fine sampling factor ('+str(self.fine_sampling_factor)+').')
        shapeout = list(shapein)
        shapeout[0] = shapeout[0] / self.fine_sampling_factor
        return tuple(shapeout)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return
        super(PacsMultiplexing, self)._validate_shape_transpose(shapeout)
        shapein = list(shapeout)
        shapein[0] = shapein[0] * self.fine_sampling_factor
        return tuple(shapein)


#-------------------------------------------------------------------------------


class Compression(AcquisitionModel):
    """
    Superclass for compressing the input signal.
    Author: P. Chanial
    """
    def __init__(self, compression_direct, compression_transpose, compression_factor, description):
        AcquisitionModel.__init__(self, description)
        self.compression_direct = compression_direct
        self.compression_transpose = compression_transpose
        self.compression_factor = compression_factor

    def direct(self, signal):
        if self.compression_factor == 1:
            return signal
        self._validate_output_direct(Tod, self._validate_input_direct(Tod, signal))
        self._output_transpose = signal
        self.compression_direct(signal, self._output_direct, self.compression_factor)
        return self._output_direct

    def transpose(self, compressed):
        if self.compression_factor == 1:
            return compressed
        self._validate_output_transpose(Tod, self._validate_input_transpose(Tod, compressed))
        self._output_direct = compressed
        self.compression_transpose(compressed, self._output_transpose, self.compression_factor)
        return self._output_transpose

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        super(Compression, self)._validate_shape_direct(shapein)
        if shapein[0] % self.compression_factor != 0:
            raise ValidationError('The input timeline size ('+str(shapein[0])+') is not an integer times the compression factor ('+str(self.compression_factor)+').')
        shapeout = list(shapein)
        shapeout[0] /= self.compression_factor
        return tuple(shapeout)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return
        super(Compression, self)._validate_shape_transpose(shapeout)
        shapein = list(shapeout)
        shapein[0] *= self.compression_factor
        return tuple(shapein)

    def __str__(self):
        return super(Compression, self).__str__()+' (x'+str(self.compression_factor)+')'
        

#-------------------------------------------------------------------------------


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, tmmf.compression_average_direct, tmmf.compression_average_transpose, compression_factor, description)
        

#-------------------------------------------------------------------------------


class Identity(AcquisitionModel):
    """
    Do nothing class.
    Author: P. Chanial
    """
    def __init__(self, description=None):
        AcquisitionModel.__init__(self, description)

    def direct(self, signal):
        return signal

    def transpose(self, signal):
        return signal
        

#-------------------------------------------------------------------------------


class Masking(AcquisitionModel):
    """
    Apply a mask. Ensure that the output is a FitsMaskedArray.
    Author: P. Chanial
    """
    def __init__(self, mask, description=None):
        AcquisitionModel.__init__(self, description)
        if mask is None:
            self.mask = numpy.ma.nomask
        else:
            self.mask = mask

    def direct(self, signal):
        if isinstance(signal, FitsMaskedArray):
            signal.mask = self.mask
        else:
            signal = FitsMaskedArray(signal, mask=self.mask)
        numpy.putmask(signal, signal.mask, 0.)
        return signal

    def transpose(self, signal):
        return self.direct(signal)
   



#-------------------------------------------------------------------------------
#
# PACS observations & simulations
#
#-------------------------------------------------------------------------------



class _Pacs():

    """
    Base class which encapsulates handy information about the PACS instrument and the processed 
    observation. See subclasses PacsObservation and PacsSimulation for details.
    Author: P. Chanial
    """

    sampling = 0.025 # 40 Hz sampling

    def __init__(self, array, pointing_time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_pixel_mask, keep_bad_detectors):

        from os import getenv
        from os.path import join
        from pyfits import fitsopen

        if array.lower() not in ('blue', 'green', 'red'):
            raise ValueError("The input array is not 'blue', 'green' nor 'red'.")

        if observing_mode is not None:
            observing_mode = observing_mode.lower()

        if observing_mode not in (None, 'prime', 'parallel', 'transparent'):
            raise ValueError("Observing mode is not 'prime', 'parallel' nor 'transparent'.")

        if compression_factor is None:
            if observing_mode is None:
                raise ValueError('The compression factor is not specified.')
            compression_factor = {'prime':4, 'parallel':8, 'transparent':1}[observing_mode]

        self.transparent_mode = observing_mode == 'transparent'

        if (fine_sampling_factor not in (1,2,4,8,16)):
            raise ValueError('Invalid fine sampling factor. It may be 1, 2, 4, 8 or 16')
 
        # 'blue', 'green' or 'red' array
        self.array = array.lower()

        # astrometry
        if pointing_time is None or ra is None or dec is None or pa is None or chop is None:
            raise ValueError('The simulated scan astrometry is not defined.')
        shape = pointing_time.shape
        if ra.shape != shape or dec.shape != shape or pa.shape != shape or chop.shape != shape:
            raise ValueError('The input time, ra, dec, pa, chope do not have the same shape.')
        self.pointing_time = pointing_time
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.chop = chop
        self.header = None

        # sampling information
        self.fine_time = fine_time
        self.nfinesamples = fine_time.size
        self.npixels_per_sample = npixels_per_sample
        self.fine_sampling_factor = fine_sampling_factor
        self.compression_factor = compression_factor
        
        if bad_pixel_mask is None:
            bad_pixel_maskFile = join(getenv('TAMASIS_DIR'),'data','PCalPhotometer_BadPixelMask_FM_v3.fits')
            bad_pixel_mask = numpy.array(fitsopen(bad_pixel_maskFile)[self.array].data, order='fortran')
        else:
            array_shapes = {'blue':(32,64), 'green':(32,64), 'red':(16,32)}
            if bad_pixel_mask.shape != array_shapes[self.array]:
                raise ValueError('Input bad pixel mask has incorrect shape '+str(bad_pixel_mask.shape)+' instead of '+str(array_shapes[self.array]))
        self.bad_pixel_mask = bad_pixel_mask.copy('fortran').astype(bool)
        if self.transparent_mode:
            self.bad_pixel_mask[:,0:16] = True
            self.bad_pixel_mask[16:32,16:32] = True
            self.bad_pixel_mask[:,32:] = True
        self.keep_bad_detectors = keep_bad_detectors

        self.ndetectors = self.bad_pixel_mask.size
        if not self.keep_bad_detectors:
            self.ndetectors -= int(numpy.sum(self.bad_pixel_mask))
        
        #XXX NOT IMPLEMENTED!
        self.ij = None


    @staticmethod
    def _str2fitsheader(string):
        """
        Convert a string into a pyfits.Header object
        All cards are extracted from the input string until the END keyword is reached.
        """
        import pyfits
        header = pyfits.Header()
        cards = header.ascardlist()
        iline = 0
        while (iline*80 < len(string)):
            line = string[iline*80:(iline+1)*80]
            if line[0:3] == 'END': break
            cards.append(pyfits.Card().fromstring(line))
            iline += 1
        return header
    

#-------------------------------------------------------------------------------


class PacsObservation(_Pacs):
    """
    Class which encapsulates handy information about the PACS instrument and the processed 
    observation. It contains the following attributes:
    - filename           : name of the file name, including the array colour, but excluding 
                           the type. Example: '1342184520_blue'.
    - first              : first sample to be read in the FITS files (starting from 0).
                           The associated time from the Time file is the time of the first projection.
    - last               : last sample to be read in the FITS files (excluded)
    - header             : pyfits header of the input sky map
    - npixels_per_sample   : number of sky pixels which intersect a PACS detector
    - observing_mode      : 'prime', 'parallel' or 'transparent'
    - fine_sampling_factor : number of frames to be processed during the 0.025s period.
    - compression_factor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - bad_pixel_mask       : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, filename, first=0, last=None, header=None, resolution=None, npixels_per_sample=9, fine_sampling_factor=1, bad_pixel_mask=None, keep_bad_detectors=False):

        import pyfits

        if last is None:
            ntotsamples = tmmf.pacs_info_nsamples(filename)
            print "This observation contains "+str(ntotsamples)+" samples."
            return

        if header is not None and resolution is not None:
            raise ValueError('It is not possible to specify a FITS header and an image resolution.')

        if filename[-4:].lower() == 'blue':
            array = 'blue'
        elif filename[-5:].lower() == 'green':
            array = 'green'
        elif filename[-3:].lower() == 'red':
            array = 'red'
        else:
            array = ''

        # reading astrometry
        time = numpy.array(pyfits.fitsopen(filename+'_Time.fits'        )[1].data[first:last], dtype='float64')*1.e-6
        ra   = numpy.array(pyfits.fitsopen(filename+'_RaArray.fits'     )[1].data[first:last], dtype='float64')
        dec  = numpy.array(pyfits.fitsopen(filename+'_DecArray.fits'    )[1].data[first:last], dtype='float64')
        pa   = numpy.array(pyfits.fitsopen(filename+'_PaArray.fits'     )[1].data[first:last], dtype='float64')
        chop = numpy.array(pyfits.fitsopen(filename+'_ChopFpuAngle.fits')[1].data[first:last], dtype='float64')

        #because telemetry drops can occur, we compute a fine time for each sample in the files
        compression = (time[1]-time[0]) / _Pacs.sampling
        compression_factor = int(compression + 0.5)
        if abs(compression_factor - compression) > 1.001 * compression:
            raise ValueError('The file '+filename+'_Time.fits is not sampled an integer time of '+str(_Pacs.sampling)+'s.')
        if compression_factor == 1:
            observing_mode = 'transparent'
        elif compression_factor == 4:
            observing_mode = 'prime'
        elif compression_factor == 8:
            observing_mode = 'parallel'
        else:
            observing_mode = None
        sampling_factor = fine_sampling_factor * compression_factor
        fine_time = numpy.zeros((last - first) * sampling_factor, dtype='float64')
        for i in range(sampling_factor):
            fine_time[i::sampling_factor] = time + i * _Pacs.sampling / fine_sampling_factor

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_pixel_mask, keep_bad_detectors)

        if resolution is None:
            resolution = 3.
        if header is None:
            header = self._str2fitsheader(tmmf.pacs_map_header(array, time, ra, dec, pa, chop, fine_time, self.transparent_mode, self.bad_pixel_mask.astype('int8'), keep_bad_detectors, resolution))

        self.header = header

        self.filename = filename
        self.first = first
        self.last = last
        
        #XXX REMOVE ME
        self.ij2 = tmmf.pacs_info_ij(filename, True, self.ndetectors)
 
    def get_tod(self):
        """
        Returns the signal and mask timelines as fortran arrays.
        """
        signal, mask = tmmf.pacs_timeline(self.filename, self.array, self.first+1, self.last, self.ndetectors, self.transparent_mode, self.bad_pixel_mask, self.keep_bad_detectors)
        mask = numpy.array(mask, dtype=bool)
        if self.keep_bad_detectors:
           for idetector, badpixel in enumerate(self.bad_pixel_mask.flat):
               if badpixel:
                   mask[:,idetector] = True
        return Tod(signal, mask=mask)
   

#-------------------------------------------------------------------------------


class PacsSimulation(_Pacs):
    """
    Class which encapsulates handy information describing the simulation and the setting of the PACS instrument. 
    It contains the following attributes:
    - time               : time vector. time[0] is the time of the first sky projection
    - ra                 : boresight right ascension vector
    - dec                : boresight declination vector
    - pa                 : position Angle vector
    - chop               : chop angle vector
    - header             : pyfits header of the input sky map
    - ndetectors         : number of detectors involved in the processing
    - npixels_per_sample   : number of sky pixels which intersect a PACS detector
    - observing_mode      : 'prime', 'parallel' or 'transparent'
    - fine_sampling_factor : number of frames to be processed during the 1/40Hz period.
    - compression_factor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - bad_pixel_mask       : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, inputmap, time=None, ra=None, pa=None, chop=None, array='blue', npixels_per_sample=9, observing_mode='prime', fine_sampling_factor=1, compression_factor=None, bad_pixel_mask=None, keep_bad_detectors=False):

        if not isinstance(inputmap, Map):
            raise TypeError('The input is not a Map.')

        if inputmap.header is None:
            raise TypeError('The input map header is not known.')

        fine_time = numpy.arange(time[0], floor((time[-1] - time[0]) / pacs.sampling) * fine_sampling_factor, _Pacs.sampling/fine_sampling_factor)

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_pixel_mask, keep_bad_detectors)
 
        self.inputmap
        self.header = inputmap.header

    def get_tod(self, model, noise=None):
        """
        Returns simulated signal and mask (=None) timelines.
        """
        signal = model.direct(self.inputmap)
        if self.keep_bad_detectors:
            mask = numpy.zeros(signal.shape, dtype=bool, order='fortran')
            for idetector, badpixel in enumerate(self.bad_pixel_mask.flat):
                if badpixel:
                    mask[:,idetector] = True
        else:
            mask = None

        return Tod(signal, mask=mask)
    



#-------------------------------------------------------------------------------
#
# Map & Tod class
#
#-------------------------------------------------------------------------------


class FitsMaskedArray(numpy.ma.MaskedArray):

    def __new__(cls, data, dtype=None, copy=False, mask=numpy.ma.nomask, header=None):
        result = numpy.ma.MaskedArray(data, dtype=dtype, copy=copy, mask=mask)
        result = result.view(cls)
        if header is not None:
            result.header = header
        elif hasattr(data, 'header'):
            result.header = data.header
        else:
            result.header = None
        return result

    def __array_finalize__(self, obj):
        numpy.ma.MaskedArray.__array_finalize__(self, obj)
        if obj is None: return
        self.header = getattr(obj, '_header', None)

    def __array_wrap__(self, obj, context=None):
        result = numpy.ma.MaskedArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.header = self.header
        return result

    @staticmethod
    def empty(shape, dtype='float64', order='C', header=None):
        return FitsMaskedArray(numpy.ma.empty(shape, dtype, order), header=header)

    @staticmethod
    def ones(shape, dtype='float64', order='C', header=None):
        return FitsMaskedArray(numpy.ma.ones(shape, dtype, order), header=header)

    @staticmethod
    def zeros(shape, dtype='float64', order='C', header=None):
        return FitsMaskedArray(numpy.ma.zeros(shape, dtype, order), header=header)

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        import pyfits
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header ('+str(type(header))+').') 
        self._header = header

    def tofile(self, fid, sep="", format="%s"):
        super(numpy.ma.MaskedArray,self).tofile(fid, sep, format)

    def writefits(self, filename):
        """Save a FitsMaskedArray instance to a fits file given a filename
        
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


#-------------------------------------------------------------------------------


class Map(FitsMaskedArray):

    def __new__(cls, data, dtype=None, copy=False, mask=numpy.ma.nomask, header=None):
        result = FitsMaskedArray(data, dtype=dtype, copy=copy, mask=mask, header=header)
        result = result.view(cls)
        return result

    def __array_finalize__(self, obj):
        FitsMaskedArray.__array_finalize__(self, obj)
        if obj is None: return
        self.mytoddata = getattr(obj, 'mytoddata', None)

    def __array_wrap__(self, obj, context=None):
        result = FitsMaskedArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.header = self.header
        return result

    @staticmethod
    def empty(shape, dtype='float64', order='C', header=None):
        return Map(FitsMaskedArray.empty(shape, dtype, order, header=header))

    @staticmethod
    def ones(shape, dtype='float64', order='C', header=None):
        return Map(FitsMaskedArray.ones(shape, dtype, order, header=header))

    @staticmethod
    def zeros(shape, dtype='float64', order='C', header=None):
        return Map(FitsMaskedArray.zeros(shape, dtype, order, header=header))

    def imshow(self, num=None, axis=True, title=None):
        """A simple graphical display function for the Map class"""
        from matplotlib.pyplot import gray, figure, imshow, colorbar, \
            draw, show, xlabel, ylabel, title as pyplottitle

        if self.mask.all():
            raise ValueError('All pixels are masked.')

        #gray()
        figure(num=num)
        imshow(self, 
               interpolation='nearest', 
               origin='lower')
        xlabel('Right Ascension')
        ylabel("Declination")
        if title is not None:
            pyplottitle(title)
        colorbar()
        draw()


#-------------------------------------------------------------------------------


class Tod(FitsMaskedArray):

    def __new__(cls, data, dtype=None, copy=False, mask=numpy.ma.nomask, header=None):
        result = FitsMaskedArray(data, dtype=dtype, copy=copy, mask=mask, header=header)
        result = result.view(cls)
        return result

    def __array_finalize__(self, obj):
        FitsMaskedArray.__array_finalize__(self, obj)
        if obj is None: return
        self.mytoddata = getattr(obj, 'mytoddata', None)

    def __array_wrap__(self, obj, context=None):
        result = FitsMaskedArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.header = self.header
        return result

    @staticmethod
    def empty(shape, dtype='float64', order='C', header=None):
        return Tod(FitsMaskedArray.empty(shape, dtype, order, header=header))

    @staticmethod
    def ones(shape, dtype='float64', order='C', header=None):
        return Tod(FitsMaskedArray.ones(shape, dtype, order, header=header))

    @staticmethod
    def zeros(shape, dtype='float64', order='C', header=None):
        return Tod(FitsMaskedArray.zeros(shape, dtype, order, header=header))




#-------------------------------------------------------------------------------
#
# Miscellaneous routines
#
#-------------------------------------------------------------------------------



def hcss_photproject(pacs):
    """
    Returns a map, as calculated by HCSS's PhotProject
    """
    from copy import copy
    if pacs.fine_sampling_factor != 1:
        raise ValueError('Fine sampling factor should be 1 for hcssPhotProject.') # or add decimation
    tod = pacs.get_tod()
    model = Masking(tod.mask) * PacsProjectionSharpEdges(pacs)
    mymap = copy(model.transpose(tod))
    tod[:] = 1.
    weights = model.transpose(tod)
    mymap /= weights
    return mymap

 
