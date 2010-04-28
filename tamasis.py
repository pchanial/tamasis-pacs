import tamasisfortran as tmf
import numpy
import scipy.sparse.linalg    


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
    def direct(self, data, reusein=False, copyout=True):
        for i, model in enumerate(self):
            data = model.direct(data, reusein=reusein or i != 0, copyout=copyout and i==len(self)-1)
        return data

    #abstractmethod
    def transpose(self, data, reusein=False, copyout=True):
        for i, model in enumerate(reversed(self)):
            data = model.transpose(data, reusein=reusein or i != 0, copyout=copyout and i==len(self)-1)
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

    def _validate_input_direct(self, cls, data, reusein):
        if not isinstance(data, cls):
            raise TypeError("The input of '"+self.__class__.__name__+"' has an invalid type '"+data.__class__.__name__+"' instead of '"+cls.__name__+"'.")
        shapeout = self._validate_shape_direct(data.shape)
        if reusein: self._output_transpose = data
        return data, shapeout

    def _validate_input_transpose(self, cls, data, reusein):
        if not isinstance(data, cls):
            raise TypeError("The input of the transpose of '"+self.__class__.__name__+"' has an invalid type '"+data.__class__.__name__+"' instead of '"+cls.__name__+"'.")
        shapein = self._validate_shape_transpose(data.shape)
        if reusein: self._output_direct = data
        return data, shapein

    # allocate memory for the output of the direct operation.
    # **options should not depend on the input array, since it might only be used for the first call and it would not be
    # propagated during the next calls
    def _validate_output_direct(self, cls, shapeout, **options):
        if shapeout is None:
            raise ValueError('The shape of the output of '+self.__class__.__name__+' is not known.')
        if self._output_direct is not None and shapeout == self._output_direct.shape:
            return
        print 'Info: Allocating '+str(numpy.product(shapeout)/2.**17)+' MiB for the output of '+self.__class__.__name__
        if cls == numpy.ndarray:
            self._output_direct = numpy.empty(shapeout, dtype=numpy.float64, **options)
        else:
            self._output_direct = cls.empty(shapeout, dtype=numpy.float64, **options)

    # allocate memory for the output of the transpose operation.
    # **options should not depend on the input array
    def _validate_output_transpose(self, cls, shapein, **options):
        if shapein is None:
            raise ValueError('The shape of the input of '+self.__class__.__name__+' is not known.')
        if self._output_transpose is not None and shapein == self._output_transpose.shape:
            return
        print 'Info: Allocating '+str(numpy.product(shapein)/2.**17)+' MiB for the output of the transpose of '+self.__class__.__name__
        if cls == numpy.ndarray:
            self._output_transpose = numpy.empty(shapein, dtype=numpy.float64, **options)
        else:
            self._output_transpose = cls.empty(shapein, dtype=numpy.float64, **options)

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
        result += '\n  '+'\n  '.join((str(block) for block in self))
        return result


#-------------------------------------------------------------------------------


class Projection(AcquisitionModel):
    """
    This class handles the direct and transpose operations by the pointing matrix
    The input observation has the following required attributes/methods: 
        - header
        - npixels_per_sample
        - nfinesamples
        - nsamples
        - ndetectors
        - get_pointing_matrix()
    Author: P. Chanial
    """
    def __init__(self, observation, description=None, finer_sampling=True):
        AcquisitionModel.__init__(self, description)
        self.npixels_per_sample = observation.npixels_per_sample
        self.header = observation.header
        self.pmatrix = observation.get_pointing_matrix(finer_sampling=finer_sampling)
        self.shapein = tuple([self.header['naxis'+str(i+1)] for i in reversed(range(self.header['naxis']))])
        self.shapeout = (observation.ndetectors, observation.nfinesamples if finer_sampling else numpy.sum(observation.nsamples)) 

    def direct(self, map2d, reusein=False, copyout=True):
        map2d, shapeout = self._validate_input_direct(Map, map2d, reusein)
        self._validate_output_direct(Tod, shapeout)
        map2d.header = self.header
        tmf.pointing_matrix_direct(self.pmatrix, map2d.T, self._output_direct.T, self.npixels_per_sample)
        if not copyout: return self._output_direct
        return self._output_direct.copy('a')

    def transpose(self, signal, reusein=False, copyout=True):
        signal, shapein = self._validate_input_transpose(Tod, signal, reusein)
        self._validate_output_transpose(Map, shapein, header=self.header)
        tmf.pointing_matrix_transpose(self.pmatrix, signal.T, self._output_transpose.T,  self.npixels_per_sample)
        if not copyout: return self._output_transpose
        return self._output_transpose.copy('a')

    def __str__(self):
        return super(Projection, self).__str__()


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

    def direct(self, signal, reusein=False, copyout=True):
        signal, shapeout = self._validate_input_direct(Tod, signal, reusein)
        self._validate_output_direct(Tod, shapeout)
        self._output_direct.nsamples = tuple(numpy.divide(signal.nsamples, self.fine_sampling_factor))
        tmf.pacs_multiplexing_direct(signal.T, self._output_direct.T, self.fine_sampling_factor, self.ij)
        if not copyout: return self._output_direct
        return self._output_direct.copy('a')

    def transpose(self, signal, reusein=False, copyout=True):
        signal, shapein = self._validate_input_transpose(Tod, signal, reusein)
        self._validate_output_transpose(Tod, shapein)
        self._output_transpose.nsamples = tuple(numpy.multiply(signal.nsamples, self.fine_sampling_factor))
        tmf.pacs_multiplexing_transpose(signal.T, self._output_transpose.T, self.fine_sampling_factor, self.ij)
        if not copyout: return self._output_transpose
        return self._output_transpose.copy('a')

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        super(PacsMultiplexing, self)._validate_shape_direct(shapein)
        if shapein[1] % self.fine_sampling_factor != 0:
            raise ValidationError('The input timeline size ('+str(shapein[1])+') is not an integer times the fine sampling factor ('+str(self.fine_sampling_factor)+').')
        shapeout = list(shapein)
        shapeout[1] = shapeout[1] / self.fine_sampling_factor
        return tuple(shapeout)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return
        super(PacsMultiplexing, self)._validate_shape_transpose(shapeout)
        shapein = list(shapeout)
        shapein[1] = shapein[1] * self.fine_sampling_factor
        return tuple(shapein)


#-------------------------------------------------------------------------------


class Compression(AcquisitionModel):
    """
    Superclass for compressing the input signal.
    Author: P. Chanial
    """
    def __init__(self, compression_direct, compression_transpose, compression_factor, description):
        AcquisitionModel.__init__(self, description)
        self._compression_direct = compression_direct
        self._compression_transpose = compression_transpose
        self.compression_factor = compression_factor

    def direct(self, signal, reusein=False, copyout=True):
        if self.compression_factor == 1:
            if copyout: return signal.copy('a')
            return signal
        signal, shapeout = self._validate_input_direct(Tod, signal, reusein)
        self._validate_output_direct(Tod, shapeout)
        self._output_direct.nsamples = tuple(numpy.divide(signal.nsamples, self.compression_factor))
        self._compression_direct(signal, self._output_direct, self.compression_factor)
        if not copyout: return self._output_direct
        return self._output_direct.copy('a')

    def transpose(self, compressed, reusein=False, copyout=True):
        if self.compression_factor == 1:
            if copyout: return compressed.copy('a')
            return compressed
        compressed, shapein = self._validate_input_transpose(Tod, compressed, reusein)
        self._validate_output_transpose(Tod, shapein)
        self._output_transpose.nsamples = tuple(numpy.multiply(compressed.nsamples, self.compression_factor))
        self._compression_transpose(compressed, self._output_transpose, self.compression_factor)
        if not copyout: return self._output_transpose
        return self._output_transpose.copy('a')

    def _validate_shape_direct(self, shapein):
        if shapein is None:
            return None
        super(Compression, self)._validate_shape_direct(shapein)
        if shapein[1] % self.compression_factor != 0:
            raise ValidationError('The input timeline size ('+str(shapein[1])+') is not an integer times the compression factor ('+str(self.compression_factor)+').')
        shapeout = list(shapein)
        shapeout[1] /= self.compression_factor
        return tuple(shapeout)

    def _validate_shape_transpose(self, shapeout):
        if shapeout is None:
            return
        super(Compression, self)._validate_shape_transpose(shapeout)
        shapein = list(shapeout)
        shapein[1] *= self.compression_factor
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
        Compression.__init__(self, self.compression_average_direct, self.compression_average_transpose, compression_factor, description)

    def compression_average_direct(self, signal, output, compression_factor):
        tmf.compression_average_direct(signal.T, output.T, compression_factor)
         
    def compression_average_transpose(self, signal, output, compression_factor):
        tmf.compression_average_transpose(signal.T, output.T, compression_factor)
         

#-------------------------------------------------------------------------------


class DownSampling(Compression):
    """
    Downsample the input signal by picking up one sample out of a number specified by the compression factor
    Author: P. Chanial
    """
    def __init__(self, compression_factor, description=None):
        Compression.__init__(self, tmf.downsampling_direct, tmf.downsampling_transpose, compression_factor, description)
       

#-------------------------------------------------------------------------------


class Identity(AcquisitionModel):
    """
    Do nothing class.
    Author: P. Chanial
    """
    def __init__(self, description=None):
        AcquisitionModel.__init__(self, description)

    def direct(self, signal, reusein=False, copyout=True):
        self._validate_input_direct(numpy.ndarray, signal, False)
        if copyout: return signal.copy('a')
        return signal

    def transpose(self, signal, reusein=False, copyout=True):
        self._validate_input_transpose(numpy.ndarray, signal, False)
        if copyout: return signal.copy('a')
        return signal
        

#-------------------------------------------------------------------------------


class Masking(AcquisitionModel):
    """
    Apply a mask: where the mask is non-null or True, the input values are set to zero
    Author: P. Chanial
    """
    def __init__(self, mask, description=None):
        AcquisitionModel.__init__(self, description)
        if mask is None:
            self.mask = None
            self.shapein = None
            self.shapeout = None
        else:
            self.mask = mask
            self.shapein = mask.shape
            self.shapeout = mask.shape

    def direct(self, signal, reusein=False, copyout=True):
        self._validate_input_direct(numpy.ndarray, signal, False)
        if copyout: signal = signal.copy('a')
        status = tmf.masking(signal.T, self.mask.T)
        if status != 0: raise RuntimeError()
        return signal

    def transpose(self, signal, reusein=False, copyout=True):
        return self.direct(signal, reusein=reusein, copyout=copyout)
   
    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = None
        elif not isinstance(mask, numpy.ndarray):
            raise TypeError('Incorrect type for the input mask ('+str(type(mask))+').')
        elif numpy.rank(mask) == 0:
            raise ValueError('The input mask should not be scalar.')
        elif self.shapein is not None and self.shapein != mask.shapein:
            raise ValueError('The input mask has an incompatible shape '+str(mask.shape)+' instead of '+str(self.shapein)+'.')
        else:
            self._mask = numpy.array(mask, dtype='int8', copy=False)


#-------------------------------------------------------------------------------


class Packing(AcquisitionModel):
    """
    Convert 2d map into 1d map, under the control of a mask (true means observed)
    Author: P. Chanial
    """
    def __init__(self, mask, field=numpy.nan, description=None):
        AcquisitionModel.__init__(self, description)
        self.mask = numpy.array(mask, dtype='int8', copy=False)
        self.field = field
        self.shapein  = tuple(self.mask.shape)
        self.shapeout = (numpy.sum(self.mask == 0),)

    def direct(self, unpacked, reusein=False, copyout=True):
        unpacked, shapeout = self._validate_input_direct(Map, unpacked, reusein)
        self._validate_output_direct(numpy.ndarray, shapeout)
        tmf.unpack_transpose(unpacked.T, self.mask.T, self._output_direct.T)
        if copyout: return self._output_direct.copy('a')
        return self._output_direct

    def transpose(self, packed, reusein=False, copyout=True):
        packed, shapein = self._validate_input_transpose(Map, packed, reusein)
        self._validate_output_transpose(numpy.ndarray, shapein)
        tmf.unpack_direct(packed.T, self.mask.T, self._output_transpose.T, self.field)
        if copyout: return self._output_transpose.copy('a')
        return self._output_transpose


#-------------------------------------------------------------------------------


class Unpacking(AcquisitionModel):
    """
    Convert 1d map into 2d map, under the control of a mask (true means non-observed)
    Author: P. Chanial
    """
    def __init__(self, mask, field=0., description=None):
        AcquisitionModel.__init__(self, description)
        self.mask = numpy.array(mask, dtype='int8', copy=False)
        self.field = field
        self.shapein  = (numpy.sum(self.mask == 0),)
        self.shapeout = tuple(self.mask.shape)

    def direct(self, packed, reusein=False, copyout=True):
        packed, shapeout = self._validate_input_direct(numpy.ndarray, packed, reusein)
        self._validate_output_direct(Map, shapeout)
        tmf.unpack_direct(packed.T, self.mask.T, self._output_direct.T, self.field)
        if copyout: return self._output_direct.copy('a')
        return self._output_direct

    def transpose(self, unpacked, reusein=False, copyout=True):
        unpacked, shapein = self._validate_input_transpose(Map, unpacked, reusein)
        self._validate_output_transpose(numpy.ndarray, shapein)
        tmf.unpack_transpose(unpacked.T, self.mask.T, self._output_transpose.T)
        if copyout: return self._output_transpose.copy('a')
        return self._output_transpose


#-------------------------------------------------------------------------------


class Reshaping(AcquisitionModel):
    """
    Reshape arrays
    Author: P. Chanial
    """
    def __init__(self, shapein, shapeout, description=None):
        AcquisitionModel.__init__(self, description)
        self.shapein  = shapein
        self.shapeout = shapeout

    def direct(self, array, reusein=False, copyout=True):
        self._validate_input_direct(numpy.ndarray, array, False)
        output = array.reshape(self.shapeout)
        if copyout: return output.copy('a')
        return output

    def transpose(self, array, reusein=False, copyout=True):
        self._validate_input_transpose(numpy.ndarray, array, False)
        output = array.reshape(self.shapein)
        if copyout: return output.copy('a')
        return output



#-------------------------------------------------------------------------------
#
# PACS observations & simulations
#
#-------------------------------------------------------------------------------



class _Pacs():

    """
    Base class which encapsulates handy information about the PACS instrument and the processed 
    observation. See subclasses PacsObservation and PacsSimulation for details.
    - observing_mode      : 'prime', 'parallel' or 'transparent'
    - fine_sampling_factor : number of frames to be processed during the 0.025s period.
    - compression_factor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    Author: P. Chanial
    """

    sampling = 0.025 # 40 Hz sampling xXX not exact value! rather 0.024996

    def __init__(self, array, pointing_time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_detector_mask, keep_bad_detectors):

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
        
        if bad_detector_mask is None:
            bad_detector_maskFile = join(getenv('TAMASIS_DIR'),'data','PCalPhotometer_BadPixelMask_FM_v3.fits')
            bad_detector_mask = numpy.array(fitsopen(bad_detector_maskFile)[self.array].data, dtype='int8')
        else:
            array_shapes = {'blue':(32,64), 'green':(32,64), 'red':(16,32)}
            if bad_detector_mask.shape != array_shapes[self.array]:
                raise ValueError('Input bad pixel mask has incorrect shape '+str(bad_detector_mask.shape)+' instead of '+str(array_shapes[self.array]))
            bad_detector_mask = numpy.array(bad_detector_mask, dtype='int8', copy=False)
        self.bad_detector_mask = bad_detector_mask
        if self.transparent_mode:
            self.bad_detector_mask[:,0:16] = True
            self.bad_detector_mask[16:32,16:32] = True
            self.bad_detector_mask[:,32:] = True
        self.keep_bad_detectors = keep_bad_detectors

        self.ndetectors = self.bad_detector_mask.size
        if not self.keep_bad_detectors:
            self.ndetectors -= int(numpy.sum(self.bad_detector_mask))
        
        #XXX NOT IMPLEMENTED!
        self.ij = None


#-------------------------------------------------------------------------------


class PacsObservation(_Pacs):
    """
    Class which encapsulates handy information about the PACS instrument and the processed 
    observation.
    It is assumed that the detectors have sharp edges, so that the integration is simply equal 
    to the sum of the sky pixels values weighted by their intersections with the detector surface.
    It contains the following attributes:
    - filename           : name of the file name, including the array colour, but excluding 
                           the type. Example: '1342184520_blue'.
    - header             : pyfits header of the sky map
    - resolution         : resolution of the sky map, in arcsec
    - npixels_per_sample : number of sky pixels which intersect a PACS detector
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - bad_detector_mask  : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, filename, header=None, resolution=None, npixels_per_sample=None, fine_sampling_factor=1, bad_detector_mask=None, keep_bad_detectors=False):

        import pyfits

        if header is not None and resolution is not None:
            raise ValueError('It is not possible to specify a FITS header and an image resolution.')

        filename_, nfilenames = self._files2tmf(filename)

        channel, status = tmf.pacs_info_channel(filename_, nfilenames)
        if status != 0: raise RuntimeError()

        if npixels_per_sample is None:
            npixels_per_sample = 11 if channel == 'r' else 5

        use_mask = bad_detector_mask is not None
        if use_mask:
            if channel == 'r':
                if bad_detector_mask.shape != (16,32):
                    raise ValueError('Invalid shape of the input red bad detector mask: '+str(bad_detector_mask.shape)+'.')
            else:
                if bad_detector_mask.shape != (32,64):
                    raise ValueError('Invalid shape of the input blue bad detector mask: '+str(bad_detector_mask.shape)+'.')
            bad_detector_mask = numpy.array(bad_detector_mask, dtype='int8', copy=False)

        else:
            if channel == 'r':
                bad_detector_mask = numpy.zeros((16,32), dtype='int8')
            else:
                bad_detector_mask = numpy.zeros((32,64), dtype='int8')

        # retrieve information from the observations
        ndetectors, bad_detector_mask, transparent_mode, compression_factor, nsamples, status = tmf.pacs_info(filename_, nfilenames, fine_sampling_factor, keep_bad_detectors, use_mask, numpy.asfortranarray(bad_detector_mask))
        if status != 0: raise RuntimeError()

        if resolution is None:
            resolution = 3.
        if header is None:
            hstr, status = tmf.pacs_map_header(filename_, nfilenames, True, fine_sampling_factor, keep_bad_detectors, use_mask, bad_detector_mask, resolution)
            if status != 0: raise RuntimeError()
            header = str2fitsheader(hstr)

        self.filename = filename
        
        self.channel = channel
        self.nobservations = nfilenames
        self.npixels_per_sample = npixels_per_sample
        self.nsamples = tuple(nsamples)
        self.nfinesamples = numpy.sum(nsamples * compression_factor) * fine_sampling_factor
        self.ndetectors = ndetectors
        self.bad_detector_mask = numpy.ascontiguousarray(bad_detector_mask)
        self.keep_bad_detectors = keep_bad_detectors
        self.fine_sampling_factor = fine_sampling_factor
        self.transparent_mode = transparent_mode
        self.compression_factor = compression_factor
        self.header = header

         
        #_Pacs.__init__(self, channel, npixels, nsamples, ndetectors_per_sample, fine_sampling_factor, compression_factor, bad_detector_mask, keep_bad_detectors)

    def get_tod(self, do_flatfielding=True, do_subtraction_mean=True):
        """
        Returns the signal and mask timelines.
        """
        filename_, nfilenames = self._files2tmf(self.filename)
        signal, mask, status = tmf.pacs_timeline(filename_, self.nobservations, numpy.sum(self.nsamples), self.ndetectors, self.keep_bad_detectors, numpy.asfortranarray(self.bad_detector_mask), do_flatfielding, do_subtraction_mean)
        if status != 0: raise RuntimeError()
        
        return Tod(signal.T, mask=mask.T, nsamples=self.nsamples)

    def get_pointing_matrix(self, finer_sampling=True):
        nsamples = self.nfinesamples if finer_sampling else numpy.sum(self.nsamples)
        sizeofpmatrix = self.npixels_per_sample * nsamples * self.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        filename_, nfilenames = self._files2tmf(self.filename)
        status = tmf.pacs_pointing_matrix_filename(filename_, self.nobservations, finer_sampling, self.fine_sampling_factor, self.npixels_per_sample, nsamples, self.ndetectors, self.keep_bad_detectors, numpy.asfortranarray(self.bad_detector_mask), str(self.header).replace('\n', ''), pmatrix)
        if status != 0: raise RuntimeError()
        return pmatrix

    @staticmethod
    def _files2tmf(filename):
        if type(filename).__name__ == 'str':
            return filename, 1

        nfilenames = len(filename)
        length = max(len(f) for f in filename)
        filename_ = ''
        for f in filename:
            filename_ += f + (length-len(f))*' '
        return filename_, nfilenames


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
    - bad_detector_mask  : (nx,ny) mask of int8 values (0 or 1). 1 means dead pixel.
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, inputmap, time=None, ra=None, pa=None, chop=None, array='blue', npixels_per_sample=9, observing_mode='prime', fine_sampling_factor=1, compression_factor=None, bad_detector_mask=None, keep_bad_detectors=False):

        if not isinstance(inputmap, Map):
            raise TypeError('The input is not a Map.')

        if inputmap.header is None:
            raise TypeError('The input map header is not known.')

        fine_time = numpy.arange(time[0], floor((time[-1] - time[0]) / pacs.sampling) * fine_sampling_factor, _Pacs.sampling/fine_sampling_factor)

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fine_time, npixels_per_sample, observing_mode, fine_sampling_factor, compression_factor, bad_detector_mask, keep_bad_detectors)
 
        self.inputmap
        self.header = inputmap.header

    def get_tod(self, model, noise=None):
        """
        Returns simulated signal and mask (=None) timelines.
        """
        signal = model.direct(self.inputmap)
        if self.keep_bad_detectors:
            mask = numpy.zeros(signal.shape, dtype=bool)
            for idetector, badpixel in enumerate(self.bad_detector_mask.flat):
                if badpixel != 0:
                    mask[idetector,:] = True
        else:
            mask = None

        return Tod(signal, mask=mask, nsamples=self.nsamples)

    

#-------------------------------------------------------------------------------
#
# MADmap1 observations
#
#-------------------------------------------------------------------------------


class MadMap1Observation(object):
    """Class for the handling of an observation in the MADMAP1 format"""
    def __init__(self, todfile, invnttfile, mapmaskfile, convert, ndetectors):
        import pyfits
        import re
        self.npixels_per_sample, self.nfinesamples, status = tmf.read_madmap1_info(todfile, convert, ndetectors)
        if (status != 0): raise RuntimeError()

        self.todfile = todfile
        self.invnttfile = invnttfile
        m=re.search(r'(?P<filename>.*)\[(?P<extname>\w+)\]$', mapmaskfile)
        if m is None:
            raise ValueError('Invalid filename of the map mask.')
        filename = m.group('filename')
        extname  = m.group('extname')
        coverage = pyfits.fitsopen(filename)[extname].data
        self.mapmask = numpy.array(numpy.isnan(coverage), dtype='int8')
        self.convert = convert
        self.ndetectors = ndetectors
        self.nsamples = (self.nfinesamples,)
        header = pyfits.Header()
        header.update('simple', True)
        header.update('bitpix', -64)
        header.update('extend', True)
        header.update('naxis', 1)
        header.update('naxis1', numpy.sum(self.mapmask == 0))
        self.header = header

    def get_pointing_matrix(self, finer_sampling=False):
        tod = Tod.zeros((self.ndetectors,self.nfinesamples))
        sizeofpmatrix = self.npixels_per_sample * self.nfinesamples * self.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.read_madmap1(self.todfile, self.invnttfile, self.convert, self.npixels_per_sample, tod.T, pmatrix)
        if (status != 0): raise RuntimeError()
        return pmatrix

    def get_tod(self):
        tod = Tod.zeros((self.ndetectors,self.nfinesamples))
        sizeofpmatrix = self.npixels_per_sample * self.nfinesamples * self.ndetectors
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        status = tmf.read_madmap1(self.todfile, self.invnttfile, self.convert, self.npixels_per_sample, tod.T, pmatrix)
        if (status != 0): raise RuntimeError()
        return tod


#-------------------------------------------------------------------------------
#
# Map & Tod class
#
#-------------------------------------------------------------------------------



class FitsArray(numpy.ndarray):

    def __new__(cls, data, dtype=None, order='C', header=None, copy=True):
        from copy import copy as cp
        result = numpy.array(data, dtype=dtype, order=order, copy=copy, subok=True)
        if not isinstance(result, FitsArray):
            result = result.view(cls)
        if header is not None:
            result.header = header
        elif isinstance(data, FitsArray):
            result.header = cp(data.header) if copy else data.header
        else:
            result.header = None
        return result

    def __array_finalize__(self, obj):
        if obj is None: return
        self._header = getattr(obj, '_header', None)

    def __array_wrap__(self, obj, context=None):
        result = numpy.ndarray.__array_wrap__(self, obj, context).view(type(self))
        result._header = self._header
        return result

    @staticmethod
    def empty(shape, dtype='float64', order='C', header=None):
        return FitsArray(numpy.empty(shape, dtype, order), header=header, copy=False)

    @staticmethod
    def ones(shape, dtype='float64', order='C', header=None):
        return FitsArray(numpy.ones(shape, dtype, order), header=header, copy=False)

    @staticmethod
    def zeros(shape, dtype='float64', order='C', header=None):
        return FitsArray(numpy.zeros(shape, dtype, order), header=header, copy=False)

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        import pyfits
        if header is not None and not isinstance(header, pyfits.Header):
            raise TypeError('Incorrect type for the input header ('+str(type(header))+').') 
        self._header = header

    def tofile(self, fid, sep='', format='%s'):
        super(FitsArray,self).tofile(fid, sep, format)

    def writefits(self, filename):
        """Save a FitsArray instance to a fits file given a filename
        
        If the same file already exist it overwrites it.
        It should correspond to the input required by the PACS simulator
        """
        from pyfits import PrimaryHDU
        from os.path import isfile
        from os import remove
   
        hdu = PrimaryHDU(self.T if self.flags.f_contiguous else self, self.header)
        try:
            remove(filename)
        except OSError:
            pass
        try:
            hdu.writeto(filename)
        except IOError:
            pass


#-------------------------------------------------------------------------------


class Map(FitsArray):

    def __new__(cls, data, dtype=None, order='C', header=None, copy=True):
        result = FitsArray(data, dtype=dtype, order=order, copy=copy, header=header)
        if not isinstance(data,Map): result = result.view(cls)
        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None: return

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.header = self.header
        return result

    @staticmethod
    def empty(shape, dtype='float64', order='C', header=None):
        return Map(FitsArray.empty(shape, dtype, order, header), copy=False)

    @staticmethod
    def ones(shape, dtype='float64', order='C', header=None):
        return Map(FitsArray.ones(shape, dtype, order, header), copy=False)

    @staticmethod
    def zeros(shape, dtype='float64', order='C', header=None):
        return Map(FitsArray.zeros(shape, dtype, order, header), copy=False)

    def copy(self, order='C'):
        return Map(self, order=order, copy=True)

    def imshow(self, num=None, axis=True, title=None):
        """A simple graphical display function for the Map class"""
        from matplotlib.pyplot import gray, figure, imshow, colorbar, \
            draw, show, xlabel, ylabel, title as pyplottitle

        #gray()
        figure(num=num)

        # HACK around bug in imshow !!!
        imshow(self.T if self.flags.f_contiguous else self, 
               interpolation='nearest', 
               origin='lower')
        xlabel('Right Ascension')
        ylabel("Declination")
        if title is not None:
            pyplottitle(title)
        colorbar()
        draw()


#-------------------------------------------------------------------------------


class Tod(FitsArray):

    def __new__(cls, data, dtype=None, order='C', header=None, mask=None, nsamples=None, copy=False):
        result = FitsArray(data, dtype=dtype, order=order, copy=copy, header=header)
        if not isinstance(data, Tod):
            result = result.view(cls)
        if mask is None and isinstance(data, Tod):
            mask = data.mask
        if mask is not None:
            mask = numpy.array(mask, dtype='int8', copy=copy, order=order)
        result.mask = mask
        if nsamples is None:
            if isinstance(data, Tod):
                result.nsamples = data.nsamples
            else:
                result.nsamples = (data.shape[-1],)
        else:
            if numpy.sum(nsamples) != data.shape[-1]:
                raise ValueError('The sum of the slice sizes is not equal to the number of samples in the input argument).')
            result.nsamples = tuple(nsamples)
        
        return result

    def __array_finalize__(self, obj):
        FitsArray.__array_finalize__(self, obj)
        if obj is None: return
        self.mask = getattr(obj, 'mask', None)
        self.nsamples = getattr(obj, 'nsamples', None)

    def __array_wrap__(self, obj, context=None):
        result = FitsArray.__array_wrap__(self, obj, context=context).view(type(self))
        result.mask = self.mask
        result.nsamples = self.nsamples
        return result

    @staticmethod
    def empty(shape, dtype='float64', order='C', header=None, mask=None, nsamples=None):
        return Tod(FitsArray.empty(shape, dtype, order, header), mask=mask, nsamples=nsamples, copy=False)

    @staticmethod
    def ones(shape, dtype='float64', order='C', header=None, mask=None, nsamples=None):
        return Tod(FitsArray.ones(shape, dtype, order, header), mask=mask, nsamples=nsamples, copy=False)

    @staticmethod
    def zeros(shape, dtype='float64', order='C', header=None, mask=None, nsamples=None):
        return Tod(FitsArray.zeros(shape, dtype, order, header), mask=mask, nsamples=nsamples, copy=False)
    
    def copy(self, order='C'):
        return Tod(self, order=order, copy=True)

    def imshow(self, num=None, axis=True, title=None):
        """A simple graphical display function for the Map class"""
        from matplotlib.pyplot import gray, figure, imshow, colorbar, \
            draw, show, xlabel, ylabel, title as pyplottitle

        gray()
        figure(num=num)
        imshow(self, 
               aspect='auto', 
               interpolation='nearest', 
               origin='lower')
        xlabel('Detector number')
        ylabel("Signal")
        if title is not None:
            pyplottitle(title)
        colorbar()
        draw()

    def __str__(self):
        nslices = len(self.nsamples)
        output = 'Tod ['+str(self.shape[-2])+' detector'
        if self.shape[-2] > 1:
            output += 's'
        strsamples = ','.join((str(self.nsamples[i]) for i in range(nslices)))
        if nslices > 1:
            output += ', (' + strsamples + ') samples in ' + str(nslices) + ' slices'
        else:
            output += ', ' + strsamples + ' sample'
            if self.shape[-1] > 1:
                output += 's'
        return output + ']'



#-------------------------------------------------------------------------------
#
# Linear operator
#
#-------------------------------------------------------------------------------


class LeastSquareMatvec(object):
    def __init__(self, model, unpacking=None):
        if unpacking is None:
            self.unpacking = Reshaping(numpy.product(model.shapein), model.shapein)
        else:
            self.unpacking = unpacking
        self.model = model * unpacking
    def __call__(self, xvec):
        xout = self.model.transpose(self.model.direct(Map(xvec, copy=False), copyout=False))
        return xout

class PcgCallback():
    def __init__(self):
        self.count = 1
    def __call__(self, x):
        print 'PCG Iteration '+str(self.count)
        self.count += 1

class RegularizedLeastSquareMatvec(LeastSquareMatvec):
    def __init__(self, model, unpacking, hyper):
        LeastSquareMatvec.__init__(self, model, unpacking)
        # to enforce 2d hyper
        if numpy.size(hyper) == 1:
            hyper = 2 * (hyper,)
        self.hyper = hyper
        # define smooth prior
        self.prior_direct = [Smooth(i) for i in xrange(len(hyper))]
        self.prior_transpose = [SmoothT(i) for i in xrange(len(hyper))]
    def __call__(self, xvec):
        # likelihood
        xout = super(RegularizedLeastSquareMatvec, self).__call__(xvec)
        x = self.unpacking.direct(Map(xvec, copy=False), copyout=False)
        # add prior component
        for i in xrange(len(self.hyper)):
            dTd = Map(self.prior_transpose[i](self.prior_direct[i](x)), copy=False)
            dTd *= self.hyper[i]
            xout += self.unpacking.transpose(dTd, copyout=False)
        return xout

# Priors
class Smooth():
    """Define smoothness priors along specific axis"""
    def __init__(self, axis):
        self.axis = axis
    def __call__(self, arr):
        return numpy.diff(arr, axis=self.axis)

class SmoothT(Smooth):
    """Define smoothness priors along specific axis
    This is the transpose model
    """
    def __call__(self, arr):
        return diffT(arr, axis=self.axis)

def diffT(arr, axis=-1):
    """transpose of diff for n=1"""
    shape = list(arr.shape)
    shape[axis] = 1
    border = numpy.zeros(shape)
    out = numpy.concatenate((border, arr, border), axis=axis)
    return - numpy.diff(out, axis=axis)



#-------------------------------------------------------------------------------
#
# Miscellaneous routines
#
#-------------------------------------------------------------------------------


def naive_mapper(tod, model, weights=False):
    """
    Returns a naive map, i.e.: model.transpose(tod) / model.transpose(1)
    """
    mymap = model.transpose(tod)
    unity = tod.copy()
    unity[:] = 1.
    map_weights = model.transpose(unity)
    mymap /= map_weights
    if weights:
        return mymap, map_weights
    return mymap

 
#-------------------------------------------------------------------------------


def deglitch_l2std(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the standard deviation to the mean.
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    tmf.deglitch_l2b_std(projection.pmatrix, nx, ny, tod.T, tod.mask.T, nsigma, npixels_per_sample)

 
#-------------------------------------------------------------------------------


def deglitch_l2mad(tod, projection, nsigma=5.):
    """
    Second level deglitching. Each frame is projected onto the sky. The values associated
    to a sky map pixel are then sigma clipped, using the MAD (median absolute deviation to the median)
    """
    nx = projection.header['naxis1']
    ny = projection.header['naxis2']
    npixels_per_sample = projection.npixels_per_sample
    tmf.deglitch_l2b_mad(projection.pmatrix, nx, ny, tod.T, tod.mask.T, nsigma, npixels_per_sample)

 
#-------------------------------------------------------------------------------


def filter_median(tod, length=10):
    """
    Median filtering, O(1) in window length
    """
    filtered = tod.copy()
    status = tmf.filter_median(filtered.T, length, tod.nsamples)
    if status != 0:
        raise RuntimeError()
    return filtered

 
#-------------------------------------------------------------------------------


def str2fitsheader(string):
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


def any_neq(a,b, precision):
    mask = numpy.isnan(a)
    if numpy.any(mask != numpy.isnan(b)):
        return True
    return numpy.any((abs(a-b) > 10.**(-precision) * abs(a))[numpy.isnan(a).__neg__()])
