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
    def __init__(self):
        pass

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

    #abstractmethod
    def validateDirect(self, shapeIn):
        self.shapeOut = shapeIn
        for model in self:
            self.shapeOut = model.validateDirect(self.shapeOut)

    #abstractmethod
    def validateTranspose(self, shapeOut):
        self.shapeIn = shapeOut
        for model in reversed(self):
            self.shapeIn = model.validateTranspose(self.shapeIn)

    def validate(self):
        self.validateDirect(None)
        self.validateTranspose(None)

    def _validateShapeDirect(self, shapeIn):
        if shapeIn is None or self.shapeIn is None:
            return
        if shapeIn != self.shapeIn:
            raise ValidationError('The input of '+self.__class__.__name__+' has incompatible shape '+str(shapeIn)+' instead of '+str(self.shapeIn)+'.')

    def _validateShapeTranspose(self, shapeOut):
        if shapeOut is None or self.shapeOut is None:
            return
        if shapeOut != self.shapeOut:
            raise ValidationError('The input of '+self.__class__.__name__+' transpose has incompatible shape '+str(shapeOut)+' instead of '+str(self.shapeOut)+'.')

    def _ensureAllocatedDirect(self, shapeOut):
        if shapeOut is None:
            raise ValueError('The shape of the output of '+self.__class__.__name__+' is not known.')
        if self.outputDirect is not None and shapeOut == self.outputDirect.shape:
            return
        print 'Info: allocating '+str(numpy.product(shapeOut)/2.**17)+' MiB for the output of '+self.__class__.__name__
        self.outputDirect = numpy.zeros(shapeOut, dtype=numpy.float64, order='fortran')

    def _ensureAllocatedTranspose(self, shapeIn):
        if shapeIn is None:
            raise ValueError('The shape of the input of '+self.__class__.__name__+' is not known.')
        if self.outputTranspose is not None and shapeIn == self.outputTranspose.shape:
            return
        print 'Info: allocating '+str(numpy.product(shapeIn)/2.**17)+' MiB for the output of the transpose of '+self.__class__.__name__
        self.outputTranspose = numpy.zeros(shapeIn, dtype=numpy.float64, order='fortran')

    shapeIn         = None   # input is unconstrained
    shapeOut        = None   # output is unconstrained
    blocks          = []     # components are ordered as in the direct order
    outputDirect    = None   # stores the input of the transpose model whose memory allocation is re-used as the output of the direct model
    outputTranspose = None   # stores the input of the direct model whose memory allocation is re-used as the output of the transpose model

    def __mul__(self, other):
        newModel = AcquisitionModel()
        newModel.blocks = []
        for block in other:
            newModel.blocks.append(block)
        for block in self:
            newModel.blocks.append(block)
        newModel.validate()
        return newModel

    def __getitem__(self, index):
        if len(self) == 0:
            raise IndexError('This acquisition model has no component.')
        return self.blocks[index]

    def __setitem__(self, index, value):
        if len(self) == 0:
            raise IndexError('This acquisition model has no component.')

        try:
            oldvalue = self[index]
            oldshapeIn = self.shapeIn
            oldshapeOut = self.shapeOut
        except:
            raise IndexError('Only substitutions are allowed in an acquisition model.')

        self.blocks[index] = value
        try:
            self.validate()
        except ValidationError:
            self.blocks[index] = oldvalue
            self.shapeIn = oldshapeIn
            self.shapeOut = oldshapeOut
            raise ValidationError('Incompatible component.')

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
        result = self.__class__.__name__ + ' [input:'
        result += 'unconstrained' if self.shapeIn is None else str(self.shapeIn)
        result += ', output:'
        result += 'unconstrained' if self.shapeOut is None else str(self.shapeOut)
        result += ']'
        if len(self) == 0:
            return result
        return result+'\n  '+'\n  '.join((str(block) for block in reversed(self)))


#-------------------------------------------------------------------------------


class PacsProjectionSharpEdges(AcquisitionModel):
    """
    This class handles the integration of the sky inside the PACS detectors.
    It is assumed that the detectors have sharp edges, so that the integration is simply equal 
    to the sum of the sky pixels values weighted by their intersections with the detector surface.
    Author: P. Chanial
    """
    def __init__(self, pacs):
        AcquisitionModel.__init__(self)
        self.npixelsPerSample = pacs.npixelsPerSample
        sizeofpmatrix = pacs.npixelsPerSample * pacs.nfineSamples * pacs.ndetectors
        print 'Info: allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        self.pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        self.header = pacs.header
        self.shapeIn = (pacs.header["naxis1"], pacs.header["naxis2"])
        self.shapeOut = (pacs.nfineSamples, pacs.ndetectors) 
        tmmf.pacs_pointing_matrix_array(pacs.pointingTime, pacs.ra, pacs.dec, pacs.pa, pacs.chop, pacs.fineTime, pacs.npixelsPerSample, pacs.ndetectors, pacs.array, pacs.transparentMode, str(pacs.header).replace('\n', ''), self.pmatrix)

    def direct(self, map2d):
        self._ensureAllocatedDirect(self.validateDirect(map2d.shape))
        self.outputTranspose = map2d
        tmmf.pacs_projection_sharp_edges_direct(self.pmatrix, map2d, self.outputDirect, self.npixelsPerSample)
        return self.outputDirect

    def transpose(self, signal):
        self._ensureAllocatedTranspose(self.validateTranspose(signal.shape))
        self.outputDirect = signal
        tmmf.pacs_projection_sharp_edges_transpose(self.pmatrix, signal, self.outputTranspose,  self.npixelsPerSample)
        return self.outputTranspose

    def validateDirect(self, shapeIn):
        self._validateShapeDirect(shapeIn)
        return self.shapeOut

    def validateTranspose(self, shapeOut):
        self._validateShapeTranspose(shapeOut)
	return self.shapeIn

    def __str__(self):
        return super(PacsProjectionSharpEdges, self).__str__()#+' => '+self.filename+' ['+str(self.first)+','+str(self.last)+']'


#-------------------------------------------------------------------------------


class PacsMultiplexing(AcquisitionModel):
    """
    Performs the multiplexing of the PACS subarrays. The subarray columns are read one after the
    other, in a 0.025s cycle (40Hz).
    Author: P. Chanial
    """
    def __init__(self, pacs):
        AcquisitionModel.__init__(self)
        self.fineSamplingFactor = pacs.fineSamplingFactor
        self.ij = pacs.ij

    def direct(self, signal):
        self._ensureAllocatedDirect(self.validateDirect(signal.shape))
        self.outputTranspose = signal
        tmmf.pacs_multiplexing_direct(signal, self.outputDirect, self.fineSamplingFactor, self.ij)
        return self.outputDirect

    def transpose(self, signal):
        self._ensureAllocatedTranspose(self.validateTranspose(signal.shape))
        self.outputDirect = signal
        tmmf.pacs_multiplexing_transpose(signal, self.outputTranspose, self.fineSamplingFactor, self.ij)
        return self.outputTranspose

    def validateDirect(self, shapeIn):
        if shapeIn is None:
            return None
        if shapeIn[0] % self.fineSamplingFactor != 0:
            raise ValidationError('The input timeline size ('+str(shapeIn[0])+') is not an integer times the fine sampling factor ('+str(self.fineSamplingFactor)+').')
        shapeOut = list(shapeIn)
        shapeOut[0] = shapeOut[0] / self.fineSamplingFactor
        return tuple(shapeOut)

    def validateTranspose(self, shapeOut):
        if shapeOut is None:
            return
        shapeIn = list(shapeOut)
        shapeIn[0] = shapeIn[0] * self.fineSamplingFactor
        return tuple(shapeIn)


#-------------------------------------------------------------------------------


class CompressionAverage(AcquisitionModel):
    """
    Compress the input signal by averaging blocks of specified size.
    Author: P. Chanial
    """
    def __init__(self, compressionFactor):
        AcquisitionModel.__init__(self)
        self.compressionFactor = compressionFactor

    def direct(self, signal):
        if self.compressionFactor == 1:
            return signal
        self._ensureAllocatedDirect(self.validateDirect(signal.shape))
        self.outputTranspose = signal
        tmmf.compression_average_direct(signal, self.outputDirect, self.compressionFactor)
        return self.outputDirect

    def transpose(self, compressed):
        if self.compressionFactor == 1:
            return compressed
        self._ensureAllocatedTranspose(self.validateTranspose(compressed.shape))
        self.outputDirect = compressed
        tmmf.compression_average_transpose(compressed, self.outputTranspose, self.compressionFactor)
        return self.outputTranspose

    def validateDirect(self, shapeIn):
        if shapeIn is None:
            return None
        if shapeIn[0] % self.compressionFactor != 0:
            raise ValidationError('The input timeline size ('+str(shapeIn[0])+') is not an integer times the compression factor ('+str(self.compressionFactor)+').')
        shapeOut = list(shapeIn)
        shapeOut[0] /= self.compressionFactor
        return tuple(shapeOut)

    def validateTranspose(self, shapeOut):
        if shapeOut is None:
            return
        shapeIn = list(shapeOut)
        shapeIn[0] *= self.compressionFactor
        return tuple(shapeIn)

    def __str__(self):
        return super(CompressionAverage, self).__str__()+' x'+str(self.compressionFactor)
        

#-------------------------------------------------------------------------------


class Identity(AcquisitionModel):
    """
    Do nothing class.
    Author: P. Chanial
    """
    def __init__(self):
        AcquisitionModel.__init__(self)

    def direct(self, signal):
        return signal

    def transpose(self, signal):
        return signal

    def validateDirect(self, shape):
        return shape

    def validateTranspose(self, shape):
        return shape
    



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

    def __init__(self, array, pointingTime, ra, dec, pa, chop, fineTime, npixelsPerSample, observingMode, fineSamplingFactor, compressionFactor, badPixelMask, keepBadPixels):

        from os import getenv
        from os.path import join
        from pyfits import fitsopen

        if array.lower() not in ('blue', 'green', 'red'):
            raise ValueError("The input array is not 'blue', 'green' nor 'red'.")

        if observingMode is not None:
            observingMode = observingMode.lower()

        if observingMode not in (None, 'prime', 'parallel', 'transparent'):
            raise ValueError("Observing mode is not 'prime', 'parallel' nor 'transparent'.")

        if compressionFactor is None:
            if observingMode is None:
                raise ValueError('The compression factor is not specified.')
            compressionFactor = {'prime':4, 'parallel':8, 'transparent':1}[observingMode]

        self.transparentMode = observingMode == 'transparent'

        if (fineSamplingFactor not in (1,2,4,8,16)):
            raise ValueError('Invalid fine sampling factor. It may be 1, 2, 4, 8 or 16')
 
        # 'blue', 'green' or 'red' array
        self.array = array.lower()

        # astrometry
        if pointingTime is None or ra is None or dec is None or pa is None or chop is None:
            raise ValueError('The simulated scan astrometry is not defined.')
        shape = pointingTime.shape
        if ra.shape != shape or dec.shape != shape or pa.shape != shape or chop.shape != shape:
            raise ValueError('The input time, ra, dec, pa, chope do not have the same shape.')
        self.pointingTime = pointingTime
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.chop = chop
        self.header = None

        # sampling information
        self.fineTime = fineTime
        self.nfineSamples = fineTime.size
        self.npixelsPerSample = npixelsPerSample
        self.fineSamplingFactor = fineSamplingFactor
        self.compressionFactor = compressionFactor
        
        if badPixelMask is None:
            badPixelMaskFile = join(getenv('TAMASIS_DIR'),'data','PCalPhotometer_BadPixelMask_FM_v3.fits')
            badPixelMask = fitsopen(badPixelMaskFile)[self.array].data
        else:
            arrayShapes = {'blue':(32,64), 'green':(32,64), 'red':(16,32)}
            if badPixelMask.shape != arrayShapes[self.array]:
                raise ValueError('Input bad pixel mask has incorrect shape '+str(badPixelMask.shape)+' instead of '+str(arrayShapes[self.array]))
        self.badPixelMask = badPixelMask.copy()
        if self.transparentMode:
            self.badPixelMask[:,0:16] = 1
            self.badPixelMask[16:32,16:32] = 1
            self.badPixelMask[:,32:] = 1
        self.keepBadPixels = keepBadPixels

        self.ndetectors = self.badPixelMask.size
        if not self.keepBadPixels:
            self.ndetectors -= int(numpy.sum(self.badPixelMask))
        
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
            cards.append(pyfits.Card().fromstring(line))
            if line[0:3] == 'END': break
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
    - npixelsPerSample   : number of sky pixels which intersect a PACS detector
    - observingMode      : 'prime', 'parallel' or 'transparent'
    - fineSamplingFactor : number of frames to be processed during the 0.025s period.
    - compressionFactor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - badPixelMask       : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - keepBadPixels      : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, filename, first=0, last=None, header=None, resolution=None, npixelsPerSample=9, fineSamplingFactor=1, badPixelMask=None, keepBadPixels=False):

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
        compressionFactor = int(compression + 0.5)
        if abs(compressionFactor - compression) > 1.001 * compression:
            raise ValueError('The file '+filename+'_Time.fits is not sampled an integer time of '+str(_Pacs.sampling)+'s.')
        if compressionFactor == 1:
            observingMode = 'transparent'
        elif compressionFactor == 4:
            observingMode = 'prime'
        elif compressionFactor == 8:
            observingMode = 'parallel'
        else:
            observingMode = None
        samplingFactor = fineSamplingFactor * compressionFactor
        fineTime = numpy.zeros((last - first) * samplingFactor, dtype='float64')
        for i in range(samplingFactor):
            fineTime[i::samplingFactor] = time + i * _Pacs.sampling / fineSamplingFactor

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fineTime, npixelsPerSample, observingMode, fineSamplingFactor, compressionFactor, badPixelMask, keepBadPixels)

        if resolution is None:
            resolution = 3.
        if header is None:
            header = self._str2fitsheader(tmmf.pacs_map_header_file(filename, self.transparentMode, first, last, self.fineSamplingFactor, self.compressionFactor, resolution))
        self.header = header

        self.filename = filename
        self.first = first
        self.last = last
        
        #XXX REMOVE ME
        self.ij2 = tmmf.pacs_info_ij(filename, True, self.ndetectors)
 
    def get_timeline(self):
        """
        Returns the signal and mask timelines as fortran arrays.
        """
        signal, mask = tmmf.pacs_timeline(self.filename, self.transparentMode, self.first+1, self.last, self.ndetectors)
        mask = numpy.array(mask, dtype=bool)
        return signal, mask
   

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
    - npixelsPerSample   : number of sky pixels which intersect a PACS detector
    - observingMode      : 'prime', 'parallel' or 'transparent'
    - fineSamplingFactor : number of frames to be processed during the 1/40Hz period.
    - compressionFactor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - badPixelMask       : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - keepBadPixels      : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, inputmap, time=None, ra=None, pa=None, chop=None, array='blue', npixelsPerSample=9, observingMode='prime', fineSamplingFactor=1, compressionFactor=None, badPixelMask=None, keepBadPixels=False):

        if not isinstance(inputmap, Map):
            raise TypeError('The input is not a Map.')

        if inputmap.header is None:
            raise TypeError('The input map header is not known.')

        fineTime = numpy.arange(time[0], floor((time[-1] - time[0]) / pacs.sampling) * fineSamplingFactor, _Pacs.sampling/fineSamplingFactor)

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fineTime, npixelsPerSample, observingMode, fineSamplingFactor, compressionFactor, badPixelMask, keepBadPixels)
 
        self.inputmap
        self.header = inputmap.header

    def get_timeline(self, model, noise=None):
        """
        Returns simulated signal and mask (=None) timelines.
        """
        return model.direct(self.inputmap), None

    


#-------------------------------------------------------------------------------
#
# Miscellaneous methods
#
#-------------------------------------------------------------------------------



def applymask(signal, mask):
    """
    Sets masked values of input signal to zero (in-place operation) and returns it.
    Author: P. Chanial
    """
    if mask is None:
        return signal
    if (signal.shape != mask.shape):
        raise ValueError("Signal and mask are not conformant.")
    signal[mask] = 0
    return signal



