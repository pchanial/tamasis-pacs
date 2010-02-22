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
        except ValidationError as inst:
            self.blocks[index] = oldvalue
            self.shapeIn = oldshapeIn
            self.shapeOut = oldshapeOut
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
        tmmf.pacs_pointing_matrix(pacs.array, pacs.pointingTime, pacs.ra, pacs.dec, pacs.pa, pacs.chop, pacs.fineTime, pacs.npixelsPerSample, pacs.ndetectors, pacs.transparentMode, pacs.badPixelMask.astype('int8'), pacs.keepBadPixels, str(pacs.header).replace('\n', ''), self.pmatrix)

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
            badPixelMask = numpy.array(fitsopen(badPixelMaskFile)[self.array].data, order='fortran')
        else:
            arrayShapes = {'blue':(32,64), 'green':(32,64), 'red':(16,32)}
            if badPixelMask.shape != arrayShapes[self.array]:
                raise ValueError('Input bad pixel mask has incorrect shape '+str(badPixelMask.shape)+' instead of '+str(arrayShapes[self.array]))
        self.badPixelMask = badPixelMask.copy('fortran').astype(bool)
        if self.transparentMode:
            self.badPixelMask[:,0:16] = True
            self.badPixelMask[16:32,16:32] = True
            self.badPixelMask[:,32:] = True
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
            header = self._str2fitsheader(tmmf.pacs_map_header(array, time, ra, dec, pa, chop, fineTime, self.transparentMode, self.badPixelMask.astype('int8'), keepBadPixels, resolution))

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
        signal, mask = tmmf.pacs_timeline(self.filename, self.array, self.first+1, self.last, self.ndetectors, self.transparentMode, self.badPixelMask, self.keepBadPixels)
        mask = numpy.array(mask, dtype=bool)
        if self.keepBadPixels:
           for idetector, badPixel in enumerate(self.badPixelMask.flat):
               if badPixel:
                   mask[:,idetector] = True
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
        signal = model.direct(self.inputmap)
        if self.keepBadPixels:
            mask = numpy.zeros(signal.shape, dtype=bool, order='fortran')
            for idetector, badPixel in enumerate(self.badPixelMask.flat):
                if badPixel:
                    mask[:,idetector] = True
        else:
            mask = None

        return signal, mask
    



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



