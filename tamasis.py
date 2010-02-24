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

    def validateChainDirect(self, shapeIn):
        self.shapeOut = shapeIn
        for model in self:
            self.shapeOut = model.validateShapeDirect(self.shapeOut)

    def validateChainTranspose(self, shapeOut):
        self.shapeIn = shapeOut
        for model in reversed(self):
            self.shapeIn = model.validateShapeTranspose(self.shapeIn)

    def validateChain(self):
        self.validateChainDirect(None)
        self.validateChainTranspose(None)

    def validateShapeDirect(self, shapeIn):
        if shapeIn is None or self.shapeIn is None:
            if self.shapeOut is not None:
               return self.shapeOut
            else:
               return shapeIn
        if shapeIn != self.shapeIn:
            raise ValidationError('The input of '+self.__class__.__name__+' has incompatible shape '+str(shapeIn)+' instead of '+str(self.shapeIn)+'.')
        return self.shapeOut

    def validateShapeTranspose(self, shapeOut):
        if shapeOut is None or self.shapeOut is None:
            if self.shapeIn is not None:
                return self.shapeIn
            else:
                return shapeOut
        if shapeOut != self.shapeOut:
            raise ValidationError('The input of '+self.__class__.__name__+' transpose has incompatible shape '+str(shapeOut)+' instead of '+str(self.shapeOut)+'.')
        return self.shapeIn

    def validateInputDirect(self, cls, data):
        if not isinstance(data, cls):
            raise TypeError("The input of '+self.__class__.__name+' has an invalid type '"+cls.__name__+"'.")
        return self.validateShapeDirect(data.shape)

    def validateInputTranspose(self, cls, data):
        if not isinstance(data, cls):
            raise TypeError("The input of '+self.__class__.__name+' transpose has an invalid type '"+cls.__name__+"'.")
        return self.validateShapeTranspose(data.shape)

    def validateOutputDirect(self, cls, shapeOut, **options):
        if shapeOut is None:
            raise ValueError('The shape of the output of '+self.__class__.__name__+' is not known.')
        if self.outputDirect is not None and shapeOut == self.outputDirect.shape:
            return
        print 'Info: allocating '+str(numpy.product(shapeOut)/2.**17)+' MiB for the output of '+self.__class__.__name__
        self.outputDirect = numpy.ndarray.__new__(cls, shapeOut, dtype=numpy.float64, order='fortran')
        self.outputDirect.__init__(options)

    def validateOutputTranspose(self, cls, shapeIn, **options):
        if shapeIn is None:
            raise ValueError('The shape of the input of '+self.__class__.__name__+' is not known.')
        if self.outputTranspose is not None and shapeIn == self.outputTranspose.shape:
            return
        print 'Info: allocating '+str(numpy.product(shapeIn)/2.**17)+' MiB for the output of the transpose of '+self.__class__.__name__
        self.outputTranspose = numpy.ndarray.__new__(cls, shapeIn, dtype=numpy.float64, order='fortran')
        self.outputTranspose.__init__(options)

    shapeIn         = None   # input is unconstrained
    shapeOut        = None   # output is unconstrained
    outputDirect    = None   # stores the input of the transpose model. Its memory allocation is re-used as the output of the direct model
    outputTranspose = None   # stores the input of the direct model. Its memory allocation is re-used as the output of the transpose model
    blocks          = []     # components are ordered as in the direct order

    def __mul__(self, other):
        newModel = AcquisitionModel(None)
        newModel.blocks = []
        for block in other:
            newModel.blocks.append(block)
        for block in self:
            newModel.blocks.append(block)
        newModel.validateChain()
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
        result = self.__class__.__name__
        if self.description is not None:
            result = self.description+' ('+self.__class__.__name__+')'
        else:
            result = self.__class__.__name__
        if self.shapeIn is not None or self.shapeOut is not None:
            result += ' [input:'
            result += 'unconstrained' if self.shapeIn is None else str(self.shapeIn).replace(' ','')
            result += ', output:'
            result += 'unconstrained' if self.shapeOut is None else str(self.shapeOut).replace(' ','')
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
        self.npixelsPerSample = pacs.npixelsPerSample
        sizeofpmatrix = pacs.npixelsPerSample * pacs.nfineSamples * pacs.ndetectors
        print 'Info: allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        self.pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        self.header = pacs.header
        self.shapeIn = (pacs.header["naxis1"], pacs.header["naxis2"])
        self.shapeOut = (pacs.nfineSamples, pacs.ndetectors) 
        tmmf.pacs_pointing_matrix(pacs.array, pacs.pointingTime, pacs.ra, pacs.dec, pacs.pa, pacs.chop, pacs.fineTime, pacs.npixelsPerSample, pacs.ndetectors, pacs.transparentMode, pacs.badPixelMask.astype('int8'), pacs.keepBadDetectors, str(pacs.header).replace('\n', ''), self.pmatrix)

    def direct(self, map2d):
        self.validateOutputDirect(Tod, self.validateInputDirect(Map, map2d))
        map2d.header = self.header
        self.outputTranspose = map2d
        tmmf.pacs_projection_sharp_edges_direct(self.pmatrix, map2d, self.outputDirect, self.npixelsPerSample)
        return self.outputDirect

    def transpose(self, signal):
        self.validateOutputTranspose(Map, self.validateInputTranspose(Tod, signal), header=self.header)
        self.outputDirect = signal
        tmmf.pacs_projection_sharp_edges_transpose(self.pmatrix, signal, self.outputTranspose,  self.npixelsPerSample)
        return self.outputTranspose

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
        self.fineSamplingFactor = pacs.fineSamplingFactor
        self.ij = pacs.ij

    def direct(self, signal):
        self.validateOutputDirect(Tod, self.validateInputDirect(Tod, signal))
        self.outputTranspose = signal
        tmmf.pacs_multiplexing_direct(signal, self.outputDirect, self.fineSamplingFactor, self.ij)
        return self.outputDirect

    def transpose(self, signal):
        self.validateOutputTranspose(Tod, self.validateInputTranspose(Tod, signal))
        self.outputDirect = signal
        tmmf.pacs_multiplexing_transpose(signal, self.outputTranspose, self.fineSamplingFactor, self.ij)
        return self.outputTranspose

    def validateShapeDirect(self, shapeIn):
        if shapeIn is None:
            return None
        super(PacsMultiplexing, self).validateShapeDirect(shapeIn)
        if shapeIn[0] % self.fineSamplingFactor != 0:
            raise ValidationError('The input timeline size ('+str(shapeIn[0])+') is not an integer times the fine sampling factor ('+str(self.fineSamplingFactor)+').')
        shapeOut = list(shapeIn)
        shapeOut[0] = shapeOut[0] / self.fineSamplingFactor
        return tuple(shapeOut)

    def validateShapeTranspose(self, shapeOut):
        if shapeOut is None:
            return
        super(PacsMultiplexing, self).validateShapeTranspose(shapeOut)
        shapeIn = list(shapeOut)
        shapeIn[0] = shapeIn[0] * self.fineSamplingFactor
        return tuple(shapeIn)


#-------------------------------------------------------------------------------


class Compression(AcquisitionModel):
    """
    Superclass for compressing the input signal.
    Author: P. Chanial
    """
    def __init__(self, compressionDirect, compressionTranspose, compressionFactor, description):
        AcquisitionModel.__init__(self, description)
        self.compressionDirect = compressionDirect
        self.compressionTranspose = compressionTranspose
        self.compressionFactor = compressionFactor

    def direct(self, signal):
        if self.compressionFactor == 1:
            return signal
        self.validateOutputDirect(Tod, self.validateInputDirect(Tod, signal))
        self.outputTranspose = signal
        self.compressionDirect(signal, self.outputDirect, self.compressionFactor)
        return self.outputDirect

    def transpose(self, compressed):
        if self.compressionFactor == 1:
            return compressed
        self.validateOutputTranspose(Tod, self.validateInputTranspose(Tod, compressed))
        self.outputDirect = compressed
        self.compressionTranspose(compressed, self.outputTranspose, self.compressionFactor)
        return self.outputTranspose

    def validateShapeDirect(self, shapeIn):
        if shapeIn is None:
            return None
        super(Compression, self).validateShapeDirect(shapeIn)
        if shapeIn[0] % self.compressionFactor != 0:
            raise ValidationError('The input timeline size ('+str(shapeIn[0])+') is not an integer times the compression factor ('+str(self.compressionFactor)+').')
        shapeOut = list(shapeIn)
        shapeOut[0] /= self.compressionFactor
        return tuple(shapeOut)

    def validateShapeTranspose(self, shapeOut):
        if shapeOut is None:
            return
        super(Compression, self).validateShapeTranspose(shapeOut)
        shapeIn = list(shapeOut)
        shapeIn[0] *= self.compressionFactor
        return tuple(shapeIn)

    def __str__(self):
        return super(Compression, self).__str__()+' (x'+str(self.compressionFactor)+')'
        

#-------------------------------------------------------------------------------


class CompressionAverage(Compression):
    """
    Compress the input signal by averaging blocks of specified size.
    Author: P. Chanial
    """
    def __init__(self, compressionFactor, description=None):
        Compression.__init__(self, tmmf.compression_average_direct, tmmf.compression_average_transpose, compressionFactor, description)
        

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

    def __init__(self, array, pointingTime, ra, dec, pa, chop, fineTime, npixelsPerSample, observingMode, fineSamplingFactor, compressionFactor, badPixelMask, keepBadDetectors):

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
        self.keepBadDetectors = keepBadDetectors

        self.ndetectors = self.badPixelMask.size
        if not self.keepBadDetectors:
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
    - npixelsPerSample   : number of sky pixels which intersect a PACS detector
    - observingMode      : 'prime', 'parallel' or 'transparent'
    - fineSamplingFactor : number of frames to be processed during the 0.025s period.
    - compressionFactor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - badPixelMask       : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - keepBadDetectors : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, filename, first=0, last=None, header=None, resolution=None, npixelsPerSample=9, fineSamplingFactor=1, badPixelMask=None, keepBadDetectors=False):

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

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fineTime, npixelsPerSample, observingMode, fineSamplingFactor, compressionFactor, badPixelMask, keepBadDetectors)

        if resolution is None:
            resolution = 3.
        if header is None:
            header = self._str2fitsheader(tmmf.pacs_map_header(array, time, ra, dec, pa, chop, fineTime, self.transparentMode, self.badPixelMask.astype('int8'), keepBadDetectors, resolution))

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
        signal, mask = tmmf.pacs_timeline(self.filename, self.array, self.first+1, self.last, self.ndetectors, self.transparentMode, self.badPixelMask, self.keepBadDetectors)
        mask = numpy.array(mask, dtype=bool)
        if self.keepBadDetectors:
           for idetector, badPixel in enumerate(self.badPixelMask.flat):
               if badPixel:
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
    - npixelsPerSample   : number of sky pixels which intersect a PACS detector
    - observingMode      : 'prime', 'parallel' or 'transparent'
    - fineSamplingFactor : number of frames to be processed during the 1/40Hz period.
    - compressionFactor  : number of frames which are collapsed during the compression stage
    - array              : 'blue', 'green' or 'red' array name
    - ndetectors         : number of detectors involved in the processing
    - badPixelMask       : (nx,ny) mask of uint8 values (0 or 1). 1 means dead pixel.
    - keepBadDetectors : if set to True, force the processing to include all pixels
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, inputmap, time=None, ra=None, pa=None, chop=None, array='blue', npixelsPerSample=9, observingMode='prime', fineSamplingFactor=1, compressionFactor=None, badPixelMask=None, keepBadDetectors=False):

        if not isinstance(inputmap, Map):
            raise TypeError('The input is not a Map.')

        if inputmap.header is None:
            raise TypeError('The input map header is not known.')

        fineTime = numpy.arange(time[0], floor((time[-1] - time[0]) / pacs.sampling) * fineSamplingFactor, _Pacs.sampling/fineSamplingFactor)

        _Pacs.__init__(self, array, time, ra, dec, pa, chop, fineTime, npixelsPerSample, observingMode, fineSamplingFactor, compressionFactor, badPixelMask, keepBadDetectors)
 
        self.inputmap
        self.header = inputmap.header

    def get_tod(self, model, noise=None):
        """
        Returns simulated signal and mask (=None) timelines.
        """
        signal = model.direct(self.inputmap)
        if self.keepBadDetectors:
            mask = numpy.zeros(signal.shape, dtype=bool, order='fortran')
            for idetector, badPixel in enumerate(self.badPixelMask.flat):
                if badPixel:
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

    def __init__(self, data, header=None, **options):
        self.header = header
        self.set_fill_value(0)

    def copy(self, order='any'):
        newheader = None if self.header is None else self.header.copy()
        return FitsMaskedArray(super(FitsMaskedArray,self).copy(order), header=newheader)

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

     def __init__(self, bad_detectors=None, **options):
         self.bad_detectors = bad_detectors




#-------------------------------------------------------------------------------
#
# Miscellaneous routines
#
#-------------------------------------------------------------------------------



def hcssPhotProject(pacs):
    """
    Returns a map, as calculated by HCSS's PhotProject
    """
    if pacs.fineSamplingFactor != 1:
        raise ValueError('Fine sampling factor should be 1 for hcssPhotProject.') # or add decimation
    tod = pacs.get_tod()
    model = Masking(tod.mask) * PacsProjectionSharpEdges(pacs)
    mymap = model.transpose(tod).copy('any')
    tod[:] = 1.
    weights = model.transpose(tod)
    mymap /= weights
    return mymap

 
