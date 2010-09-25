import glob
import matplotlib
import numpy
import os
import pyfits
import re
import tamasisfortran as tmf

from acquisitionmodels import AcquisitionModel, CompressionAverage, Masking, Projection, ValidationError
from config import tamasis_dir
from datatypes import *
from mappers import mapper_naive
from mpi4py import MPI
from processing import deglitch_l2mad, filter_median
from unit import Quantity
from utils import MaskPolicy

__all__ = [ 'PacsObservation', 'pacs_plot_scan', 'pacs_preprocess', 'pacs_status' ]


class PacsObservation(object):
    """
    Class which encapsulates handy information about the PACS instrument and the processed
    observation.
    It is assumed that the detectors have sharp edges, so that the integration is simply equal
    to the sum of the sky pixels values weighted by their intersections with the detector surface.
    It contains the following attributes:
    - filename           : name of the file name, including the array colour, but excluding
                           the type. Example: '1342184520_blue'.
    - npixels_per_sample : number of sky pixels which intersect a PACS detector
    - keep_bad_detectors : if set to True, force the processing to include all pixels
    - detector_mask  : (nrows,ncolumns) mask of uint8 values (0 or 1). 1 means dead pixel.
    - ij(2,ndetectors)   : the row and column number (starting from 0) of the detectors
    Author: P. Chanial
    """
    def __init__(self, filename, fine_sampling_factor=1, detector_mask=None, reject_bad_line=False, frame_policy_inscan='keep', frame_policy_turnaround='keep', frame_policy_other='remove', frame_policy_invalid='mask'):

        if type(filename) == str:
            filename = (filename,)
        filename_, nfilenames = _files2tmf(filename)

        channel, transparent_mode, status = tmf.pacs_info_channel(filename_, nfilenames)
        if status != 0: raise RuntimeError()
        channel = channel.strip()

        nrows, ncolumns = (16,32) if channel == 'Red' else (32,64)

        # get the detector mask, before distributing to the processors
        if detector_mask is not None:
            if detector_mask.shape != (nrows, ncolumns):
                raise ValueError('Invalid shape of the input detector mask: ' + str(detector_mask.shape) + ' for the ' + channel + ' channel.')
            detector_mask = numpy.array(detector_mask, dtype='int8', copy=False)

        else:
            detector_mask, status = tmf.pacs_info_bad_detector_mask(tamasis_dir, channel, transparent_mode, reject_bad_line, nrows, ncolumns)
            if status != 0: raise RuntimeError()
            detector_mask = numpy.ascontiguousarray(detector_mask)

        # get the observations and detector mask for the current processor
        slice_observation, slice_detector = _split_observation(nfilenames, int(numpy.sum(detector_mask == 0)))
        filename = filename[slice_observation]
        igood = numpy.where(detector_mask.flat == 0)[0]
        detector_mask = numpy.ones(detector_mask.shape, dtype='int8')
        detector_mask.flat[igood[slice_detector]] = 0
        filename_, nfilenames = _files2tmf(filename)
        ndetectors = int(numpy.sum(detector_mask == 0))

        # frame policy
        frame_policy = MaskPolicy('inscan,turnaround,other,invalid'.split(','), (frame_policy_inscan, frame_policy_turnaround, frame_policy_other, frame_policy_invalid), 'Frame Policy')

        # retrieve information from the observations
        nthreads = tmf.info_nthreads()
        print 'Info: ' + MPI.Get_processor_name() + ', ' + str(nthreads) + ' core' + ('s' if nthreads > 1 else '') + ' (' + str(ndetectors) + ')'
        compression_factor, nsamples, unit, responsivity, detector_area, dflat, oflat, status = tmf.pacs_info(tamasis_dir, filename_, nfilenames, transparent_mode, fine_sampling_factor, numpy.array(frame_policy), numpy.asfortranarray(detector_mask))
        if status != 0: raise RuntimeError()

        self.filename = filename
        self.channel = channel
        self.nrows = nrows
        self.ncolumns = ncolumns
        self.default_npixels_per_sample = 11 if channel == 'Red' else 6
        self.default_resolution = 6.4 if channel == 'Red' else 3.2
        self.nobservations = nfilenames
        self.nsamples = tuple(nsamples)
        self.nfinesamples = tuple(nsamples * compression_factor * fine_sampling_factor)
        self.ndetectors = ndetectors
        self.detector_mask = detector_mask
        self.reject_bad_line = reject_bad_line
        self.frame_policy = frame_policy
        self.fine_sampling_factor = fine_sampling_factor
        self.transparent_mode = False if transparent_mode == 0 else True
        self.compression_factor = compression_factor
        self.unit = unit.strip()
        if self.unit.find('/') == -1:
            self.unit += ' / detector'
        self.responsivity = Quantity(responsivity, 'V/Jy')
        self.detector_area = Map(detector_area, unit='arcsec^2/detector')
        self.flatfield = {
            'total'   : Map(dflat*oflat),
            'detector': Map(numpy.ascontiguousarray(dflat)),
            'optical' : Map(numpy.ascontiguousarray(oflat))
            }

    def get_map_header(self, resolution=None, oversampling=True):
        if MPI.COMM_WORLD.Get_size() > 1:
            raise ErrorNotImplemented('The common map header should be specified if more than one job is running.')
        if resolution is None:
            resolution = self.default_resolution
        filename_, nfilenames = _files2tmf(self.filename)
        header, status = tmf.pacs_map_header(tamasis_dir, filename_, nfilenames, oversampling, self.fine_sampling_factor, numpy.array(self.frame_policy), numpy.asfortranarray(self.detector_mask), resolution)
        if status != 0: raise RuntimeError()
        header = _str2fitsheader(header)
        return header
   
    def get_tod(self, unit=None, flatfielding=True, subtraction_mean=True):
        """
        Returns the signal and mask timelines.
        """
        filename_, nfilenames = _files2tmf(self.filename)
        signal, mask, status = tmf.pacs_timeline(tamasis_dir, filename_, self.nobservations, numpy.sum(self.nsamples), self.ndetectors, numpy.array(self.frame_policy), numpy.asfortranarray(self.detector_mask), flatfielding, subtraction_mean)
        if status != 0: raise RuntimeError()
       
        tod = Tod(signal.T, mask.T, nsamples=self.nsamples, unit=self.unit)

        # the flux calibration has been done by using HCSS photproject and assuming that the central detectors had a size of 3.2x3.2
        # squared arcseconds. To be consistent with detector sharp edge model, we need to adjust the Tod.
        if tod.unit == 'Jy / detector' or tod.unit == 'V / detector':
            i = self.nrows    / 2 - 1
            j = self.ncolumns / 2 - 1
            tod *= numpy.mean(self.detector_area[i:i+2,j:j+2]) / Quantity(3.2**2, 'arcsec^2/detector')

        if unit is None:
            return tod

        newunit = Quantity(1., unit)
        newunit_si = newunit.SI._unit
        if 'sr' in newunit_si and newunit_si['sr'] == -1:
            area = self.detector_area[self.detector_mask == 0].reshape((self.ndetectors,1))
            tod /= area
           
        if 'V' in tod._unit and tod._unit['V'] == 1 and 'V' not in newunit_si:
            tod /= self.responsivity

        tod.unit = newunit._unit
        return tod

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, method=None, oversampling=True):
        if method is None:
            method = 'sharp edges'
        method = method.lower().replace('_', ' ').strip()
        if method not in ('nearest neighbour', 'sharp edges'):
            raise ValueError("Invalid method '" + method + "'. Valids methods are 'nearest neighbour' or 'sharp edges'")

        nsamples = self.nfinesamples if oversampling else self.nsamples
        if npixels_per_sample is None:
            npixels_per_sample = self.default_npixels_per_sample if method != 'nearest neighbour' else 1
        if header is None:
            if MPI.COMM_WORLD.Get_size() > 1:
                raise ValueError('In parallel mode, the map header must be speficied.')
            header = self.get_map_header(resolution, oversampling)
        elif isinstance(header, str):
            header = _str2fitsheader(header)

        filename_, nfilenames = _files2tmf(self.filename)
        sizeofpmatrix = npixels_per_sample * numpy.sum(nsamples) * self.ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.zeros(sizeofpmatrix, dtype=numpy.int64)
        
        status = tmf.pacs_pointing_matrix_filename(tamasis_dir, filename_, self.nobservations, method, oversampling, self.fine_sampling_factor, npixels_per_sample, numpy.sum(nsamples), self.ndetectors, numpy.array(self.frame_policy), numpy.asfortranarray(self.detector_mask), str(header).replace('\n', ''), pmatrix)
        if status != 0: raise RuntimeError()

        return pmatrix, header, self.ndetectors, nsamples, npixels_per_sample

    def get_filter_uncorrelated(self):
        """
        Reads an inverse noise time-time correlation matrix from a calibration file, in PACS-DP format.
        """
        ncorrelations, status = tmf.pacs_read_filter_calibration_ncorrelations(tamasis_dir, self.channel)
        if status != 0: raise RuntimeError()

        data, status = tmf.pacs_read_filter_calibration(tamasis_dir, self.channel, ncorrelations, self.ndetectors, numpy.asfortranarray(self.detector_mask))
        if status != 0: raise RuntimeError()

        return data.T

   
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

    def direct(self, signal, reusein=False, reuseout=False):
        signal, shapeout = self.validate_input_direct(Tod, signal, reusein)
        output = self.validate_output_direct(Tod, shapeout, reuseout)
        output.nsamples = tuple(numpy.divide(signal.nsamples, self.fine_sampling_factor))
        tmf.pacs_multiplexing_direct(signal.T, output.T, self.fine_sampling_factor, self.ij)
        return output

    def transpose(self, signal, reusein=False, reuseout=False):
        signal, shapein = self.validate_input_transpose(Tod, signal, reusein)
        output = self.validate_output_transpose(Tod, shapein, reuseout)
        output.nsamples = tuple(numpy.multiply(signal.nsamples, self.fine_sampling_factor))
        tmf.pacs_multiplexing_transpose(signal.T, output.T, self.fine_sampling_factor, self.ij)
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        if shapein[1] % self.fine_sampling_factor != 0:
            raise ValidationError('The input timeline size ('+str(shapein[1])+') is not an integer times the fine sampling factor ('+str(self.fine_sampling_factor)+').')
        shapeout = list(shapein)
        shapeout[1] = shapeout[1] / self.fine_sampling_factor
        return tuple(shapeout)

    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return
        super(PacsMultiplexing, self).validate_shapeout(shapeout)
        shapein = list(shapeout)
        shapein[1] = shapein[1] * self.fine_sampling_factor
        return tuple(shapein)


#-------------------------------------------------------------------------------


def pacs_plot_scan(patterns):
    if type(patterns) not in (tuple, list):
        patterns = (patterns,)

    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    for file in files:
        try:
            status = pacs_status(file)
        except IOError as error:
            raise IOError("Cannot extract status from file '"+file+"': "+str(error))
        matplotlib.pyplot.plot(status['ra'], status['dec'])


#-------------------------------------------------------------------------------


class pacs_status(object):
    def __init__(self, filename):
        hdu = pyfits.open(filename)[2]
        while True:
            try:
                self.status = hdu.data
                break
            except IndexError, errmsg:
                pass

    def __getitem__(self, key):
        if key == 'ra':
            return Quantity(self.status.field('RaArray'), 'deg')
        if key == 'dec':
            return Quantity(self.status.field('DecArray'), 'deg')
        if key == 'pa':
            return Quantity(self.status.field('PaArray'), 'deg')
        if key == 'time':
            return Quantity(self.status.field('FineTime')*1.e-6, 's')
        if key == 'velocity':
            ra = self['ra']
            dec = self['dec']
            time = self['time']
            dra  = numpy.diff(ra)
            ddec = numpy.diff(dec)
            dtime = numpy.diff(time)
            vel = numpy.sqrt((dra*numpy.cos(dec[0:-1].inunit('rad')))**2 + ddec**2) / dtime
            vel.unit = 'arcsec/s'
            u = vel._unit
            vel = numpy.append(vel, vel[-1])
            # BUG: append eats the unit...
            vel._unit = u
            return vel
        return self.status.field(key)
    
    def __str__(self):
        names = ['ra', 'dec', 'pa', 'time']
        for n in self.status.names:
            names.append(n)
        names.insert(0, 'ra')
        names.insert(1, 'dec')
        names.insert(2, 'pa')
        names.insert(3, 'time')
        names.insert(3, 'velocity')
        return 'PACS status: ' + str(names)
        

#-------------------------------------------------------------------------------


def pacs_preprocess(obs, projection_method='sharp edges', oversampling=True, npixels_per_sample=None, deglitching_hf_length=20, compression_factor=1, nsigma=5., hf_length=30000, tod=None):
    """
    deglitch, filter and potentially compress if the observation is in transparent mode
    """
    if tod is None:
        projection = Projection(obs, method='sharp edges', oversampling=False, npixels_per_sample=npixels_per_sample)
        tod = obs.get_tod()
        tod.mask[:] = 0
        tod_filtered = filter_median(tod, deglitching_hf_length)
        tod.mask = deglitch_l2mad(tod_filtered, projection, nsigma=nsigma)
        tod = filter_median(tod, hf_length)
        masking = Masking(tod.mask)
        tod = masking(tod)
    else:
        projection = None
        masking = Masking(tod.mask)
    
    # get the proper projector
    if projection is None or projection_method != 'sharp edges' or oversampling and numpy.any(obs.compression_factor * obs.fine_sampling_factor > 1):
        projection = Projection(obs, method=projection_method, oversampling=oversampling, npixels_per_sample=npixels_per_sample)

    # bail out if not in transparent mode
    if not obs.transparent_mode or compression_factor == 1:
        model = CompressionAverage(obs.compression_factor) * projection
        map_mask = model.T(Tod(tod.mask, nsamples=tod.nsamples))
        model = masking * model
        return tod, model, mapper_naive(tod, model), map_mask

    # compress the transparent observation
    compression = CompressionAverage(compression_factor)
    todc = compression(tod)
    mask = compression(tod.mask)
    mask[mask != 0] = 1
    todc.mask = numpy.array(mask, dtype='int8')
    maskingc = Masking(todc.mask)

    model = compression * projection
    map_mask = model.T(tod.mask)
    model = masking * model
    print 'XXX PACS PREPROCESS CHECK compression(tod.mask) and model.T(mask)'

    return todc, model, mapper_naive(todc, model), map_mask


#-------------------------------------------------------------------------------


def _files2tmf(filename):
    nfilenames = len(filename)
    length = max(len(f) for f in filename)
    filename_ = ''
    for f in filename:
        filename_ += f + (length-len(f))*' '
    return filename_, nfilenames


#-------------------------------------------------------------------------------


def _split_observation(nobservations, ndetectors):
    nnodes  = MPI.COMM_WORLD.Get_size()
    nthreads = tmf.info_nthreads()

    # number of observations. They should approximatively be of the same length
    nx = nobservations

    # number of detectors, grouped by the number of cpu cores
    ny = int(numpy.ceil(float(ndetectors) / nthreads))

    # we start with the miminum blocksize and increase it until we find a configuration that covers all the observations
    blocksize = int(numpy.ceil(float(nx * ny) / nnodes))
    while True:
        # by looping over x first, we favor larger number of detectors and fewer number of observations per processor, to minimise inter-processor communication in case of correlations between detectors
        for xblocksize in range(1, blocksize+1):
            if float(blocksize) / xblocksize != blocksize // xblocksize:
                continue
            yblocksize = int(blocksize // xblocksize)
            nx_block = int(numpy.ceil(float(nx) / xblocksize))
            ny_block = int(numpy.ceil(float(ny) / yblocksize))
            if nx_block * ny_block <= nnodes:
                break
        if nx_block * ny_block <= nnodes:
            break
        blocksize += 1

    rank = MPI.COMM_WORLD.Get_rank()

    ix = rank // ny_block
    iy = rank %  ny_block

    # check that the processor has something to do
    if ix >= nx_block:
        iobservation = slice(0,0)
        idetector = slice(0,0)
    else:
        iobservation = slice(ix * xblocksize, (ix+1) * xblocksize)
        idetector    = slice(iy * yblocksize * nthreads, (iy+1) * yblocksize * nthreads)

    return iobservation, idetector
        

#-------------------------------------------------------------------------------


def _str2fitsheader(string):
    """
    Convert a string into a pyfits.Header object
    All cards are extracted from the input string until the END keyword is reached.
    """
    header = pyfits.Header()
    cards = header.ascardlist()
    iline = 0
    while (iline*80 < len(string)):
        line = string[iline*80:(iline+1)*80]
        if line[0:3] == 'END': break
        cards.append(pyfits.Card().fromstring(line))
        iline += 1
    return header


