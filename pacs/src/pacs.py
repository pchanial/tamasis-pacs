import glob
import kapteyn
import numpy
import os
import pyfits
import re
import tempfile

from matplotlib import pyplot
from mpi4py import MPI
from tamasis import numpyutils
from tamasis.core import *
from tamasis.config import __version__ as tamasis_version
from tamasis.observations import *

__all__ = [ 'PacsObservation', 'PacsSimulation', 'pacs_plot_scan', 'pacs_preprocess' ]

DEFAULT_RESOLUTION = {'blue':3.2, 'green':3.2, 'red':6.4}
DEFAULT_NPIXELS_PER_SAMPLE = {'blue':6, 'green':6, 'red':11}

class _Pacs(Observation):

    def get_pointing_matrix(self, header, resolution, npixels_per_sample, method=None, oversampling=True):
        if method is None:
            method = 'sharp'
        method = method.lower()
        if method not in ('nearest', 'sharp'):
            raise ValueError("Invalid method '" + method + "'. Valids methods are 'nearest' or 'sharp'")

        nsamples = self.get_nfinesamples() if oversampling else self.get_nsamples()
        if npixels_per_sample is None:
            npixels_per_sample = DEFAULT_NPIXELS_PER_SAMPLE[self.instrument.band] if method != 'nearest' else 1
        if header is None:
            if MPI.COMM_WORLD.Get_size() > 1:
                raise ValueError('With MPI, the map header must be specified.')
            header = self.get_map_header(resolution, oversampling)
        elif isinstance(header, str):
            header = _str2fitsheader(header)

        ndetectors = self.get_ndetectors()
        nvalids = int(numpy.sum(nsamples))
        sizeofpmatrix = npixels_per_sample * nvalids * ndetectors
        print 'Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the pointing matrix.'
        pmatrix = numpy.empty(sizeofpmatrix, dtype=numpy.int64)
        
        status = tmf.pacs_pointing_matrix(self.instrument.band,
                                          nvalids,
                                          numpy.ascontiguousarray(self.slice.nsamples_all, dtype='int32'),
                                          numpy.ascontiguousarray(self.slice.compression_factor, dtype='int32'),
                                          self.instrument.fine_sampling_factor,
                                          oversampling,
                                          numpy.ascontiguousarray(self.pointing.time),
                                          numpy.ascontiguousarray(self.pointing.ra),
                                          numpy.ascontiguousarray(self.pointing.dec),
                                          numpy.ascontiguousarray(self.pointing.pa),
                                          numpy.ascontiguousarray(self.pointing.chop),
                                          numpy.ascontiguousarray(self.pointing.masked, dtype='int8'),
                                          numpy.ascontiguousarray(self.pointing.removed,dtype='int8'),
                                          method,
                                          numpy.asfortranarray(self.instrument.detector_mask),
                                          self.get_ndetectors(),
                                          self.instrument.detector_center.base.base.swapaxes(0,1).copy().T,
                                          self.instrument.detector_corner.base.base.swapaxes(0,1).copy().T,
                                          numpy.asfortranarray(self.instrument.detector_area),
                                          self.instrument.distortion_yz.base.base.T,
                                          npixels_per_sample,
                                          str(header).replace('\n',''),
                                          pmatrix)
        if status != 0: raise RuntimeError()

        return pmatrix, header, ndetectors, nsamples, npixels_per_sample

    def get_filter_uncorrelated(self):
        """
        Read an inverse noise time-time correlation matrix from a calibration file, in PACS-DP format.
        """
        ncorrelations, status = tmf.pacs_read_filter_calibration_ncorrelations(tamasis_dir, self.instrument.band)
        if status != 0: raise RuntimeError()

        data, status = tmf.pacs_read_filter_calibration(tamasis_dir, self.instrument.band, ncorrelations, self.get_ndetectors(), numpy.asfortranarray(self.instrument.detector_mask))
        if status != 0: raise RuntimeError()

        return data.T

    @property
    def status(self):
        if self._status is not None:
            return self._status
        
        status = []
        for filename in self.slice.filename:
            hdu = pyfits.open(filename)['STATUS']
            while True:
                try:
                    s = hdu.data
                    break
                except IndexError, errmsg:
                    pass
            
            status.append(s.base)
        self._status = numpy.concatenate(status).view(numpy.recarray)
        return self._status

    def __str__(self):
        nthreads = tmf.info_nthreads()
        ndetectors = self.get_ndetectors()
        sp = len('Info: ')*' '
        unit = 'unknown' if self.slice[0].unit == '' else self.slice[0].unit
        if MPI.COMM_WORLD.Get_size() > 1:
            mpistr = 'Process '+str(MPI.COMM_WORLD.Get_rank()+1) + '/' + str(MPI.COMM_WORLD.Get_size()) + \
                ' on node ' + MPI.Get_processor_name() + ', '
        else:
            mpistr = ''
        result = '\nInfo: ' + mpistr + str(nthreads) + ' core' + ('s' if nthreads > 1 else '') + ' handling ' + str(ndetectors) + ' detector' + ('s' if ndetectors > 1 else '') + '\n'
        result += sp + self.instrument.band.capitalize() + ' band, unit is ' + unit + '\n'
        dest = 0
        for islice, slice in enumerate(self.slice):

            p = self.pointing[dest:dest+slice.nsamples_all]

            def _str_policy(array):
                if array.size != p.size:
                    raise ValueError('Size problem.')
                
                nkept    = numpy.sum(numpy.logical_and(array, numpy.logical_and(~p.masked, ~p.removed)))
                nmasked  = numpy.sum(numpy.logical_and(array, numpy.logical_and(p.masked, ~p.removed)))
                nremoved = numpy.sum(numpy.logical_and(array, p.removed))

                if nkept+nmasked+nremoved != numpy.sum(array):
                    raise ValueError('This case should not happen.')
                if nkept+nmasked+nremoved == 0:
                    return 'None'
                result = []
                if nkept != 0:
                    result.append(str(nkept)+ ' kept')
                if nmasked != 0:
                    result.append(str(nmasked)+ ' masked')
                if nremoved != 0:
                    result.append(str(nremoved)+ ' removed')
                return ', '.join(result)

            result += sp + 'Observation'
            if self.slice.size > 1:
                result += ' #' + str(islice+1)
            result += ' ' + slice.filename
            first = numpy.argmin(p.removed) + 1
            last = p.size - numpy.argmin(p.removed[::-1]) + 1
            result += '[' + str(first) + ':' + str(last) + ']:\n'
            result += sp + '      Compression: x' + str(slice.compression_factor) + '\n'
            result += sp + '      In-scan:     ' + _str_policy(p.info == Pointing.INSCAN) + '\n'
            result += sp + '      Turnaround:  ' + _str_policy(p.info == Pointing.TURNAROUND) + '\n'
            result += sp + '      Other:       ' + _str_policy(p.info == Pointing.OTHER)
            if islice + 1 < self.slice.size:
                result += '\n'

        return result

    def _get_detector_mask(self, band, detector_mask, transparent_mode, reject_bad_line):
        shape = (16,32) if band == 'red' else (32,64)
        if type(detector_mask) is str:
            if detector_mask == 'calibration':
                detector_mask, status = tmf.pacs_info_detector_mask(tamasis_dir, band, shape[0], shape[1])
                if status != 0: raise RuntimeError()
                detector_mask = numpy.ascontiguousarray(detector_mask)
            else:
                raise ValueError('Invalid specification for the detector_mask.')
        elif detector_mask is None:
            detector_mask = numpy.zeros(shape, dtype='uint8')
        else:
            if detector_mask.shape != shape:
                raise ValueError('Invalid shape of the input detector mask: ' + str(detector_mask.shape) + ' for the ' + band + ' band.')
            detector_mask = numpy.array(detector_mask, dtype='uint8', copy=False)

        # mask non-transmitting detectors in transparent mode
        if transparent_mode:
            if band == 'red':
                detector_mask[0:8,0:8] = 1
                detector_mask[0:8,16:] = 1
                detector_mask[8:,:]    = 1
            else:
                detector_mask[0:16,0:16] = 1
                detector_mask[0:16,32:]  = 1
                detector_mask[16:,:]     = 1

        # mask erratic line
        if reject_bad_line and band != 'red':
            detector_mask[11,16:32] = 1

        return detector_mask


#-------------------------------------------------------------------------------


class PacsObservation(_Pacs):
    """
    Class which encapsulates handy information about the PACS instrument and 
    the observations to be processed.
    """
    def __init__(self, filename, fine_sampling_factor=1, detector_mask='calibration', reject_bad_line=False, policy_inscan='keep', policy_turnaround='keep', policy_other='remove', policy_invalid='mask'):
        """
        Parameters
        ----------
        filename: string or array of string
              This argument indicates the filenames of the observations
              to be processed. The files must be FITS files, as saved by
              the HCSS method FitsArchive.
        fine_sampling_factor: integer
              Set this value to a power of two, which will increase the
              acquisition model's sampling frequency beyond 40Hz
        detector_mask: None, 'calibration' or int8 array
              If None, no detector will be filtered out. if 'calibration',
              the detector mask will be read from a calibration file.
              Otherwise, it must be of the same shape as the camera:
              (16,32) for the red band and (32,64) for the others.
              Use 0 to keep a detector and 1 to filter it out.
        reject_bad_line: boolean
              If True, the erratic line [11,16:32] (starting from 0) in
              the blue channel will be filtered out
        policy_inscan: 'keep', 'mask' or 'remove'
              This sets the policy for the in-scan frames
        policy_turnaround: 'keep', 'mask' or 'remove'
              This sets the policy for the turnaround frames
        policy_other: 'keep', 'mask' or 'remove'
              This sets the policy for the other frames
        policy_invalid: 'keep', 'mask' or 'remove'
              This sets the policy for the invalid frames

        Returns
        -------
        The returned object contains the following attributes:
        - instrument: information about the instrument
        - pointing: information about the pointings, as a recarray.
        Pointings for which the policy is 'remove' are still available.
        - slice: information about the observations as a recarray
        - status: the HCSS Frames' Status ArrayDataset, as a recarray
        Like the pointing attribute, it also contains the pointings
        for which the policy is 'remove'
        - policy: frame policy

        """
        if type(filename) == str:
            filename = (filename,)
        filename_, nfilenames = _files2tmf(filename)

        band, transparent_mode, nsamples_all, status = tmf.pacs_info_init(filename_, nfilenames)
        if status != 0: raise RuntimeError()
        band = band.strip()

        # get the detector mask, before distributing the detectors to the processors
        detector_mask = self._get_detector_mask(band, detector_mask, transparent_mode, reject_bad_line)

        # get the observations and detector mask for the current processor
        slice_observation, slice_detector = _split_observation(nfilenames, int(numpy.sum(detector_mask == 0)))
        filename = filename[slice_observation]
        nsamples_all = nsamples_all[slice_observation]
        nfilenames = len(filename)
        igood = numpy.where(detector_mask.flat == 0)[0]
        detector_mask = numpy.ones(detector_mask.shape, dtype='uint8')
        detector_mask.flat[igood[slice_detector]] = 0
        filename_, nfilenames = _files2tmf(filename)
        ndetectors = int(numpy.sum(detector_mask == 0))
        ftype = get_default_dtype_float()

        # frame policy
        policy = MaskPolicy('inscan,turnaround,other,invalid', (policy_inscan, policy_turnaround, policy_other, policy_invalid), 'Frame Policy')

        # retrieve information from the observation
        mode, compression_factor, unit, ra, dec, cam_angle, scan_angle, scan_length, scan_speed, scan_step, scan_nlegs, frame_time, frame_ra, frame_dec, frame_pa, frame_chop, frame_info, frame_masked, frame_removed, status = tmf.pacs_info_observation(filename_, nfilenames, numpy.array(policy, dtype='int32'), numpy.sum(nsamples_all))
        flen_value = len(unit) // nfilenames
        mode = [mode[i*flen_value:(i+1)*flen_value].strip() for i in range(nfilenames)]
        unit = [unit[i*flen_value:(i+1)*flen_value].strip() for i in range(nfilenames)]
        unit = map(lambda x: x if x.find('/') != -1 else x + ' / detector', unit)

        # Store instrument information
        detector_center, detector_corner, detector_area, distortion_yz, oflat, dflat, responsivity, status = tmf.pacs_info_instrument(tamasis_dir, band, numpy.asfortranarray(detector_mask))
        if status != 0: raise RuntimeError()
        self.instrument = Instrument('PACS/' + band.capitalize(), detector_mask)
        self.instrument.band = band
        self.instrument.reject_bad_line = reject_bad_line
        self.instrument.fine_sampling_factor = fine_sampling_factor
        self.instrument.detector_center = detector_center.T.swapaxes(0,1).copy().view(dtype=[('u',ftype),('v',ftype)]).view(numpy.recarray)
        self.instrument.detector_corner = detector_corner.T.swapaxes(0,1).copy().view(dtype=[('u',ftype),('v',ftype)]).view(numpy.recarray)
        self.instrument.detector_area = Map(numpy.ascontiguousarray(detector_area), unit='arcsec^2/detector', origin='upper')
        self.instrument.distortion_yz = distortion_yz.T.view(dtype=[('y',ftype), ('z',ftype)]).view(numpy.recarray)
        self.instrument.flatfield = FlatField(oflat, dflat)
        self.instrument.responsivity = Quantity(responsivity, 'V/Jy')

        # Store slice information
        self.slice = numpy.recarray(nfilenames, dtype=[('filename', 'S256'), ('nsamples_all', int), ('mode', 'S32'), ('compression_factor', int), ('unit', 'S32'), ('ra', float), ('dec', float), ('cam_angle', float), ('scan_angle', float), ('scan_length', float), ('scan_nlegs', int), ('scan_step', float), ('scan_speed', float)])

        regex = re.compile(r'(.*?)(\[[0-9]*:?[0-9]*\])? *$')
        for ifile, file in enumerate(filename):
            match = regex.match(file)
            self.slice[ifile].filename = match.group(1)
        self.slice.nsamples_all = nsamples_all
        self.slice.nfinesamples = nsamples_all * compression_factor * fine_sampling_factor
        self.slice.mode = mode
        self.slice.compression_factor = compression_factor
        self.slice.unit = unit
        self.slice.ra = ra
        self.slice.dec = dec
        self.slice.cam_angle = cam_angle
        self.slice.scan_angle = scan_angle
        self.slice.scan_length = scan_length
        self.slice.scan_nlegs = scan_nlegs
        self.slice.scan_step = scan_step
        self.slice.scan_speed = scan_speed
        
        # Store pointing information
        self.pointing = PacsPointing(frame_time, frame_ra, frame_dec, frame_pa, frame_chop, frame_info, frame_masked, frame_removed, nsamples=self.slice.nsamples_all)

        # Store frame policy
        self.policy = policy

        # Status
        self._status = None

        print self

    def get_tod(self, unit=None, flatfielding=True, subtraction_mean=True, raw_data=False):
        """
        Returns the signal and mask timelines.
        """

        if raw_data:
            flatfielding = False
            subtraction_mean = False

        signal, mask, status = tmf.pacs_tod(self.instrument.band,
                                            _files2tmf(self.slice.filename)[0],
                                            numpy.asarray(self.slice.nsamples_all, dtype='int32'),
                                            numpy.asarray(self.slice.compression_factor, dtype='int32'),
                                            self.instrument.fine_sampling_factor,
                                            numpy.ascontiguousarray(self.pointing.time),
                                            numpy.ascontiguousarray(self.pointing.ra),
                                            numpy.ascontiguousarray(self.pointing.dec),
                                            numpy.ascontiguousarray(self.pointing.pa),
                                            numpy.ascontiguousarray(self.pointing.chop),
                                            numpy.ascontiguousarray(self.pointing.masked, dtype='int8'),
                                            numpy.ascontiguousarray(self.pointing.removed,dtype='int8'),
                                            numpy.asfortranarray(self.instrument.detector_mask),
                                            numpy.asfortranarray(self.instrument.flatfield.detector),
                                            flatfielding,
                                            subtraction_mean,
                                            int(numpy.sum(self.get_nsamples())),
                                            self.get_ndetectors())
        if status != 0: raise RuntimeError()
       
        tod = Tod(signal.T, mask.T, nsamples=self.get_nsamples(), unit=self.slice[0].unit, copy=False)

        # the flux calibration has been done by using HCSS photproject and assuming that the central detectors had a size of 3.2x3.2
        # squared arcseconds. To be consistent with detector sharp edge model, we need to adjust the Tod.
        if not raw_data and tod.unit in ('Jy / detector', 'V / detector'):
            i = self.instrument.shape[0] / 2 - 1
            j = self.instrument.shape[1] / 2 - 1
            nominal_area = (6.4 if self.instrument.band == 'red' else 3.2)**2
            tod *= numpy.mean(self.instrument.detector_area[i:i+2,j:j+2]) / Quantity(nominal_area, 'arcsec^2/detector')

        if unit is None:
            return tod

        newunit = Quantity(1., unit)
        newunit_si = newunit.SI._unit
        if 'sr' in newunit_si and newunit_si['sr'] == -1:
            area = self.instrument.detector_area[self.instrument.detector_mask == 0].reshape((self.get_ndetectors(),1))
            tod /= area
           
        if 'V' in tod._unit and tod._unit['V'] == 1 and 'V' not in newunit_si:
            tod /= self.instrument.responsivity

        tod.unit = newunit._unit
        return tod

    def get_map_header(self, resolution=None, oversampling=True):
        if MPI.COMM_WORLD.Get_size() > 1:
            raise NotImplementedError('The common map header should be specified if more than one job is running.')
        if resolution is None:
            resolution = DEFAULT_RESOLUTION[self.instrument.band]

        header, status = tmf.pacs_map_header(self.instrument.band,
                                             numpy.ascontiguousarray(self.slice.nsamples_all, dtype='int32'),
                                             numpy.ascontiguousarray(self.slice.compression_factor, dtype='int32'),
                                             self.instrument.fine_sampling_factor,
                                             oversampling,
                                             numpy.ascontiguousarray(self.pointing.time),
                                             numpy.ascontiguousarray(self.pointing.ra),
                                             numpy.ascontiguousarray(self.pointing.dec),
                                             numpy.ascontiguousarray(self.pointing.pa),
                                             numpy.ascontiguousarray(self.pointing.chop),
                                             numpy.ascontiguousarray(self.pointing.masked, dtype='int8'),
                                             numpy.ascontiguousarray(self.pointing.removed,dtype='int8'),
                                             numpy.asfortranarray(self.instrument.detector_mask),
                                             self.instrument.detector_corner.base.base.swapaxes(0,1).copy().T,
                                             self.instrument.distortion_yz.base.base.T,
                                             resolution)
        if status != 0: raise RuntimeError()
        header = _str2fitsheader(header)
        return header
   

#-------------------------------------------------------------------------------


class PacsPointing(Pointing):
    """
    This class subclasses tamasis.Pointing, by adding the chopper information.
    """
    def __new__(cls, time, ra, dec, pa, chop, info=None, masked=None, removed=None, nsamples=None):
        if nsamples is None:
            nsamples = (time.size,)
        else:
            if numpyutils._my_isscalar(nsamples):
                nsamples = (nsamples,)
            else:
                nsamples = tuple(nsamples)
        nsamples_tot = numpy.sum(nsamples)
        if numpy.any(numpy.array([ra.size,dec.size,pa.size,info.size,masked.size,removed.size]) != nsamples_tot):
            print numpy.array([ra.size,dec.size,pa.size,info.size,masked.size,removed.size]), nsamples_tot, nsamples
            raise ValueError('The pointing inputs do not have the same size.')

        if info is None:
            info = Pointing.INSCAN

        if masked is None:
            masked = False

        if removed is None:
            removed = False

        result = numpy.recarray(nsamples_tot, dtype=PacsPointing._dtype)
        result.time = time
        result.ra   = ra
        result.dec  = dec
        result.pa   = pa
        result.chop = chop
        result.info = info
        result.masked = masked
        result.removed = removed
        result = result.view(cls)
        result.nsamples = nsamples
        return result

    _dtype = [('time', get_default_dtype_float()), ('ra', get_default_dtype_float()), ('dec', get_default_dtype_float()),
              ('pa', get_default_dtype_float()), ('chop', get_default_dtype_float()), ('info', numpy.int64), 
              ('masked', numpy.bool_), ('removed', numpy.bool_)]


#-------------------------------------------------------------------------------


class PacsSimulation(_Pacs):
    """
    This class creates a simulated PACS observation.
    """
    def __init__(self, band, center, mode='prime', cam_angle=0., scan_angle=0., scan_length=30., scan_nlegs=3, scan_step=20., scan_speed=10., fine_sampling_factor=1, detector_mask='calibration', reject_bad_line=False, policy_inscan='keep', policy_turnaround='keep', policy_other='remove', policy_invalid='mask'):
        band = band.lower()
        if band not in ('blue', 'green', 'red'):
            raise ValueError("Band is not 'blue', 'green', nor 'red'.")

        mode = mode.lower()
        if mode not in ('prime', 'parallel', 'transparent'):
            raise ValueError("Observing mode is not 'prime', 'parallel', nor 'transparent'.")
        
        compression_factor = 8 if mode == 'parallel' and band != 'red' else 1 if mode == 'transparent' else 4
        detector_mask = self._get_detector_mask(band, detector_mask, mode == 'transparent', reject_bad_line)
        ftype = get_default_dtype_float()

        # Store pointing information
        self.pointing = _generate_pointing(center[0], center[1], cam_angle, scan_angle, scan_length, scan_nlegs, scan_step,
                                           scan_speed, compression_factor, dtype=PacsPointing._dtype)
        self.pointing.chop = 0.
        self.pointing.removed = policy_inscan == 'remove' and self.pointing.info == Pointing.INSCAN or \
                                policy_turnaround == 'remove' and self.pointing.info == Pointing.TURNAROUND

        # Store instrument information
        detector_center, detector_corner, detector_area, distortion_yz, oflat, dflat, responsivity, status = tmf.pacs_info_instrument(tamasis_dir, band, numpy.asfortranarray(detector_mask))
        if status != 0: raise RuntimeError()
        self.instrument = Instrument('PACS/'+band.capitalize(),detector_mask)
        self.instrument.band = band
        self.instrument.reject_bad_line = reject_bad_line
        self.instrument.fine_sampling_factor = fine_sampling_factor
        self.instrument.detector_center = detector_center.T.swapaxes(0,1).copy().view(dtype=[('u',ftype),('v',ftype)]).view(numpy.recarray)
        self.instrument.detector_corner = detector_corner.T.swapaxes(0,1).copy().view(dtype=[('u',ftype),('v',ftype)]).view(numpy.recarray)
        self.instrument.detector_area = Map(numpy.ascontiguousarray(detector_area), unit='arcsec^2/detector', origin='upper')
        self.instrument.distortion_yz = distortion_yz.T.view(dtype=[('y',ftype), ('z',ftype)]).view(numpy.recarray)
        self.instrument.flatfield = FlatField(oflat, dflat)
        self.instrument.responsivity = Quantity(responsivity, 'V/Jy')

        self.slice = numpy.recarray(1, dtype=[('filename', 'S256'), ('nsamples_all', int), ('mode', 'S32'), ('compression_factor', int), ('unit', 'S32'), ('ra', float), ('dec', float), ('cam_angle', float), ('scan_angle', float), ('scan_length', float), ('scan_nlegs', int), ('scan_step', float), ('scan_speed', float)])
        self.slice.filename = ''
        self.slice.nsamples_all = self.pointing.size
        self.slice.mode = mode
        self.slice.compression_factor = compression_factor
        self.slice.unit = ''
        self.slice.ra = center[0]
        self.slice.dec = center[1]
        self.slice.cam_angle = cam_angle
        self.slice.scan_angle = scan_angle
        self.slice.scan_length = scan_length
        self.slice.scan_nlegs = scan_nlegs
        self.slice.scan_step = scan_step
        self.slice.scan_speed = scan_speed
        self.slice.ninscans = numpy.sum(self.pointing.info == Pointing.INSCAN)
        self.slice.nturnarounds = numpy.sum(self.pointing.info == Pointing.TURNAROUND)
        self.slice.nothers = numpy.sum(self.pointing.info == Pointing.OTHER)
        self.slice.ninvalids = 0
        
        # Store policy
        self.policy = MaskPolicy('inscan,turnaround,other,invalid', 'keep,keep,remove,mask', 'Frame Policy')

        self._status = _write_status(self)

        print self

    def save(self, filename, tod):
        
        if numpy.rank(tod) != 3:
            tod = self.unpack(tod)

        nsamples = numpy.sum(self.slice.nsamples_all)
        if nsamples != tod.shape[-1]:
            raise ValueError("The input Tod has a number of samples'" + str(tod.shape[-1]) + "' incompatible with that of this observation '" + str(nsamples) + "'.")

        _write_status(self, filename)
        if tod.header is None:
            header = create_fitsheader(tod, extname='Signal')
        else:
            header = tod.header.copy()
            header.update('EXTNAME', 'Signal')

        if tod.unit != '':
            header.update('QTTY____', tod.unit)

        pyfits.append(filename, tod, header)
        if tod.mask is not None:
            mask = numpy.abs(tod.mask).view('uint8')
            header = create_fitsheader(mask, extname='Mask')
            pyfits.append(filename, mask, header)
        
   
#-------------------------------------------------------------------------------


class PacsMultiplexing(AcquisitionModel):
    """
    Performs the multiplexing of the PACS subarrays. The subarray columns are read one after the
    other, in a 0.025s cycle (40Hz).
    Author: P. Chanial
    """
    def __init__(self, obs, description=None):
        AcquisitionModel.__init__(self, description)
        self.fine_sampling_factor = obs.instrument.fine_sampling_factor
        self.ij = obs.instrument.ij

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


def pacs_plot_scan(patterns, title=None, new_figure=True):
    if type(patterns) not in (tuple, list):
        patterns = (patterns,)

    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    for ifile, file in enumerate(files):
        try:
            hdu = pyfits.open(file)['STATUS']
        except Exception as error:
            print "Warning: Cannot extract status from file '"+file+"': "+str(error)
            continue

        while True:
            try:
                status = hdu.data
                break
            except IndexError:
                pass

        if ifile == 0:
            image = plot_scan(status.RaArray, status.DecArray, title=title, new_figure=new_figure)
        else:
            x, y = image.topixel(status.RaArray, status.DecArray)
            p = pyplot.plot(x, y, linewidth=2)
            pyplot.plot(x[0], y[0], 'o', color = p[0]._color)


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
    if projection is None or projection_method != 'sharp edges' or oversampling and numpy.any(obs.slice.compression_factor * obs.instrument.fine_sampling_factor > 1):
        projection = Projection(obs, method=projection_method, oversampling=oversampling, npixels_per_sample=npixels_per_sample)

    # bail out if not in transparent mode
    if not obs.slice[0].mode == 'transparent' or compression_factor == 1:
        model = CompressionAverage(obs.slice.compression_factor) * projection
        map_mask = model.T(Tod(tod.mask, nsamples=tod.nsamples))
        model = masking * model
        return tod, model, mapper_naive(tod, model), map_mask

    # compress the transparent observation
    compression = CompressionAverage(compression_factor)
    todc = compression(tod)
    mask = compression(tod.mask)
    mask[mask != 0] = 1
    todc.mask = numpy.array(mask, dtype='uint8')
    maskingc = Masking(todc.mask)

    model = compression * projection
    map_mask = model.T(tod.mask)
    model = masking * model
    print 'XXX PACS PREPROCESS CHECK compression(tod.mask) and model.T(mask)'

    return todc, model, mapper_naive(todc, model), map_mask


#-------------------------------------------------------------------------------


def _change_coord(ra0, dec0, pa0, lon, lat):
    """
    Transforms the longitude and latitude coordinates expressed in the
    native coordinate system attached to a map into right ascension and
    declination.
    Author: R. Gastaud
    """
    from numpy import arcsin, arctan2, cos, deg2rad, mod, sin, rad2deg

    # Arguments in radian
    lambd = deg2rad(lon)
    beta  = deg2rad(lat) 
    alpha = deg2rad(ra0)
    delta = deg2rad(dec0)
    eps   = deg2rad(pa0)

    # Cartesian coordinates from longitude and latitude
    z4 = sin(beta)
    y4 = sin(lambd)*cos(beta)
    x4 = cos(lambd)*cos(beta)

    # rotation about the x4 axe of -eps 
    x3 =  x4
    y3 =  y4*cos(eps) + z4*sin(eps)
    z3 = -y4*sin(eps) + z4*cos(eps)
    
    # rotation about the axis Oy2, angle delta
    x2 = x3*cos(delta) - z3*sin(delta)
    y2 = y3
    z2 = x3*sin(delta) + z3*cos(delta)

    # rotation about the axis Oz1, angle alpha
    x1 = x2*cos(alpha) - y2*sin(alpha)
    y1 = x2*sin(alpha) + y2*cos(alpha)
    z1 = z2

    # Compute angles from cartesian coordinates
    # it is the only place where we can get nan with arcsinus
    dec = rad2deg(arcsin(numpy.clip(z1, -1., 1.)))
    ra  = mod(rad2deg(arctan2(y1, x1)), 360.)

    return ra, dec


#-------------------------------------------------------------------------------


def _files2tmf(filename):
    nfilenames = len(filename)
    length = max(len(f) for f in filename)
    filename_ = ''
    for f in filename:
        filename_ += f + (length-len(f))*' '
    return filename_, nfilenames


#-------------------------------------------------------------------------------


def _generate_pointing(ra0, dec0, pa0, scan_angle, scan_length=30., scan_nlegs=3, scan_step=20., scan_speed=10., 
                       compression_factor=4, dtype=None):
    """
    compute the pointing timeline of the instrument reference point
    from the description of a scan map
    Authors: R. Gastaud
    """
    
    # some info
    gamma = 4.
    sampling_pacs = 0.024996

    scan_angle -= 90.

    scan_length = float(scan_length)
    if scan_length <= 0:
        raise ValueError('Input scan_length must be strictly positive.')

    scan_nlegs = int(scan_nlegs)
    if scan_nlegs <= 0:
        raise ValueError('Input scan_nlegs must be strictly positive.')
    
    scan_step = float(scan_step)
    if scan_step <= 0: 
        raise ValueError('Input scan_step must be strictly positive.')

    scan_speed = float(scan_speed)
    if scan_speed <= 0:
        raise ValueError('Input scan_speed must be strictly positive.')

    compression_factor = int(compression_factor)
    if compression_factor not in (1, 4, 8):
        raise ValueError("Input compression_factor must be 1, 4 or 8.")
    sampling_frequency = 1. / sampling_pacs / compression_factor
    sampling_period    = 1. / sampling_frequency

    # compute the different times and the total number of points
    # acceleration time at the beginning of a leg, and deceleration time
    # at the end
    extra_time1 = scan_speed / gamma
    # corresponding length 
    extralength = 0.5 * gamma * extra_time1 * extra_time1
    # Time needed to go from a scan line to the next 
    extra_time2 = numpy.sqrt(scan_step / gamma)
    # Time needed to turn around (not used)
    turnaround_time = 2 * (extra_time1 + extra_time2)

    # Time needed to go along the scanline at constant speed
    line_time = scan_length / scan_speed
    # Total time for a scanline
    full_line_time=extra_time1 + line_time + extra_time1 + extra_time2 + extra_time2 
    # Total duration of the observation
    total_time = full_line_time * scan_nlegs - 2 * extra_time2

    # Number of samples
    nsamples = int(numpy.ceil(total_time * sampling_frequency))

    # initialization
    time          = numpy.zeros(nsamples)
    latitude      = numpy.zeros(nsamples)
    longitude     = numpy.zeros(nsamples)
    infos         = numpy.zeros(nsamples, dtype=int)
    line_counters = numpy.zeros(nsamples, dtype=int)

    # Start of computations, alpha and delta are the longitude and
    # latitide in arc seconds in the referential of the map.
    signe = 1
    delta = -extralength - scan_length/2.
    alpha = -scan_step * (scan_nlegs-1)/2.
    alpha0 = alpha
    line_counter = 0
    working_time = 0.

    for i in xrange(nsamples):
        info = 0
   
        # check if new line
        if working_time > full_line_time:
            working_time = working_time - full_line_time
            signe = -signe
            line_counter = line_counter + 1
            alpha = -scan_step * (scan_nlegs-1)/2. + line_counter * scan_step
            alpha0 = alpha
   
        # acceleration at the beginning of a scan line to go from 0 to the
        # scan_speed. 
        if working_time < extra_time1:
            delta = -signe*(extralength + scan_length/2) + signe * 0.5 * gamma * working_time * working_time
            info  = Pointing.TURNAROUND
   
        # constant speed
        if working_time >=  extra_time1 and working_time < extra_time1+line_time:
            delta = signe*(-scan_length/2+ (working_time-extra_time1)*scan_speed)
            info  = Pointing.INSCAN
 
        # Deceleration at then end of the scanline to stop
        if working_time >= extra_time1+line_time and working_time < extra_time1+line_time+extra_time1:
            dt = working_time - extra_time1 - line_time
            delta = signe * (scan_length/2 + scan_speed*dt - 0.5 * gamma*dt*dt)
            info  = Pointing.TURNAROUND
  
        # Acceleration to go toward the next scan line
        if working_time >= 2*extra_time1+line_time and working_time < 2*extra_time1+line_time+extra_time2:
            dt = working_time-2*extra_time1-line_time
            alpha = alpha0 + 0.5*gamma*dt*dt
            info  = Pointing.TURNAROUND
   
        # Deceleration to stop at the next scan line
        if working_time >= 2*extra_time1+line_time+extra_time2 and working_time < full_line_time:
            dt = working_time-2*extra_time1-line_time-extra_time2
            speed = gamma*extra_time2
            alpha = (alpha0+scan_step/2.) + speed*dt - 0.5*gamma*dt*dt
            info  = Pointing.TURNAROUND

        time[i] = i / sampling_frequency
        infos[i] = info
        latitude[i] = delta
        longitude[i] = alpha
        line_counters[i] = line_counter
        working_time = working_time + sampling_period

    # Convert the longitude and latitude *expressed in degrees) to ra and dec
    ra, dec = _change_coord(ra0, dec0, scan_angle, longitude/3600., latitude/3600.)

    return Pointing(time, ra, dec, pa0, infos, dtype=dtype)


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


#-------------------------------------------------------------------------------


def _write_status(obs, filename=None):

    s = obs.slice[0]

    if any(obs.slice.compression_factor != obs.slice[0].compression_factor):
        raise ValueError('Unable to save into a single file. The observations do not have the same compression factor.')
    compression_factor = s.compression_factor

    if any(obs.slice.mode != obs.slice[0].mode):
        raise ValueError('Unable to save into a single file. The observations do not have the same observing mode.')
    mode = s.mode
    band_status = {'blue':'BS', 'green':'BL', 'red':'R '}[obs.instrument.band]
    band_type   = {'blue':'BS', 'green':'BL', 'red':'RS'}[obs.instrument.band]

    cusmode = {'prime':'PacsPhoto', 'parallel':'SpirePacsParallel', 'transparent':'__PacsTranspScan', 'unknown':'__Calibration'}
    if mode == 'prime' or obs.instrument.band == 'red':
        comp_mode = 'Photometry Default Mode'
    elif mode == 'parallel':
        comp_mode = 'Photometry Double Compression Mode'
    elif mode == 'transparent':
        comp_mode = 'Photometry Lossless Compression Mode'
    else:
        comp_mode = ''

    if obs.slice[0].scan_speed == 10.:
        scan_speed = 'low'
    elif obs.slice[0].scan_speed == 20.:
        scan_speed = 'medium'
    elif obs.slice[0].scan_speed == 60.:
        scan_speed = 'high'
    else:
        scan_speed = str(obs.slice[0].scan_speed)
        
    p = obs.pointing
    if p.size != numpy.sum(obs.slice.nsamples_all):
        raise ValueError('The pointing and slice attribute are incompatible. This case should not happen.')

    # get status
    table = numpy.recarray(p.size, dtype=[('BBID', numpy.int64), ('FINETIME', numpy.int64), ('BAND', 'S2'), ('CHOPFPUANGLE', numpy.float64), ('RaArray', numpy.float64), ('DecArray', numpy.float64), ('PaArray', numpy.float64)])
    table.BAND     = band_status
    table.FINETIME = numpy.round(p.time*1000000.)
    table.RaArray  = p.ra
    table.DecArray = p.dec
    table.PaArray  = p.pa
    table.CHOPFPUANGLE = 0. if not hasattr(p, 'chop') else p.chop
    table.BBID = 0
    #XXX Numpy ticket #1645
    table.BBID[p.info == Pointing.INSCAN] = 0xcd2 << 16
    table.BBID[p.info == Pointing.TURNAROUND] = 0x4000 << 16

    if filename is None:
        return table

    fits = pyfits.HDUList()

    # Primary header
    cc = pyfits.createCard
    
    header = pyfits.Header([
            cc('simple', True), 
            cc('BITPIX', 32), 
            cc('NAXIS', 0), 
            cc('EXTEND', True, 'May contain datasets'), 
            cc('TYPE', 'HPPAVG'+band_type, 'Product Type Identification'), 
            cc('CREATOR', 'TAMASIS v' + tamasis_version, 'Generator of this file'), 
            cc('INSTRUME', 'PACS', 'Instrument attached to this file'), 
            cc('TELESCOP', 'Herschel Space Observatory', 'Name of telescope'),
            cc('OBS_MODE', 'Scan map', 'Observation mode name'),
            cc('RA', s.ra),
            cc('DEC', s.dec),
            cc('EQUINOX', 2000., 'Equinox of celestial coordinate system'),
            cc('RADESYS', 'ICRS', 'Coordinate reference frame for the RA and DEC'),
            cc('CUSMODE', cusmode[mode], 'CUS observation mode'),
            cc('META_0', obs.instrument.shape[0]),
            cc('META_1', obs.instrument.shape[1]), 
            cc('META_2', obs.instrument.band.title()+' Photometer'), 
            cc('META_3', ('Floating Average  : ' + str(compression_factor)) if compression_factor > 1 else 'None'), 
            cc('META_4', comp_mode),
            cc('META_5', s.scan_angle),
            cc('META_6', s.scan_length),
            cc('META_7', s.scan_nlegs),
            cc('META_8', scan_speed),
            cc('HIERARCH key.META_0', 'detRow'), 
            cc('HIERARCH key.META_1', 'detCol'), 
            cc('HIERARCH key.META_2', 'camName'), 
            cc('HIERARCH key.META_3', 'algorithm'), 
            cc('HIERARCH key.META_4', 'compMode'),
            cc('HIERARCH key.META_5', 'mapScanAngle'),
            cc('HIERARCH key.META_6', 'mapScanLegLength'),
            cc('HIERARCH key.META_7', 'mapScanNumLegs'),
            cc('HIERARCH key.META_8', 'mapScanSpeed'),
            ])

    hdu = pyfits.PrimaryHDU(None, header)
    fits.append(hdu)
    
    status = pyfits.BinTableHDU(table, None, 'STATUS')
    fits.append(status)
    fits.writeto(filename, clobber=True)

    return table
