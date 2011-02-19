import glob
import kapteyn
import numpy as np
import os
import pyfits
import re
import tempfile

from . import var
from matplotlib import pyplot as P
from mpi4py import MPI
from tamasis.core import *
from tamasis.mpiutils import split_observation
from tamasis.observations import Observation, Instrument, FlatField, create_scan
from tamasis.stringutils import strenum, strplural
from tamasis.wcsutils import combine_fitsheader

__all__ = [ 'PacsObservation',
            'PacsSimulation',
            'pacs_create_scan',
            'pacs_get_psf',
            'pacs_plot_scan',
            'pacs_preprocess' ]

CALIBFILE_DTC = var.path + '/pacs/PCalPhotometer_PhotTimeConstant_FM_v2.fits'

class PacsBase(Observation):

    ACCELERATION = Quantity(4., '"/s^2')
    DEFAULT_RESOLUTION = {'blue':3.2, 'green':3.2, 'red':6.4}
    PSF_FWHM = {'blue':5.2, 'green':7.7, 'red':12.}
    POINTING_DTYPE = [('time', var.FLOAT_DTYPE), ('ra', var.FLOAT_DTYPE),
                      ('dec', var.FLOAT_DTYPE), ('pa', var.FLOAT_DTYPE),
                      ('chop', var.FLOAT_DTYPE), ('info', np.int64),
                      ('masked', np.bool8), ('removed', np.bool8)]
    SAMPLING_PERIOD = Quantity(0.024996, 's')


    def get_derived_units(self):
        volt = Quantity(1./self.instrument.responsivity, 'Jy')
        return (
            {
                'detector_reference': Quantity(1 / \
                    self.instrument.active_fraction, 'detector'),
                'detector' : self.get_detector_area(),
                'V' : volt,
            },{
                'V' : volt
            }
            )

    def get_detector_area(self):
        """
        Return the area of the non masked detectors.
        """
        return self.instrument.detector_area[~self.instrument.detector_mask]

    def get_detector_time_constant(self):
        """
        Return the time constants of the non masked detectors.
        """
        return self.instrument \
                   .detector_time_constant[~self.instrument.detector_mask]

    def get_filter_uncorrelated(self):
        """
        Read an inverse noise time-time correlation matrix from a calibration
        file, in PACS-DP format.
        """
        ncorrelations, status = tmf.pacs_read_filter_calibration_ncorrelations(
            self.instrument.band)
        if status != 0: raise RuntimeError()

        data, status = tmf.pacs_read_filter_calibration(
            self.instrument.band, ncorrelations, self.get_ndetectors(),
            np.asfortranarray(self.instrument.detector_mask, np.int8))
        if status != 0: raise RuntimeError()

        return data.T

    def get_map_header(self, resolution=None, oversampling=True):
        if var.mpi_comm.Get_size() > 1:
            raise NotImplementedError('The common map header should be '
                'specified if more than one job is running.')
        if resolution is None:
            resolution = self.DEFAULT_RESOLUTION[self.instrument.band]

        header, status = tmf.pacs_map_header(
            self.instrument.band,
            np.ascontiguousarray(self.slice.nsamples_all, np.int32),
            np.ascontiguousarray(self.slice.compression_factor, np.int32),
            np.ascontiguousarray(self.slice.delay),
            self.instrument.fine_sampling_factor,
            oversampling,
            np.ascontiguousarray(self.pointing.time),
            np.ascontiguousarray(self.pointing.ra),
            np.ascontiguousarray(self.pointing.dec),
            np.ascontiguousarray(self.pointing.pa),
            np.ascontiguousarray(self.pointing.chop),
            np.ascontiguousarray(self.pointing.masked, np.int8),
            np.ascontiguousarray(self.pointing.removed,np.int8),
            np.asfortranarray(self.instrument.detector_mask, np.int8),
            self.instrument.detector_corner.base.base.swapaxes(0,1).copy().T,
            self.instrument.distortion_yz.base.base.T,
            resolution)
        if status != 0: raise RuntimeError()
        headers = var.mpi_comm.allgather(_str2fitsheader(header))
        return combine_fitsheader(headers)
    get_map_header.__doc__ = Observation.get_map_header.__doc__
    
    def get_pointing_matrix(self, header, resolution, npixels_per_sample=0,
                            method=None, oversampling=True):
        if method is None:
            method = 'sharp'
        method = method.lower()
        choices = ('nearest', 'sharp')
        if method not in choices:
            raise ValueError("Invalid method '" + method + \
                "'. Expected values are " + strenum(choices, 'or') + '.')

        nsamples = self.get_nfinesamples() if oversampling else \
            self.get_nsamples()
        if header is None:
            if var.mpi_comm.Get_size() > 1:
                raise ValueError('With MPI, the map header must be specified.')
            header = self.get_map_header(resolution, oversampling)
        elif isinstance(header, str):
            header = _str2fitsheader(header)

        ndetectors = self.get_ndetectors()
        nvalids = int(np.sum(nsamples))
        if npixels_per_sample != 0:
            sizeofpmatrix = npixels_per_sample * nvalids * ndetectors
            print('Info: Allocating '+str(sizeofpmatrix/2.**17)+' MiB for the '
                  'pointing matrix.')
        else:
            sizeofpmatrix = 1
        pmatrix = np.empty(sizeofpmatrix, dtype=np.int64)
        
        new_npixels_per_sample, status = tmf.pacs_pointing_matrix(
            self.instrument.band,
            nvalids,
            np.ascontiguousarray(self.slice.nsamples_all, np.int32),
            np.ascontiguousarray(self.slice.compression_factor, np.int32),
            np.ascontiguousarray(self.slice.delay),
            self.instrument.fine_sampling_factor,
            oversampling,
            np.ascontiguousarray(self.pointing.time),
            np.ascontiguousarray(self.pointing.ra),
            np.ascontiguousarray(self.pointing.dec),
            np.ascontiguousarray(self.pointing.pa),
            np.ascontiguousarray(self.pointing.chop),
            np.ascontiguousarray(self.pointing.masked, np.int8),
            np.ascontiguousarray(self.pointing.removed,np.int8),
            method,
            np.asfortranarray(self.instrument.detector_mask, np.int8),
            self.get_ndetectors(),
            self.instrument.detector_center.base.base.swapaxes(0,1).copy().T,
            self.instrument.detector_corner.base.base.swapaxes(0,1).copy().T,
            np.asfortranarray(self.instrument.detector_area),
            self.instrument.distortion_yz.base.base.T,
            npixels_per_sample,
            str(header).replace('\n',''),
            pmatrix)
        if status != 0: raise RuntimeError()

        # the number of pixels per sample is now known, do the real computation
        if npixels_per_sample == 0:
            return self.get_pointing_matrix(header, resolution,
                new_npixels_per_sample, method, oversampling)

        return pmatrix, header, ndetectors, nsamples, npixels_per_sample, \
               ('/detector', '/pixel'), self.get_derived_units()
    get_pointing_matrix.__doc__ = Observation.get_pointing_matrix.__doc__

    def get_random(self, flatfielding=True, subtraction_mean=True):
        """
        Return noise data from a random slice of a real pointed observation.

        Currently, the required duration must be less than that of the 
        real observation (around 3 hours).
        """
        if any([slice.compression_factor not in [4,8] for slice in self.slice]):
            raise NotImplementedError('The compression factor must be 4 or 8.')

        path = os.path.join(var.path, 'pacs')
        files = { 'blue' : ('1342182424_blue_PreparedFrames.fits.gz',
                            10000,100000),
                  'green': ('1342182427_green_PreparedFrames.fits.gz',
                            10000, 100000),
                  'red': ('1342182427_red_PreparedFrames.fits.gz',
                          7400, 116400)}
        file, istart, iend = files[self.instrument.band]
        data = pyfits.open(os.path.join(path, file))['Signal'] \
                     .data[:,:,istart:iend]
        data /= self.instrument.active_fraction * self.instrument.responsivity
        ndetectors = self.get_ndetectors()
        nsamples = self.get_nsamples()
        nsamples_tot = np.sum(nsamples*self.slice.compression_factor/4)
        if nsamples_tot > data.shape[-1]:
            raise ValueError('There is not enough noise data for this observa' \
                             'tion.')

        irandom = np.random.random_integers(0,iend-istart-nsamples_tot)
        result = Tod(data[:,:,irandom:irandom+nsamples_tot],
                     unit='Jy/detector',
                     derived_units=self.get_derived_units()[0],
                     nsamples=nsamples)
        result = self.pack(result)
        if subtraction_mean:
            result.T[:] -= np.mean(result, axis=1)
        if flatfielding:
            result.T[:] /= self.instrument.flatfield \
                               .detector[~self.instrument.detector_mask]

        compression = CompressionAverage(self.slice.compression_factor/4)
        result = compression(result)
        return result

    def pack(self, tod):
        tod = Observation.pack(self, tod)
        if 'detector' in tod.derived_units:
            tod.derived_units['detector'] = self.get_detector_area()
        return tod
    pack.__doc__ = Observation.pack.__doc__
    
    def unpack(self, tod):
        tod = Observation.unpack(self, tod)
        if 'detector' in tod.derived_units:
            tod.derived_units['detector'] = self.instrument.detector_area
        return tod
    unpack.__doc__ = Observation.unpack.__doc__
    
    def save(self, filename, tod):
        
        if np.rank(tod) != 3:
            tod = self.unpack(tod)

        nsamples = np.sum(~self.pointing.removed)
        if nsamples != tod.shape[-1]:
            raise ValueError("The input Tod has a number of samples'" + \
                str(tod.shape[-1]) + "' incompatible with that of this observ" \
                "ation '" + str(nsamples) + "'.")

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
            _write_mask(self, tod.mask, filename)
    save.__doc__ = Observation.save.__doc__
        
    def __str__(self):

        # some helpers
        def same(a):
            return np.all(a == a[0])
        def ad_masks(slice, islice):
            if 'nmasks' not in slice.dtype.names:
                a = []
                d = []
            else:
                masks  = slice[islice].mask_name[0:slice[islice].nmasks]
                active = slice[islice].mask_activated[0:slice[islice].nmasks]
                a = [m for m, act in zip(masks, active) if act]
                d = [m for m, act in zip(masks, active) if not act]
            return (strplural(s + 'activated mask', len(m), False, \
                ': ').capitalize() + ', '.join(m) for s,m in (('',a), ('de',d)))

        nthreads = tmf.info_nthreads()
        ndetectors = self.get_ndetectors()
        sp = len('Info ')*' '
        unit = 'unknown' if self.slice[0].unit == '' else self.slice[0].unit
        if var.mpi_comm.Get_size() > 1:
            mpistr = 'Process '+str(var.mpi_comm.Get_rank()+1) + '/' + \
                     str(var.mpi_comm.Get_size()) + ', '
        else:
            mpistr = ''
        
        # print general information
        result = '\nInfo ' + MPI.Get_processor_name() + ': ' + mpistr + \
                 strplural('core', nthreads) + ' handling ' + \
                 strplural('detector', ndetectors) + '\n'
        info = [ self.instrument.band.capitalize() + ' band' ]

        # observing mode
        print_each_mode = not same(self.slice.mode)
        if not print_each_mode:
            info.append(self.slice[0].mode + ' mode')

        # compression factor
        print_each_compression = not same(self.slice.compression_factor)
        if not print_each_compression:
            info.append('compression factor x' + \
                        str(self.slice[0].compression_factor))

        # scan speed
        scan_speed = []
        scan_speed_reliable = []
        isample = 0
        for slice in self.slice:
            p = self.pointing[isample:isample+slice.nsamples_all]
            isample += slice.nsamples_all
            nsamples = np.sum(~p.removed)
            inscan = (p.info == self.pointing.INSCAN) * ~p.removed
            ninscans = np.sum(inscan)
            if ninscans == 0:
                inscan[:] = True
            scan_speed_reliable.append('' if ninscans >= 0.25 * nsamples and \
                                       nsamples > 10 else '~')
            scan_speed.append(np.round(np.median(p.velocity[inscan]),1))
        print_each_scan_speed = not same(scan_speed)
        if not print_each_scan_speed:
            info.append('scan speed ' + scan_speed_reliable[0] + \
                        str(scan_speed[0]))
        
        # unit
        if self.__class__.__name__ == 'PacsObservation':
            info.append('unit is ' + unit)

        result += sp + ', '.join(info) + '\n'

        # mask information
        if self.__class__.__name__ == 'PacsObservation':        
            homogeneous = 'nmasks' not in self.slice.dtype.names or \
                same(self.slice.nmasks) and all([same(a) for a in \
                self.slice.mask_name.T])
            if homogeneous:
                (a,d) = ad_masks(self.slice, 0)
                result += sp + a + '\n'
                result += sp + d + '\n'
            else:
                result += sp + 'The masks of the observations are heterogeneo' \
                          'us\n'

        # print slice-specific information
        isample = 0
        for islice, slice in enumerate(self.slice):

            p = self.pointing[isample:isample+slice.nsamples_all]
            isample += slice.nsamples_all

            def _str_policy(array):
                if array.size != p.size:
                    raise ValueError('Size problem.')
                
                nkept    = np.sum(array * ~p.masked * ~p.removed)
                nmasked  = np.sum(array * p.masked * ~p.removed)
                nremoved = np.sum(array * p.removed)

                if nkept+nmasked+nremoved != np.sum(array):
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
            first = np.argmin(p.removed) + 1
            last = p.size - np.argmin(p.removed[::-1])
            result += '[' + str(first) + ':' + str(last) + ']:\n'

            if self.__class__.__name__ == 'PacsObservation' and not homogeneous:
                (a,d) = ad_masks(self.slice, islice)
                result += sp + '      ' + a + '\n'
                result += sp + '      ' + d + '\n'

            if print_each_mode:
                result += sp + '     Mode:        ' + slice.mode + '\n'
            if print_each_compression:
                result += sp + '     Compression: x' + \
                    str(slice.compression_factor) + '\n'
            if print_each_scan_speed:
                result += sp + '     Scan speed:  ' + \
                    scan_speed_reliable[islice] + str(scan_speed[islice]) + '\n'
            result += sp + '     In-scan:     ' + \
                _str_policy(p.info == Pointing.INSCAN) + '\n'
            result += sp + '     Turnaround:  ' + \
                _str_policy(p.info == Pointing.TURNAROUND) + '\n'
            result += sp + '     Other:       ' + \
                _str_policy(p.info == Pointing.OTHER)
            if islice + 1 < self.slice.size:
                result += '\n'

        return result

    def _get_detector_mask(self, band, detector_mask, transparent_mode,
                           reject_bad_line):
        shape = (16,32) if band == 'red' else (32,64)
        if type(detector_mask) is str:
            if detector_mask == 'calibration':
                detector_mask, status = tmf.pacs_info_detector_mask(
                    band, shape[0], shape[1])
                if status != 0: raise RuntimeError()
                detector_mask = np.ascontiguousarray(detector_mask, np.bool8)
            else:
                raise ValueError('Invalid specification for the detector_mask.')
        elif detector_mask is None:
            detector_mask = np.zeros(shape, np.bool8)
        else:
            if detector_mask.shape != shape:
                raise ValueError('Invalid shape of the input detector mask: ' +\
                    str(detector_mask.shape) + ' for the ' + band + ' band.')
            detector_mask = np.array(detector_mask, np.bool8, copy=False)

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


class PacsObservation(PacsBase):
    """
    Class which encapsulates handy information about the PACS instrument and 
    the observations to be processed.
    """
    def __init__(self, filename, fine_sampling_factor=1,
                 detector_mask='calibration', reject_bad_line=False,
                 policy_inscan='keep', policy_turnaround='keep',
                 policy_other='remove', policy_invalid='mask',
                 active_fraction=0, delay=0.):
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
        detector_mask: None, 'calibration' or boolean array
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
        active_fraction: ratio of the geometric and reference detector area
              If set to 0, the value from the calibration file will be used
        delay: instrument clock lag wrt the spacecraft clock, in ms

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

        band, transparent_mode, nsamples_all, status = tmf.pacs_info_init(
            filename_, nfilenames)
        if status != 0: raise RuntimeError()
        band = band.strip()

        # get the detector mask, before distributing the detectors to the
        # processors
        detector_mask = self._get_detector_mask(band, detector_mask,
            transparent_mode, reject_bad_line)

        # get the observations and detector mask for the current processor
        slice_observation, slice_detector = split_observation(var.mpi_comm,
             int(np.sum(detector_mask == 0)), nfilenames)
        filename = filename[slice_observation]
        nsamples_all = nsamples_all[slice_observation]
        nfilenames = len(filename)
        igood = np.where(detector_mask.flat == 0)[0]
        detector_mask = np.ones(detector_mask.shape, np.bool8)
        detector_mask.flat[igood[slice_detector]] = 0
        filename_, nfilenames = _files2tmf(filename)

        # frame policy
        policy = MaskPolicy('inscan,turnaround,other,invalid', (policy_inscan,
            policy_turnaround, policy_other, policy_invalid), 'Frame Policy')

        # store observation information
        mode, compression_factor, unit, ra, dec, cam_angle, scan_angle, \
            scan_length, scan_speed, scan_step, scan_nlegs, frame_time, \
            frame_ra, frame_dec, frame_pa, frame_chop, frame_info, \
            frame_masked, frame_removed, nmasks, mask_name_flat, \
            mask_activated, status = tmf.pacs_info_observation(filename_,
                nfilenames, np.array(policy, np.int32),np.sum(nsamples_all))
        if status != 0: raise RuntimeError()

        flen_value = len(unit) // nfilenames
        mode = [mode[i*flen_value:(i+1)*flen_value].strip() \
                for i in range(nfilenames)]
        unit = [unit[i*flen_value:(i+1)*flen_value].strip() \
                for i in range(nfilenames)]

        # store instrument information
        detector_center, detector_corner, detector_area, distortion_yz, oflat, \
            dflat, responsivity, active_fraction, status = \
            tmf.pacs_info_instrument(band, np.asfortranarray(detector_mask,
                np.int8), float(active_fraction))
        if status != 0: raise RuntimeError()

        self.instrument = Instrument('PACS/' + band.capitalize(), detector_mask)
        self.instrument.active_fraction = active_fraction
        self.instrument.band = band
        self.instrument.reject_bad_line = reject_bad_line
        self.instrument.fine_sampling_factor = fine_sampling_factor
        self.instrument.detector_center = detector_center.T.swapaxes(0,1) \
            .copy().view(dtype=[('u',var.FLOAT_DTYPE),('v',var.FLOAT_DTYPE)]) \
            .view(np.recarray)
        self.instrument.detector_corner = detector_corner.T.swapaxes(0,1) \
            .copy().view(dtype=[('u',var.FLOAT_DTYPE),('v',var.FLOAT_DTYPE)]) \
            .view(np.recarray)
        self.instrument.detector_area = Map(np.ascontiguousarray(detector_area),
            unit='arcsec^2', origin='upper')
        self.instrument.distortion_yz = distortion_yz.T.view(dtype=[('y', \
            var.FLOAT_DTYPE), ('z',var.FLOAT_DTYPE)]).view(np.recarray)
        self.instrument.flatfield = FlatField(oflat, dflat)
        self.instrument.responsivity = Quantity(responsivity, 'V/Jy')
        self.instrument.detector_time_constant = Map(pyfits.open(CALIBFILE_DTC)\
            [{'red':1, 'green':2, 'blue':3}[band]].data, unit='ms',
            origin='upper')

        # store slice information
        nmasks_max = np.max(nmasks)
        if nmasks_max > 0:
            mask_len_max = np.max(
                [len(mask_name_flat[(i*32+j)*70:(i*32+j+1)*70].strip()) \
                 for j in range(nmasks[i]) for i in range(nfilenames)])
        else:
            mask_len_max = 1
        # 'mask_name' and 'mask_activated' would be scalar in slice:
        nmasks_max = max(nmasks_max, 2)

        mask_name = np.ndarray((nfilenames, nmasks_max), 'S'+str(mask_len_max))
        mask_activated = mask_activated.T
        for ifile in range(nfilenames):
            for imask in range(nmasks[ifile]):
                dest = (ifile*32+imask)*70
                mask_name[ifile,imask] = \
                    mask_name_flat[dest:dest+mask_len_max].lower()
            isort = np.argsort(mask_name[ifile,0:nmasks[ifile]])
            mask_name     [ifile,0:nmasks[ifile]] = mask_name     [ifile,isort]
            mask_activated[ifile,0:nmasks[ifile]] = mask_activated[ifile,isort]

        self.slice = np.recarray(nfilenames, dtype=[
                ('filename', 'S256'),
                ('nsamples_all', int),
                ('mode', 'S32'),
                ('compression_factor', int),
                ('delay', float),
                ('unit', 'S32'),
                ('ra', float),
                ('dec', float),
                ('cam_angle', float),
                ('scan_angle', float),
                ('scan_length', float),
                ('scan_nlegs', int),
                ('scan_step', float),
                ('scan_speed', float),
                ('nmasks', int),
                ('mask_name', 'S'+str(mask_len_max), nmasks_max),
                ('mask_activated', bool, nmasks_max)
                ])

        regex = re.compile(r'(.*?)(\[[0-9]*:?[0-9]*\])? *$')
        for ifile, file in enumerate(filename):
            match = regex.match(file)
            self.slice[ifile].filename = match.group(1)

        self.slice.nsamples_all = nsamples_all
        self.slice.nfinesamples = nsamples_all * compression_factor * \
                                  fine_sampling_factor
        self.slice.mode = mode
        self.slice.compression_factor = compression_factor
        self.slice.delay = delay
        self.slice.unit = unit
        self.slice.ra = ra
        self.slice.dec = dec
        self.slice.cam_angle = cam_angle
        self.slice.scan_angle = scan_angle
        self.slice.scan_length = scan_length
        self.slice.scan_nlegs = scan_nlegs
        self.slice.scan_step = scan_step
        self.slice.scan_speed = scan_speed
        self.slice.nmasks = nmasks
        self.slice.mask_name      = mask_name
        self.slice.mask_activated = mask_activated[:,0:nmasks_max]
        
        # store pointing information
        self.pointing = Pointing(frame_time, frame_ra, frame_dec, frame_pa,
            frame_info, frame_masked, frame_removed,
            nsamples=self.slice.nsamples_all, dtype=self.POINTING_DTYPE)
        self.pointing.chop = frame_chop

        # store frame policy
        self.policy = policy

        # status
        self._status = None

        print(self)

    def get_tod(self,
                unit='Jy/detector',
                flatfielding=True,
                subtraction_mean=True,
                raw=False,
                masks='activated'):
        """
        Returns the signal and mask timelines.

        By default, if no active mask is specified, the Master mask will
        be retrieved if it exists. Otherwise, the activated masks will
        be read and combined.
        """

        if raw:
            flatfielding = False
            subtraction_mean = False

        act_masks = set([m for slice in self.slice \
                         for i, m in enumerate(slice.mask_name) \
                         if m not in ('','master') and slice.mask_activated[i]])
        dea_masks = set([m for slice in self.slice \
                         for i, m in enumerate(slice.mask_name) \
                         if m not in ('','master') and \
                            not slice.mask_activated[i]])
        all_masks = set([m for slice in self.slice for m in slice.mask_name \
                         if m != ''])

        if isinstance(masks, str):
            masks = masks.split(',')

        masks = [m.strip().lower() for m in masks]
        sel_masks = set()

        if 'none' not in masks:
            for m in masks:
                if m == 'all':
                    sel_masks |= all_masks
                elif m == 'activated':
                    sel_masks |= act_masks
                elif m == 'deactivated':
                    sel_masks |= dea_masks
                elif m == '':
                    continue
                elif m not in all_masks:
                    print("Warning: mask '" + m + "' is not found.")
                else:
                    sel_masks.add(m)

        # use 'master' if all activated masks are selected
        if all(['master' in slice.mask_name for slice in self.slice]) and \
           act_masks <= sel_masks:
            sel_masks -= act_masks
            sel_masks.add('master')

        sel_masks = ','.join(sorted(sel_masks))

        signal, mask, status = tmf.pacs_tod(
            self.instrument.band,
            _files2tmf(self.slice.filename)[0],
            np.asarray(self.slice.nsamples_all, np.int32),
            np.asarray(self.slice.compression_factor, np.int32),
            np.ascontiguousarray(self.slice.delay),
            self.instrument.fine_sampling_factor,
            np.ascontiguousarray(self.pointing.time),
            np.ascontiguousarray(self.pointing.ra),
            np.ascontiguousarray(self.pointing.dec),
            np.ascontiguousarray(self.pointing.pa),
            np.ascontiguousarray(self.pointing.chop),
            np.ascontiguousarray(self.pointing.masked, np.int8),
            np.ascontiguousarray(self.pointing.removed,np.int8),
            np.asfortranarray(self.instrument.detector_mask, np.int8),
            np.asfortranarray(self.instrument.flatfield.detector),
            flatfielding,
            subtraction_mean,
            int(np.sum(self.get_nsamples())),
            self.get_ndetectors(),
            sel_masks)
        if status != 0: raise RuntimeError()
       
        tod = Tod(signal.T, 
                  mask.T,
                  nsamples=self.get_nsamples(),
                  unit=self.slice[0].unit,
                  derived_units=self.get_derived_units()[0],
                  copy=False)

        if not raw:
            tod.inunit(unit)
        return tod

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
                except IndexError as errmsg:
                    pass
            
            status.append(s.base)

        # check number of records
        if np.any([len(s) for s in status] != self.slice.nsamples_all):
            raise ValueError("The status has a number of records '" + \
                str([len(s) for s in status]) + "' incompatible with that of " \
                "the pointings '" + str(self.slice.nsamples_all) + "'.")

        # merge status if necessary
        if any([status[i].dtype != status[0].dtype \
                for i in range(1,len(status))]):
            newdtype = []
            for d in status[0].dtype.names:
                newdtype.append((d, max(status, key=lambda s: \
                                s[d].dtype.itemsize)[d].dtype))
            self._status = np.recarray(int(np.sum(self.slice.nsamples_all)),
                                       newdtype)
            dest = 0
            for s in status:
                for n in status[0].dtype.names:
                    self._status[n][dest:dest+len(s)] = s[n]
                dest += len(s)
        else:
            self._status = np.concatenate(status).view(np.recarray)
        return self._status


#-------------------------------------------------------------------------------


class PacsSimulation(PacsBase):
    """
    This class creates a simulated PACS observation.
    """
    def __init__(self, pointing, band, mode=None, fine_sampling_factor=1,
                 detector_mask='calibration', reject_bad_line=False,
                 policy_inscan='keep', policy_turnaround='keep',
                 policy_other='keep', policy_invalid='keep', active_fraction=0,
                 delay=0.):

        pointing_ = pointing.copy()
        if pointing.header is not None:
            pointing_.header = pointing.header.copy()
        pointing = pointing_

        band = band.lower()
        choices = ('blue', 'green', 'red')
        if band not in choices:
            raise ValueError("Invalid band '" + band + \
                "'. Expected values are " + strenum(choices, 'or') + '.')

        # get compression factor from input pointing
        delta = np.median(pointing.time[1:]-pointing.time[0:-1])
        compression_factor = int(np.round(delta / self.SAMPLING_PERIOD))
        if compression_factor <= 0:
            raise ValueError('Invalid time in pointing argument. Use PacsSimu' \
                             'lation.SAMPLING_PERIOD.')
        if np.abs(delta / compression_factor - self.SAMPLING_PERIOD) > 0.01 * \
           self.SAMPLING_PERIOD:
            print('Warning: the pointing time has unexpected sampling rate. ' \
                  'Assuming a compression factor of ' + str(compression_factor)\
                  + '.')

        # observing mode
        if mode is None:
            if compression_factor == 4:
                mode = 'prime' # might not be the case for the red band
            elif compression_factor == 8:
                mode = 'calibration' if band == 'red' else 'parallel'
            else:
                mode = 'calibration'
        else:
            mode = mode.lower()
            choices = ('prime', 'parallel', 'calibration', 'transparent')
            if mode not in choices:
                raise ValueError("Invalid observing mode '" + mode + \
                    "'. Expected values are " + strenum(choices, 'or') + '.')
            if mode == 'prime' and compression_factor != 4 or \
               mode == 'parallel' and compression_factor != (4 if band == 'red'\
                       else 8) or \
               mode == 'transparent' and compression_factor != 1:
                raise ValueError("The observing mode '" + mode + "' in the " + \
                    band + " band is incompatible with the compression factor" \
                    " '" + str(compression_factor) + "'.")

        # detector mask
        detector_mask = self._get_detector_mask(band, detector_mask,
            mode == 'transparent', reject_bad_line)

        # store pointing information
        if not hasattr(pointing, 'chop'):
            pointing.chop = np.zeros(pointing.size, var.FLOAT_DTYPE)
        pointing.masked  = \
            (policy_inscan     == 'mask')   * \
                (pointing.info == Pointing.INSCAN) + \
            (policy_turnaround == 'mask')   * \
                (pointing.info == Pointing.TURNAROUND) + \
            (policy_other      == 'mask')   * \
                (pointing.info == Pointing.OTHER)
        pointing.removed = \
            (policy_inscan     == 'remove') * \
                (pointing.info == Pointing.INSCAN) + \
            (policy_turnaround == 'remove') * \
                (pointing.info == Pointing.TURNAROUND) + \
            (policy_other      == 'remove') * \
                (pointing.info == Pointing.OTHER)
        self.pointing = pointing

        # store instrument information
        detector_center, detector_corner, detector_area, distortion_yz, oflat, \
            dflat, responsivity, active_fraction, status = \
            tmf.pacs_info_instrument(band, np.asfortranarray(detector_mask, \
                np.int8), float(active_fraction))
        if status != 0: raise RuntimeError()
        self.instrument = Instrument('PACS/'+band.capitalize(),detector_mask)
        self.instrument.active_fraction = active_fraction
        self.instrument.band = band
        self.instrument.reject_bad_line = reject_bad_line
        self.instrument.fine_sampling_factor = fine_sampling_factor
        self.instrument.detector_center = detector_center.T.swapaxes(0,1) \
            .copy().view(dtype=[('u',var.FLOAT_DTYPE),('v',var.FLOAT_DTYPE)]) \
            .view(np.recarray)
        self.instrument.detector_corner = detector_corner.T.swapaxes(0,1) \
            .copy().view(dtype=[('u',var.FLOAT_DTYPE),('v',var.FLOAT_DTYPE)]) \
            .view(np.recarray)
        self.instrument.detector_area = Map(np.ascontiguousarray(detector_area),
            unit='arcsec^2', origin='upper')
        self.instrument.distortion_yz = distortion_yz.T.view(dtype=[('y',
            var.FLOAT_DTYPE), ('z',var.FLOAT_DTYPE)]).view(np.recarray)
        self.instrument.flatfield = FlatField(oflat, dflat)
        self.instrument.responsivity = Quantity(responsivity, 'V/Jy')
        self.instrument.detector_time_constant = Map(pyfits.open(CALIBFILE_DTC)\
            [{'red':1, 'green':2, 'blue':3}[band]].data, unit='ms',
            origin='upper')

        self.slice = np.recarray(1, dtype=[
                ('filename', 'S256'),
                ('nsamples_all', int),
                ('mode', 'S32'),
                ('compression_factor', int),
                ('delay', float),
                ('unit', 'S32'),
                ('ra', float),
                ('dec', float),
                ('cam_angle', float),
                ('scan_angle', float),
                ('scan_length', float),
                ('scan_nlegs', int),
                ('scan_step', float),
                ('scan_speed', float)
                ])

        self.slice.filename = ''
        self.slice.nsamples_all = self.pointing.size
        self.slice.mode = mode
        self.slice.compression_factor = compression_factor
        self.slice.delay = delay
        self.slice.unit = ''
        for field in ('ra', 'dec', 'cam_angle', 'scan_angle', 'scan_nlegs',
                      'scan_length', 'scan_step', 'scan_speed'):
            self.slice[field] = pointing.header[field] if field in \
                pointing.header else 0
        self.slice.ninscans = np.sum(self.pointing.info == Pointing.INSCAN)
        self.slice.nturnarounds = np.sum(self.pointing.info==Pointing.TURNAROUND)
        self.slice.nothers = np.sum(self.pointing.info == Pointing.OTHER)
        self.slice.ninvalids = 0
        
        # store policy
        self.policy = MaskPolicy('inscan,turnaround,other,invalid',
                                 'keep,keep,remove,mask', 'Frame Policy')

        # status
        self._status = None

        print(self)

    @property
    def status(self):
        if self._status is not None:
            return self._status

        p = self.pointing
        if p.size != np.sum(self.slice.nsamples_all):
            raise ValueError('The pointing and slice attribute are incompatib' \
                             'le. This should not happen.')

        status = np.recarray(p.size, dtype=[('BBID', np.int64),
            ('FINETIME', np.int64), ('BAND', 'S2'), ('CHOPFPUANGLE', np.float64),
            ('RaArray', np.float64), ('DecArray', np.float64),
            ('PaArray', np.float64)])
        status.BAND     = {'blue':'BS', 'green':'BL', 'red':'R '} \
            [self.instrument.band]
        status.FINETIME = np.round(p.time*1000000.)
        status.RaArray  = p.ra
        status.DecArray = p.dec
        status.PaArray  = p.pa
        status.CHOPFPUANGLE = 0. if not hasattr(p, 'chop') else p.chop
        status.BBID = 0
        #XXX N ticket #1645
        status.BBID[p.info == Pointing.INSCAN] = 0xcd2 << 16
        status.BBID[p.info == Pointing.TURNAROUND] = 0x4000 << 16
        self._status = status

        return status
   
#-------------------------------------------------------------------------------


class PacsMultiplexing(AcquisitionModelLinear):
    """
    Performs the multiplexing of the PACS subarrays.
    The subarray columns are read one after the other, in a 0.025s cycle (40Hz).
    """
    def __init__(self, obs, shapein=None, description=None):
        AcquisitionModelLinear.__init__(self,
                                        cache=True,
                                        shapein=shapein,
                                        shapeout=shapeout,
                                        typein=Tod,
                                        description=description)
        self.fine_sampling_factor = obs.instrument.fine_sampling_factor
        self.ij = obs.instrument.ij

    def direct(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_direct(input, cachein, cacheout)
        tmf.pacs_multiplexing_direct(input.T, output.T,
                                     self.fine_sampling_factor, self.ij)
        return output

    def transpose(self, input, inplace, cachein, cacheout):
        input, output = self.validate_input_transpose(input, cachein, cacheout)
        tmf.pacs_multiplexing_transpose(input.T, output.T,
                                        self.fine_sampling_factor, self.ij)
        return output

    def validate_shapein(self, shapein):
        if shapein is None:
            return None
        if shapein[1] % self.fine_sampling_factor != 0:
            raise ValidationError('The input timeline size (' + \
                str(shapein[1]) + ') is not an integer times the fine sampling'\
                ' factor ('+str(self.fine_sampling_factor)+').')
        shapeout = list(shapein)
        shapeout[1] = shapeout[1] / self.fine_sampling_factor
        return tuple(shapeout)

    def validate_shapeout(self, shapeout):
        if shapeout is None:
            return None
        super(PacsMultiplexing, self).validate_shapeout(shapeout)
        shapein = list(shapeout)
        shapein[1] = shapein[1] * self.fine_sampling_factor
        return tuple(shapein)


#-------------------------------------------------------------------------------


def pacs_plot_scan(patterns, title=None, new_figure=True, **kw):

    if type(patterns) not in (tuple, list):
        patterns = (patterns,)

    files = []
    scans = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    for ifile, file in enumerate(files):
        try:
            hdu = pyfits.open(file)['STATUS']
        except Exception as error:
            print("Warning: Cannot extract status from file '" + file + \
                  "': " + str(error))
            continue

        while True:
            try:
                status = hdu.data
                break
            except IndexError:
                pass

        scans.append((status.RaArray, status.DecArray))
    
    plot_scan(scans, **kw)


#-------------------------------------------------------------------------------


def pacs_preprocess(obs, tod,
                    projection_method='sharp',
                    header=None,
                    oversampling=True,
                    npixels_per_sample=0,
                    deglitching_hf_length=20,
                    deglitching_nsigma=5.,
                    hf_length=30000,
                    transparent_mode_compression_factor=1):
    """
    deglitch, filter and potentially compress if the observation is in
    transparent mode
    """
    projection = Projection(obs,
                            method='sharp',
                            header=header,
                            oversampling=False,
                            npixels_per_sample=npixels_per_sample)
    tod_filtered = filter_median(tod, deglitching_hf_length)
    tod.mask = deglitch_l2mad(tod_filtered,
                              projection,
                              nsigma=deglitching_nsigma)
    tod = filter_median(tod, hf_length)
    masking = Masking(tod.mask)
    tod = masking(tod)
    
    # get the proper projector if necessary
    if projection is None or projection_method != 'sharp' or \
       oversampling and np.any(obs.slice.compression_factor * \
       obs.instrument.fine_sampling_factor > 1):
        projection = Projection(obs,
                                method=projection_method,
                                oversampling=oversampling,
                                header=header,
                                npixels_per_sample=npixels_per_sample)

    # bail out if not in transparent mode
    if all(obs.slice[0].compression_factor != 1) or \
       transparent_mode_compression_factor == 1:
        if oversampling:
            compression = CompressionAverage(obs.slice.compression_factor)
            model = compression * projection
        else:
            model = projection
        map_mask = model.T(Tod(tod.mask, nsamples=tod.nsamples))
        model = masking * model
        return tod, model, mapper_naive(tod, model), map_mask

    # compress the transparent observation
    compression = CompressionAverage(transparent_mode_compression_factor)
    todc = compression(tod)
    mask = compression(tod.mask)
    mask[mask != 0] = 1
    todc.mask = np.array(mask, dtype='uint8')
    maskingc = Masking(todc.mask)

    model = compression * projection
    map_mask = model.T(Tod(tod.mask, nsamples=tod.nsamples, copy=False))
    model = masking * model

    return todc, model, mapper_naive(todc, model), map_mask


#-------------------------------------------------------------------------------


def pacs_create_scan(ra0, dec0, cam_angle=0., scan_angle=0., scan_length=30.,
                     scan_nlegs=3, scan_step=148., scan_speed=20., 
                     compression_factor=4):
    if int(compression_factor) not in (1, 4, 8):
        raise ValueError("Input compression_factor must be 1, 4 or 8.")
    scan = create_scan(ra0, dec0,
                       PacsBase.ACCELERATION,
                       PacsBase.SAMPLING_PERIOD * compression_factor,
                       scan_angle=scan_angle,
                       scan_length=scan_length,
                       scan_nlegs=scan_nlegs,
                       scan_step=scan_step,
                       scan_speed=scan_speed,
                       dtype=PacsBase.POINTING_DTYPE)
    scan.header.update('HIERARCH cam_angle', cam_angle)
    scan.chop = 0.
    return scan


def pacs_get_psf(band, resolution, kind='calibration'):
    """
    Return gaussian, airy or calibration PSFs

    Calibration PSFs are rescaled to the required pixel resolution
    by a bilinear interpolation.

    Parameters
    ----------
    band : string
        PACS band ('blue', 'green' or 'red')
    resolution : float
        Pixel resolution of the PSF
    kind : string
        Values can be 'gaussian', 'airy' or 'calibration'
    """

    band = band.lower()
    choices = ('blue', 'green', 'red')
    if band not in choices:
        raise ValueError("Invalid band '" + band + "'. Expected values are " + \
                         strenum(choices, 'or') + '.')

    kind = kind.lower()
    choices = ('airy', 'gaussian', 'calibration')
    if kind not in choices:
        raise ValueError("Invalid PSF kind '" + kind + \
            "'. Expected values are " + strenum(choices, 'or') + '.')

    if kind in ('airy', 'gaussian'):
        func = { 'airy': airy_disk, 'gaussian' : gaussian }[kind]
        fwhm = PacsBase.PSF_FWHM[band]
        size = int(np.round(10 * fwhm / resolution)) // 2 * 2 + 1
        return func((size,size), fwhm=fwhm, resolution=resolution)

    raise NotImplementedError()


#-------------------------------------------------------------------------------


def _files2tmf(filename):
    nfilenames = len(filename)
    length = max(len(f) for f in filename)
    filename_ = ''
    for f in filename:
        filename_ += f + (length-len(f))*' '
    return filename_, nfilenames


#-------------------------------------------------------------------------------


def _str2fitsheader(string):
    """
    Convert a string into a pyfits.Header object

    All cards are extracted from the input string until the END keyword is
    reached.
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


def _write_mask(obs, mask, filename):

    header = create_fitsheader(np.ones(0), extname='Mask')
    header.update('DSETS___', 1)
    header.update('DS_0', 4)
    pyfits.append(filename, None, header)

    mask = mask.view()
    mask.shape = (np.product(obs.instrument.shape),-1)
    bitmask = tmf.pacs_bitmask(mask.T.view(np.int8))
    bitmask = bitmask.T
    bitmask.shape = (obs.instrument.shape[0], obs.instrument.shape[1],-1)
    header = create_fitsheader(bitmask, extname='Tamasis')
    header.update('INFO____', 'Tamasis mask')
    pyfits.append(filename, bitmask, header)


#-------------------------------------------------------------------------------


def _write_status(obs, filename):

    s = obs.slice[0]

    if any(obs.slice.compression_factor != obs.slice[0].compression_factor):
        raise ValueError('Unable to save into a single file. The observations '\
                         'do not have the same compression factor.')
    compression_factor = s.compression_factor

    if any(obs.slice.mode != obs.slice[0].mode):
        raise ValueError('Unable to save into a single file. The observations '\
                         'do not have the same observing mode.')
    mode = s.mode
    band_type   = {'blue':'BS', 'green':'BL', 'red':'RS'}[obs.instrument.band]

    cusmode = {'prime':'PacsPhoto', 'parallel':'SpirePacsParallel',
               'transparent':'__PacsTranspScan', 'unknown':'__Calibration'}
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
        
    if obs.pointing.size != np.sum(obs.slice.nsamples_all):
        raise ValueError('The pointing and slice attribute are incompatible. ' \
                         'This should not happen.')

    status = obs.status[~obs.pointing.removed]

    fits = pyfits.HDUList()

    # Primary header
    cc = pyfits.createCard
    
    header = pyfits.Header([
            cc('simple', True), 
            cc('BITPIX', 32), 
            cc('NAXIS', 0), 
            cc('EXTEND', True, 'May contain datasets'), 
            cc('TYPE', 'HPPAVG'+band_type, 'Product Type Identification'), 
            cc('CREATOR', 'TAMASIS v' + var.VERSION, 'Generator of this file'), 
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
    
    status = pyfits.BinTableHDU(status, None, 'STATUS')
    fits.append(status)
    fits.writeto(filename, clobber=True)

