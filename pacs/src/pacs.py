# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import gc
import numpy as np
import os
import pyfits
import re
import scipy
import time

from pyoperators import (Operator, HomothetyOperator, IdentityOperator,
                         MaskOperator, RoundOperator, ClipOperator, I)
from pyoperators.utils import isscalar, strenum, strplural, openmp_num_threads
from pyoperators.utils.mpi import MPI

from . import var
from tamasis.core import (Quantity, Map, Tod, MaskPolicy, Instrument, Pointing,
                          tmf, CompressionAverageOperator, ProjectionOperator,
                          PointingMatrix, create_fitsheader, filter_nonfinite)
from tamasis.mpiutils import distribute_observation, gather_fitsheader_if_needed
from tamasis.observations import Observation

__all__ = [ 'PacsInstrument',
            'PacsObservation',
            'PacsSimulation',
            'PacsConversionAdu',
            'PacsConversionVolts',
            'pacs_compute_delay',
            'pacs_get_psf',
            'pacs_preprocess',
            'step_scanline_masking',
            'step_deglitching']

p = var.path + '/pacs/'

class PacsInfo(object):
    CALFILE_BADP = p + 'PCalPhotometer_BadPixelMask_FM_v6.fits'
    CALFILE_TCST = p + 'PCalPhotometer_PhotTimeConstant_FM_v2.fits'
    CALFILE_GAIN = p + 'PCalPhotometer_Gain_FM_v1.fits'
    CALFILE_INVN = [p + 'PCalPhotometer_Invntt{0}_FM_v1.fits'.format(b)
                    for b in ('BS', 'BL', 'Red')]
    CALFILE_RESP = p + 'PCalPhotometer_Responsivity_FM_v6.fits'
    CALFILE_STDD = p + 'PCalPhotometer_Stddev_Tamasis_v1.fits'

del p

class PacsInstrument(Instrument):

    DEFAULT_RESOLUTION = {'blue':3.2, 'green':3.2, 'red':6.4}
    DETECTOR_DTYPE = [
        ('masked', np.bool8), ('removed', np.bool8),
        ('column', int), ('row', int), ('p', int), ('q', int), ('matrix', int),
        ('center', var.FLOAT_DTYPE, 2),
        ('corner', var.FLOAT_DTYPE, (4,2)),
        ('area', var.FLOAT_DTYPE), ('time_constant', var.FLOAT_DTYPE),
        ('flat_total', var.FLOAT_DTYPE), ('flat_detector', var.FLOAT_DTYPE),
        ('flat_optical', var.FLOAT_DTYPE)]
    SAMPLING_PERIOD = 0.024996 # in seconds

    info = PacsInfo()

    def __init__(self, band, active_fraction=0, fine_sampling_factor=1,
                 policy_bad_detector='mask', reject_bad_line=True,
                 comm=MPI.COMM_WORLD):

        """
        Set up the PACS instrument.

        Parameters
        ----------
        band : 'blue', 'green' or 'red'
              The PACS filter.
        fine_sampling_factor: integer
              Set this value to a power of two, which will increase the
              acquisition model's sampling frequency beyond 40Hz
        policy_bad_detector : 'keep', 'mask' or 'remove'
              This keyword controls the handling of the bad detectors.
              'keep' will handle bad detectors like valid ones.
              'mask' will  enforce a zero data value and a True mask value in
              the Tod returned by get_tod().
              'remove' will discard the bad detectors from the Tod returned by
              get_tod().
        reject_bad_line : boolean
              If True, the erratic line [11,16:32] (starting from 0) in
              the blue channel will be filtered out
        active_fraction : ratio of the geometric and reference detector area
              If set to 0, the value from the calibration file will be used

        """
        band = band.lower().strip()
        expected = 'blue', 'green', 'red'
        if band not in expected:
            raise ValueError("Invalid band '{0}'. Expected values are {1}.".
                format(band, strenum(expected)))
        self.band = band

        policy_bad_detector = policy_bad_detector.lower().strip()
        expected = 'keep', 'mask', 'remove'
        if policy_bad_detector not in expected:
            raise ValueError("Invalid detector policy '{0}'. Expected values ar"
                "e {1}.".format(policy_bad_detector, strenum(expected)))

        # get instrument information from calibration files
        # the call to pacs_info_instrument should be avoided, so that
        # all calibration files are defined and can be modified in PacsInfo.
        # Done: -responsivity
        #       -badpixel
        # TBD:  -arrayinstrument
        #       -subarrayarray
        #       -flatfield
        shape = (16,32) if band == 'red' else (32,64)
        junk, detector_center, detector_corner, detector_area, \
        distortion_yz, oflat, dflat, active_fraction, status = \
            tmf.pacs_info_instrument(band, shape[0], shape[1], active_fraction)
        if status != 0: raise RuntimeError()

        detector_bad = self.get_calibration('badpixel')
        if policy_bad_detector == 'keep':
            detector_bad[:] = False
        elif reject_bad_line and band != 'red':
            detector_bad[11,16:32] = True

        dflat[detector_bad] = np.nan
        detector_center = detector_center.T.swapaxes(0,1)
        detector_corner = detector_corner.T.swapaxes(0,1)

        Instrument.__init__(self, 'PACS/' + band.capitalize(), shape,
            detector_corner=detector_corner, default_resolution=\
            self.DEFAULT_RESOLUTION[band], dtype=self.DETECTOR_DTYPE,
            comm=comm)
        self.active_fraction = float(active_fraction)
        self.reject_bad_line = reject_bad_line
        self.fine_sampling_factor = fine_sampling_factor
        self.detector.masked = detector_bad
        pq = np.indices(shape)
        self.detector.row = pq[0] % 16
        self.detector.column = pq[1] % 16
        self.detector.p = pq[0]
        self.detector.q = pq[1]
        self.detector.matrix = pq[1] // 16 + 1
        if band != 'red':
            self.detector.matrix[0:16,:] += 4
        else:
            self.detector.matrix = 10 - self.detector.matrix
        self.detector.center = detector_center
        self.detector.corner = detector_corner
        self.detector.time_constant = self.get_calibration('timeconstant') # s
        self.detector.time_constant[detector_bad] = np.nan
        self.detector.area = detector_area # arcsec^2
        self.detector.flat_total = dflat * oflat
        self.detector.flat_detector = dflat
        self.detector.flat_optical = oflat
        self.distortion_yz = distortion_yz.T.view(dtype=[('y', \
            var.FLOAT_DTYPE), ('z',var.FLOAT_DTYPE)]).view(np.recarray)
        self.responsivity = self.get_calibration('responsivity')
        gain, offset = self.get_calibration('gain')
        self.adu_converter_gain = {'low':gain[0], 'high':gain[1],
                                   'nominal':gain[2]}
        self.adu_converter_offset = {'direct':offset[0], 'ddcs':offset[1]}
        self.adu_converter_min = {'direct':0, 'ddcs':-32768}
        self.adu_converter_max = {'direct':65535, 'ddcs':32767}

    def get_valid_detectors(self, masked=False):
        mask = ~self.detector.removed
        if masked:
            mask &= ~self.detector.masked
        i, j = np.mgrid[0:16,0:16]
        i = i.ravel()
        j = j.ravel()
        if self.band == 'red':
            i, j = (np.hstack([i,i]),
                    np.hstack([j,j+16]))
        else:
            i, j = (np.hstack([i + (k // 4) * 16 for k in range(8)]),
                    np.hstack([j + k*16 % 64 for k in range(8)]))
        mask = mask[(i,j)]
        return i[mask], j[mask]
    get_valid_detectors.__doc__ = Instrument.get_valid_detectors.__doc__

    def uv2yz(self, coords, chop=0):
        """
        Convert coordinates in the (u,v) plane into the (y,z) plane,
        assuming a chop angle.

        Parameters
        ----------
        coords : float array (last dimension is 2, for u and v)
            Coordinates in the (u,v) plane in millimeters.
        chop : float
            Chop angle.

        Returns
        -------
        coords_converted : float array (last dimension is 2, for y and z)
            Converted coordinates, in arc seconds.
        """
        coords = np.array(coords, float, order='c', copy=False)
        yz = np.empty_like(coords)
        distortion = self.distortion_yz.base.base.base
        tmf.pacs_uv2yz(coords.reshape((-1,2)).T, distortion, chop,
                       yz.reshape((-1,2)).T)
        yz *= 3600
        return yz
                
    def uv2ad(self, coords, pointing):
        """
        Convert coordinates in the (u,v) plane into celestial coordinates,
        assuming a pointing direction and a position angle.

        Parameters
        ----------
        coords : float array (last dimension is 2, for u and v)
            Coordinates in the (u,v) plane in millimeters.
        pointing : array with flexible dtype including 'ra', 'dec', 'pa', 'chop'
            Pointing direction, in degrees.

        Returns
        -------
        coords_converted : float array (last dimension is 2, for Ra and Dec)
            Converted coordinates, in degrees.

        Notes
        -----
        The routine is not accurate at the poles.
        """
        coords = np.array(coords, float, order='c', copy=False)
        ad = np.empty(pointing.shape + coords.shape, float)
        coords = coords.reshape((-1,2))
        distortion = self.distortion_yz.base.base.base
        tmf.pacs_uv2ad(coords.T, pointing['ra'].ravel(),
                       pointing['dec'].ravel(), pointing['pa'].ravel(),
                       pointing['chop'].ravel(), distortion,
                       ad.reshape((-1,)+coords.shape).T)
        for dim in range(pointing.ndim):
            ad = np.rollaxis(ad, 0, -1)
        return ad

    def yz2ad(self, coords, pointing):
        """
        Convert coordinates in the (y,z) plane into celestial coordinates,
        assuming a pointing direction and a position angle.

        Parameters
        ----------
        coords : float array (last dimension is 2)
            Coordinates in the (y,z) plane in arc seconds.
        pointing : array with flexible dtype including 'ra', 'dec', 'pa', 'chop'
            Pointing direction, in degrees.

        Returns
        -------
        coords_converted : float array (last dimension is 2, for Ra and Dec)
            Converted coordinates, in degrees.

        Notes
        -----
        The routine is not accurate at the poles.
        """
        return super(PacsInstrument, self).instrument2ad(coords, pointing)

    instrument2ad = uv2ad

    def instrument2xy_minmax(self, coords, pointing, header):
        """
        Return the minimum and maximum sky pixel coordinate values for a set of
        coordinates specified in the (u,v) frame.
        """
        coords = np.array(coords, float, order='c', copy=False)
        xmin, ymin, xmax, ymax, status = tmf.pacs.uv2xy_minmax(
            coords.reshape((-1,2)).T, pointing['ra'].ravel(),
            pointing['dec'].ravel(), pointing['pa'].ravel(),
            pointing['chop'].ravel(), self.distortion_yz.base.base.base,
            str(header).replace('\n',''))
        if status != 0:
            raise RuntimeError()
        return xmin, ymin, xmax, ymax

    def instrument2pmatrix_sharp_edges(self, coords, pointing, header, pmatrix,
                                       npixels_per_sample):
        """
        Return the dense pointing matrix whose values are intersection between
        detectors and map pixels.
        """
        coords = coords.reshape((-1,2))
        ra = pointing['ra'].ravel()
        dec = pointing['dec'].ravel()
        pa = pointing['pa'].ravel()
        chop = pointing['chop'].ravel()
        masked = pointing['masked'].ravel().view(np.int8)
        if pmatrix.size == 0:
            # f2py doesn't accept zero-sized opaque arguments
            pmatrix = np.empty(1, np.int64)
        else:
            pmatrix = pmatrix.ravel().view(np.int64)
        header = str(header).replace('\n','')
        distortion = self.distortion_yz.base.base.base

        new_npixels_per_sample, out, status = tmf.pacs.uv2pmatrix_sharp_edges(
            coords.T, ra, dec, pa, chop, masked, distortion, header, pmatrix,
            npixels_per_sample)
        if status != 0: raise RuntimeError()

        return new_npixels_per_sample, out

    def get_pointing_matrix(self, pointing, header, npixels_per_sample=0,
                            method='sharp', downsampling=False,
                            compression_factor=None, delay=0., units=None,
                            derived_units=None, comm=MPI.COMM_WORLD):
        """
        Return the PACS-specific pointing matrix for a given set of pointings.

        Parameters
        ----------
            pointing : Pointing
                The pointing containing the astrometry of the array center.
            header : pyfits.Header
                The map FITS header
            npixels_per_sample : int
                Maximum number of sky pixels intercepted by a detector.
                By setting 0 (the default), the actual value will be determined
                automatically.
            method : string
                'sharp' : the intersection of the sky pixels and the detectors
                          is computed assuming that the transmission outside
                          the detector is zero and one otherwise (sharp edge
                          geometry)
                'nearest' : the value of the sky pixel closest to the detector
                            center is taken as the sample value, assuming
                            surface brightness conservation.
            downsampling : boolean
                If True, return a pointing matrix downsampled by the compression
                factor. Otherwise return a pointing matrix sampled at 40Hz,
                unless the instrument has been setup with a higher sampling rate
                via the 'fine_sampling_factor' keyword.
            compression_factor : int
                The compression factor, used if downsampling is True.
            delay : float
                Instrument clock lag wrt the spacecraft clock, in ms.
            comm : mpi4py.MPI.Comm
                Map communicator, only used if the input FITS header is local.

        """
        if method is None:
            method = 'sharp'
        method = method.lower()
        choices = ('nearest', 'sharp')
        if method not in choices:
            raise ValueError("Invalid method '" + method + \
                "'. Expected values are " + strenum(choices) + '.')

        header = gather_fitsheader_if_needed(header, comm=comm)
        shape_input = tuple(header['NAXIS' + str(i+1)]
                            for i in range(header['NAXIS']))[::-1]
        if np.product(shape_input) > np.iinfo(np.int32).max:
            raise RuntimeError('The map is too large: pixel indices cannot be s'
                'stored using 32 bits: {0}>{1}'.format(np.product(shape_input),
                np.iinfo(np.int32).max))

        nvalids = int(np.sum(~pointing.removed))
        if not downsampling:
            if compression_factor is None:
                raise ValueError('The compression factor is missing.')
            nvalids *= self.fine_sampling_factor * compression_factor

        ndetectors = self.get_ndetectors()
        shape = (ndetectors, nvalids, npixels_per_sample)
        info = {'header':header,
                'method':method,
                'units':units,
                'derived_units':derived_units,
                'downsampling':downsampling}

        try:
            pmatrix = PointingMatrix.empty(shape, shape_input, info=info)
        except MemoryError:
            gc.collect()
            pmatrix = PointingMatrix.empty(shape, shape_input, info=info,
                                           verbose=False)

        if pmatrix.size == 0:
            # f2py doesn't accept zero-sized opaque arguments
            pmatrix_ = np.empty(1, np.int64)
        else:
            pmatrix_ = pmatrix.ravel().view(np.int64)

        new_npixels_per_sample, status = tmf.pacs_pointing_matrix(
            self.band,
            nvalids,
            compression_factor,
            delay,
            self.fine_sampling_factor,
            not downsampling,
            np.asfortranarray(pointing.time),
            np.asfortranarray(pointing.ra),
            np.asfortranarray(pointing.dec),
            np.asfortranarray(pointing.pa),
            np.asfortranarray(pointing.chop),
            np.asfortranarray(pointing.masked, np.int8),
            np.asfortranarray(pointing.removed,np.int8),
            method,
            np.asfortranarray(self.detector.removed, np.int8),
            self.get_ndetectors(),
            self.detector.center.swapaxes(0,1).copy().T,
            self.detector.corner.swapaxes(0,1).copy().T,
            np.asfortranarray(self.detector.area),
            self.distortion_yz.base.base.base,
            npixels_per_sample,
            str(header).replace('\n',''),
            pmatrix_)
        if status != 0: raise RuntimeError()

        # if the actual number of pixels per sample is greater than
        # the specified one, redo the computation of the pointing matrix
        if new_npixels_per_sample <= npixels_per_sample:
            return pmatrix

        del pmatrix_, pmatrix
        return self.get_pointing_matrix(pointing, header,
            new_npixels_per_sample, method, compression_factor= \
            compression_factor, delay=delay, units=units,
            derived_units=derived_units, downsampling=downsampling, comm=comm)

    def get_calibration(self, name):
        name = name.lower()
        expected = ('badpixel', 'gain', 'invntt', 'responsivity',
                    'timeconstant', 'stddev')
        if name not in expected:
            raise ValueError("Invalid calibration type '{0}'. Expected values a"
                             "re {1}.".format(name, strenum(expected)))

        if  name == 'badpixel':
            hdus = self._get_calfile(self.info.CALFILE_BADP, 1.0)
            return Map(hdus[self.band].data, dtype=np.bool8, origin='upper')

        if name == 'gain':
            hdus = self._get_calfile(self.info.CALFILE_GAIN, 1.0)
            return hdus['photGain'].data, hdus['photOffset'].data

        if name == 'invntt':
            filter = {'blue':0, 'green':1, 'red':2}[self.band]
            filename = self.info.CALFILE_INVN[filter]
            hdus = self._get_calfile(filename, 1.0)
            ncorrelations = hdus[1].header['NAXIS1'] - 1
            data, status = tmf.pacs_read_filter_uncorrelated(filename,
                ncorrelations, self.get_ndetectors(),
                np.asfortranarray(self.detector.removed, np.int8))
            if status != 0: raise RuntimeError()
            return data.T

        if name == 'responsivity':
            hdus = self._get_calfile(self.info.CALFILE_RESP, 2.0)
            filter = {'blue':6, 'green':4, 'red':2}[self.band]
            return Quantity(hdus[filter].data, 'V/Jy')

        if name == 'timeconstant':
            hdus = self._get_calfile(self.info.CALFILE_TCST, None)
            return hdus[{'red':1,'green':2,'blue':3}[self.band]].data / 1000 # s

        if name == 'stddev':
            hdus = self._get_calfile(self.info.CALFILE_STDD, None)
            lengths = hdus[0].data
            stddevs = hdus[self.band].data
            return lengths, stddevs

    def _get_calfile(self, filename, versions):
        """ Return FITS file as an HDU list and check format version. """
        if isscalar(versions):
            versions = (versions,)
        if not os.path.isfile(filename):
            raise IOError("The calibration '{0}' file does not exist or is "
                          "not a valid file.".format(filename))

        hdus = pyfits.open(filename)
        header = hdus[0].header
        version = _hcss_fits_keyword(header, 'formatVersion', None)
        if version is not None:
            version = float(version)
        if version not in versions:
            raise ValueError("Invalid format version '{0}'. Expected version is"
                             " {1}.".format(version, strenum(versions)))
        return hdus

class PacsBase(Observation):
    """ Abstract class for PacsObservation and PacsSimulation. """

    ACCELERATION = Quantity(4., '"/s^2')
    PSF_FWHM = {'blue':5.2, 'green':7.7, 'red':12.}
    POINTING_DTYPE = [('ra', var.FLOAT_DTYPE), ('dec', var.FLOAT_DTYPE),
                      ('pa', var.FLOAT_DTYPE), ('chop', var.FLOAT_DTYPE),
                      ('time', var.FLOAT_DTYPE), ('info', np.int64),
                      ('masked', np.bool8), ('removed', np.bool8)]
    SAMPLING_PERIOD = PacsInstrument.SAMPLING_PERIOD

    def get_derived_units(self):
        volt = Quantity(1/self.instrument.responsivity, 'Jy')
        return (
            {
                'detector_reference': Quantity(1 / \
                    self.instrument.active_fraction, 'detector'),
                'detector[rightward]' : Quantity(self.pack(
                    self.instrument.detector.area), 'arcsec^2'),
                'V' : volt,
            },{
                'V' : volt
            }
            )

    def get_filter_uncorrelated(self):
        """
        Return the inverse noise time-time correlation coefficients.
        """
        return self.instrument.get_calibration('invntt')

    def get_pointing_matrix(self, header, npixels_per_sample=0, method=None,
                            downsampling=False, section=None, 
                            comm=MPI.COMM_WORLD):
        """
        Return the PACS-specific pointing matrix for the observation.

        If the observation has several slices, as many pointing matrices
        are returned in a list.

        Parameters
        ----------
            header : pyfits.Header
                The map FITS header
            npixels_per_sample : int
                Maximum number of sky pixels intercepted by a detector.
                By setting 0 (the default), the actual value will be determined
                automatically.
            method : string
                'sharp' : the intersection of the sky pixels and the detectors
                          is computed assuming that the transmission outside
                          the detector is zero and one otherwise (sharp edge
                          geometry)
                'nearest' : the value of the sky pixel closest to the detector
                            center is taken as the sample value, assuming
                            surface brightness conservation.
            section : PacsObservation.slice
                If specified, return the pointing matrix for a specific slice.
            downsampling : boolean
                If True, return a pointing matrix downsampled by the compression
                factor. Otherwise return a pointing matrix sampled at 40Hz,
                unless the instrument has been setup with a higher sampling rate
                via the 'fine_sampling_factor' keyword.
            comm : mpi4py.MPI.Comm
                MPI communicator, only used if the input FITS header is local.

        """
        if section is None:
            return super(PacsBase, self).get_pointing_matrix(header,
                         npixels_per_sample, method, downsampling, comm=comm)

        return super(PacsBase, self).get_pointing_matrix(header,
            npixels_per_sample, method, downsampling, compression_factor=
            section.compression_factor, delay=section.delay, units=
            ('/detector', '/pixel'), derived_units=self.get_derived_units(),
            section=section, comm=comm)

    def get_random(self, flatfielding=True, subtraction_mean=True):
        """
        Return noise data from a random slice of a real pointed observation.

        Currently, the required duration must be less than that of the 
        real observation (around 3 hours).
        """
        if any([slice.compression_factor not in [4,8] for slice in self.slice]):
            raise NotImplementedError('The compression factor must be 4 or 8.')

        path = os.path.join(var.path, 'pacs')
        files = {
            'blue' : ('1342182424_blue_PreparedFrames.fits.gz', 10000,100000),
            'green': ('1342182427_green_PreparedFrames.fits.gz', 10000, 100000),
            'red': ('1342182427_red_PreparedFrames.fits.gz', 7400, 116400) }
        file, istart, iend = files[self.instrument.band]
        data = pyfits.open(os.path.join(path, file))['Signal'] \
                     .data[:,:,istart:iend]
        data /= self.instrument.active_fraction * self.instrument.responsivity
        nsamples = self.get_nsamples()
        nsamples_tot = np.sum(nsamples*self.slice.compression_factor / 4)
        if nsamples_tot > data.shape[-1]:
            raise ValueError('There is not enough noise data for this observati'
                             'on.')

        irandom = np.random.random_integers(0,iend-istart-nsamples_tot)
        result = Tod(data[:,:,irandom:irandom+nsamples_tot],
                     unit='Jy/detector',
                     derived_units=self.get_derived_units()[0])
        if subtraction_mean:
            result.T[:] -= np.mean(result, axis=1)
        if flatfielding:
            result.T[:] /= self.pack(self.instrument.detector.flat_detector)

        compression = CompressionAverageOperator(
            self.slice.compression_factor / 4)
        result = compression(result)
        return result

    def get_detector_stddev(self, length=0):
        """Returns detector's standard deviation from calibration file"""
        lengths, stddevs = self.instrument.get_calibration('stddev')
        if length == 0:
            return self.pack(stddevs[...,-1])
        if length < lengths[0]:
            raise ValueError('The value of the median filtering length should b'
                             'e greater than ' + str(lengths[0]) + '.')
        i = 1
        while i < len(lengths) - 2 and length >= lengths[i]:
            i += 1
        d = (np.log(length) - np.log(lengths[i-1])) / \
            (np.log(lengths[i]) - np.log(lengths[i-1]))
        s = np.minimum((1-d) * stddevs[...,i-1] + d * stddevs[...,i],
                       stddevs[...,-1])
        s[np.isnan(s)] = np.inf
        return self.pack(s)

    def get_nfinesamples(self):
        """
        Return the number of valid samplings for each slice, by taking
        into account compression factor and fine sampling.
        They are those whose self.pointing is not removed.
        """
        return tuple(np.asarray(self.get_nsamples()) * \
                     self.slice.compression_factor * \
                     self.instrument.fine_sampling_factor)

    def save(self, filename, tod, fitskw={}):
        
        if np.rank(tod) != 3:
            tod = self.unpack(tod)

        nsamples = np.sum(~self.pointing.removed)
        if nsamples != tod.shape[-1]:
            raise ValueError("The input Tod has a number of samples'" + \
                str(tod.shape[-1]) + "' incompatible with that of this observat"
                "ation '" + str(nsamples) + "'.")

        _write_status(self, filename, fitskw=fitskw)

        if tod.header is None:
            header = create_fitsheader(fromdata=tod, extname='Signal')
        else:
            header = tod.header.copy()
            header.update('EXTNAME', 'Signal')

        if tod.unit != '':
            header.update('QTTY____', tod.unit)

        pyfits.append(filename, tod, header)

        if tod.mask is not None:
            _write_mask(self, tod.mask, filename)
    save.__doc__ = Observation.save.__doc__

    @classmethod
    def create_scan(cls, center, length, step=148, sampling_period=None,
                    speed=20, acceleration=ACCELERATION, nlegs=3,
                    angle=0, instrument_angle=45, compression_factor=4,
                    cross_scan=True):
        """
        Return a sky scan for the PACS instrument.

        The output is a Pointing instance that can be handed to PacsSimulation
        to create a simulation.

        Parameters
        ----------
        center : tuple
            (Right Ascension, declination) of the scan center.
        length : float
            Length of the scan lines, in arcseconds.
        step : float
            Separation between scan legs, in arcseconds.
        sampling_period : float
            Duration between two pointings.
        speed : float
            Scan speed, in arcsec/s.
        acceleration : float
            Acceleration, in arcsec/s^2.
        nlegs : integer
            Number of scan legs.
        angle : float
            Angle between the scan line direction and the North minus
            90 degrees, in degrees.
        instrument_angle : float
            Angle between the scan line direction and the PACS instrument Z-axis
            (called array-to-map angle in the PACS user's manual).
        compression_factor : int
            Compression factor.
        cross_scan : boolean
            If true, a cross-scan is appended to the pointings.
        """

        if sampling_period is None:
            if int(compression_factor) not in (1, 4, 8):
                raise ValueError("Input compression_factor must be 1, 4 or 8.")
            sampling_period = PacsInstrument.SAMPLING_PERIOD*compression_factor

        scan = Observation.create_scan(center, length, step, sampling_period,
            speed, acceleration, nlegs, angle, instrument_angle, cross_scan,
            PacsBase.POINTING_DTYPE)

        first = (scan.info == scan.INSCAN).nonzero()[0][0]
        scan.info[:first] = scan.OTHER
        return scan

    def __str__(self):

        # some helpers
        def same(a):
            a = np.asarray(a)
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
            return ((('no ' if len(m) == 0 else '') + strplural(s + 'activated '
                'mask', len(m), False, ': ' + ', '.join(m))).capitalize() \
                for s,m in (('',a), ('de',d)))

        nthreads = openmp_num_threads()
        ndetectors = self.get_ndetectors()
        sp = len('Info ')*' '
        unit = 'unknown' if self.slice[0].unit == '' else self.slice[0].unit
        if self.instrument.comm.size > 1:
            mpistr = 'Process ' + str(self.instrument.comm.rank+1) + '/' + \
                     str(self.instrument.comm.size) + ', '
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
        scan_speed = self.slice.scan_speed
        print_each_scan_speed = not same(scan_speed)
        if not print_each_scan_speed:
            info.append('scan speed ' + str(scan_speed[0]) + ' "/s')
        
        # unit
        if self.__class__.__name__ == 'PacsObservation':
            info.append('unit is ' + unit)

        result += sp + ', '.join(info) + '\n'

        # mask information
        if self.__class__.__name__ == 'PacsObservation':        
            homogeneous = 'nmasks' not in self.slice.dtype.names or \
                same([set(m) for m in self.slice.mask_name])
            if homogeneous:
                (a,d) = ad_masks(self.slice, 0)
                result += sp + a + '\n'
                result += sp + d + '\n'
            else:
                result += sp + 'The masks of the observations are heterogeneous'
                result += '\n'

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
                result += sp + '     ' + a + '\n'
                result += sp + '     ' + d + '\n'

            if print_each_mode:
                result += sp + '     Mode:        ' + slice.mode + '\n'
            if print_each_compression:
                result += sp + '     Compression: x' + \
                    str(slice.compression_factor) + '\n'
            if print_each_scan_speed:
                result += sp + '     Scan speed:  ' + str(scan_speed[islice]) +\
                          ' "/s\n'
            result += sp + '     In-scan:     ' + \
                _str_policy(p.info == Pointing.INSCAN) + '\n'
            result += sp + '     Turnaround:  ' + \
                _str_policy(p.info == Pointing.TURNAROUND) + '\n'
            result += sp + '     Other:       ' + \
                _str_policy(p.info == Pointing.OTHER)
            if islice + 1 < self.slice.size:
                result += '\n'

        return result


#-------------------------------------------------------------------------------


class PacsObservation(PacsBase):
    """
    Class which encapsulates handy information about the PACS instrument and 
    the observations to be processed.
    """
    def __init__(self, filename, fine_sampling_factor=1,
                 policy_bad_detector='mask', reject_bad_line=True,
                 policy_inscan='keep', policy_turnaround='keep',
                 policy_other='remove', policy_invalid='mask',
                 active_fraction=0, delay=0., calblock_extension_time=None,
                 comm=MPI.COMM_WORLD):
        """
        Parameters
        ----------
        filenames : string or array of string
              This argument indicates the filenames of the observations
              to be processed. The files must be FITS files, as saved by
              the HCSS method FitsArchive.
        fine_sampling_factor : integer
              Set this value to a power of two, which will increase the
              acquisition model's sampling frequency beyond 40Hz
        policy_bad_detector : 'keep', 'mask' or 'remove'
              This keyword controls the handling of the bad detectors.
              'keep' will handle bad detectors like valid ones.
              'mask' will  enforce a zero data value and a True mask value in
              the Tod returned by get_tod().
              'remove' will discard the bad detectors from the Tod returned by
              get_tod().
        reject_bad_line : boolean
              If True, the erratic line [11,16:32] (starting from 0) in
              the blue channel will be filtered out
        policy_inscan : 'keep', 'mask' or 'remove'
              This sets the policy for the in-scan frames
        policy_turnaround: 'keep', 'mask' or 'remove'
              This sets the policy for the turnaround frames
        policy_other : 'keep', 'mask' or 'remove'
              This sets the policy for the other frames
        policy_invalid : 'keep', 'mask' or 'remove'
              This sets the policy for the invalid frames
        active_fraction : ratio of the geometric and reference detector area
              If set to 0, the value from the calibration file will be used
        delay : instrument clock lag wrt the spacecraft clock, in ms
        calblock_extension_time: float
              Elapsed time in seconds during which samples after a calibration
              block are flagged as 'other', similarly to the calibration block
              itself. Default is 20s for the blue channel, 40s for the red one.
        comm : mpi4py.MPI.Comm
              The TOD communicator.

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
        if not isinstance(filename, (list,tuple)):
            filename = (filename,)
        filename_, nfilenames = _files2tmf(filename)
        active_fraction = float(active_fraction)

        band, transparent_mode, status = tmf.pacs_info_instrument_init(
            filename_, nfilenames)
        if status != 0: raise RuntimeError()

        # set up the instrument
        self.instrument = PacsInstrument(band, active_fraction=active_fraction,
            fine_sampling_factor=fine_sampling_factor,
            policy_bad_detector=policy_bad_detector,
            reject_bad_line=reject_bad_line, comm=comm)
        band = self.instrument.band

        # get work load according to detector policy and MPI process
        self.instrument.detector.removed, filename = _get_workload(band,
            filename, self.instrument.detector.masked, policy_bad_detector,
            transparent_mode, comm)
        filename_, nfilenames = _files2tmf(filename)

        # get the observations and detector mask for the current processor
        nsamples_all, status = tmf.pacs_info_observation_init(filename_,
            nfilenames)
        if status != 0: raise RuntimeError()

        # frame policy
        policy = MaskPolicy('inscan,turnaround,other,invalid', (policy_inscan,
            policy_turnaround, policy_other, policy_invalid), 'Frame Policy')

        # store observation information
        if calblock_extension_time is None:
            calblock_extension_time = 40. if band == 'red' else 20.
        obsid, mode, compression_factor, unit, ra, dec, instrument_angle, \
            scan_angle, scan_length, scan_step, scan_nlegs, frame_time, \
            frame_ra, frame_dec, frame_pa, frame_chop, frame_info, \
            frame_masked, frame_removed, nmasks, mask_name_flat, \
            mask_activated, status = tmf.pacs_info_observation(filename_,
            nfilenames, np.array(policy, np.int32),np.sum(nsamples_all),
            calblock_extension_time)
        if status != 0: raise RuntimeError()

        flen_value = len(unit) // nfilenames
        mode = [mode[i*flen_value:(i+1)*flen_value].strip() \
                for i in range(nfilenames)]
        unit = [unit[i*flen_value:(i+1)*flen_value].strip() \
                for i in range(nfilenames)]

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

        mask_name = np.zeros((nfilenames, nmasks_max), 'S'+str(mask_len_max))
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
            ('start', int),
            ('stop', int),
            ('obsid', int),
            ('mode', 'S32'),
            ('compression_factor', int),
            ('delay', float),
            ('unit', 'S32'),
            ('ra', float),
            ('dec', float),
            ('instrument_angle', float),
            ('scan_angle', float),
            ('scan_length', float),
            ('scan_nlegs', int),
            ('scan_step', float),
            ('scan_speed', float),
            ('nmasks', int),
            ('mask_name', 'S'+str(mask_len_max), nmasks_max),
            ('mask_activated', bool, nmasks_max),
            ])

        regex = re.compile(r'(.*?)(\[[0-9]*:?[0-9]*\])? *$')
        for ifile, file in enumerate(filename):
            match = regex.match(file)
            self.slice[ifile].filename = match.group(1)

        self.slice.nsamples_all = nsamples_all
        self.slice.start[0] = 0
        self.slice.start[1:] = np.cumsum(nsamples_all)[:-1]
        self.slice.stop = np.cumsum(nsamples_all)
        self.slice.nfinesamples = nsamples_all * compression_factor * \
            fine_sampling_factor
        self.slice.obsid = obsid
        self.slice.mode = mode
        self.slice.compression_factor = compression_factor
        self.slice.delay = delay
        self.slice.unit = unit
        self.slice.ra = ra
        self.slice.dec = dec
        self.slice.instrument_angle = instrument_angle
        self.slice.scan_angle = scan_angle
        self.slice.scan_length = scan_length
        self.slice.scan_nlegs = scan_nlegs
        self.slice.scan_step = scan_step
        self.slice.nmasks = nmasks
        self.slice.mask_name      = mask_name
        self.slice.mask_activated = mask_activated[:,0:nmasks_max]
        
        # store pointing information
        self.pointing = Pointing((frame_ra, frame_dec, frame_pa), frame_time,
            frame_info, frame_masked, frame_removed, dtype=self.POINTING_DTYPE)
        self.pointing.chop = frame_chop
        self.slice.scan_speed = [_scan_speed(self.pointing[s.start:s.stop]) \
                                 for s in self.slice]

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

        By default, all activated masks will be combined.
        """

        if raw:
            flatfielding = False
            subtraction_mean = False

        act_masks = set([m for slice in self.slice \
                         for i, m in enumerate(slice.mask_name) \
                         if m and slice.mask_activated[i]])
        dea_masks = set([m for slice in self.slice \
                         for i, m in enumerate(slice.mask_name) \
                         if m and not slice.mask_activated[i]])
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

        sel_masks = ','.join(sorted(sel_masks))

        ndetectors = self.get_ndetectors()
        nsamples_tot = int(np.sum(self.get_nsamples()))
        tod = Tod.empty((ndetectors, nsamples_tot),
                        mask=np.empty((ndetectors, nsamples_tot), np.bool8),
                        unit=self.slice[0].unit,
                        derived_units=self.get_derived_units()[0])

        time0 = time.time()
        first = 1
        for s in self.slice:
            masked = self.pointing.masked[s.start:s.stop]
            removed = self.pointing.removed[s.start:s.stop]
            last = first + int(np.sum(~removed)) - 1
            status = tmf.pacs_read_tod(tod.T, tod.mask.view(np.int8).T,
                first, last,
                self.instrument.band,
                s.filename,
                np.asfortranarray(masked, np.int8),
                np.asfortranarray(removed, np.int8),
                np.asfortranarray(self.instrument.detector.masked, np.int8),
                np.asfortranarray(self.instrument.detector.removed, np.int8),
                np.asfortranarray(self.instrument.detector.flat_detector),
                flatfielding,
                subtraction_mean,
                sel_masks)
            if status != 0: raise RuntimeError()
            first = last + 1
        print('Reading timeline... {0:.2f}'.format(time.time()-time0))

        # IFORT stores .true. as -1, it conflicts with numpy algebra on booleans
        np.abs(tod.mask, tod.mask)

        filter_nonfinite(tod, out=tod)

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
                except IndexError:
                    pass
            
            status.append(s.base)

        # check number of records
        if np.any([len(s) for s in status] != self.slice.nsamples_all):
            raise ValueError("The status has a number of records '" + \
                str([len(s) for s in status]) + "' incompatible with that of th"
                "e pointings '" + str(self.slice.nsamples_all) + "'.")

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
    def __init__(self, pointings, band, mode=None, fine_sampling_factor=1,
                 policy_bad_detector='mask', reject_bad_line=True,
                 policy_inscan='keep', policy_turnaround='keep',
                 policy_other='keep', policy_invalid='keep', active_fraction=0,
                 delay=0., comm=MPI.COMM_WORLD):

        if comm.size > 1:
            raise NotImplementedError()

        # get work load according to detector policy and MPI process
        self.instrument = PacsInstrument(band, active_fraction=active_fraction,
            fine_sampling_factor=fine_sampling_factor,
            policy_bad_detector=policy_bad_detector,
            reject_bad_line=reject_bad_line, comm=comm)
        band = self.instrument.band

        if not isinstance(pointings, (list,tuple)):
            pointings = (pointings,)

        # get compression factor from input pointing
        deltas = [p.time[1:] - p.time[0:-1] for p in pointings if p.size > 1]
        if len(deltas) == 0:
            delta = PacsInstrument.SAMPLING_PERIOD
        else:
            deltas = [np.median(d) for d in deltas]
            delta = np.median(np.hstack(deltas))
            if not np.allclose(deltas, delta, rtol=0.01):
                raise ValueError('The input pointings do not have the same samp'
                                 'ling period: '+str(deltas)+'.')

        compression_factor = int(np.round(delta/PacsInstrument.SAMPLING_PERIOD))
        if compression_factor <= 0:
            raise ValueError('Invalid time in pointing argument. Use PacsInstru'
                             'ment.SAMPLING_PERIOD.')
        if np.abs(delta / compression_factor - PacsInstrument.SAMPLING_PERIOD) \
           > 0.01 * PacsInstrument.SAMPLING_PERIOD:
            print('Warning: the pointing time has an unexpected sampling rate. '
                  'Assuming a compression factor of {0}.'.format(
                  compression_factor))

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
                    "'. Expected values are " + strenum(choices) + '.')
            if mode == 'prime' and compression_factor != 4 or \
               mode == 'parallel' and compression_factor != (4 if band == 'red'\
                    else 8) or \
               mode == 'transparent' and compression_factor != 1:
                raise ValueError("The observing mode '" + mode + "' in the " + \
                    band + " band is incompatible with the compression factor '"
                    "" + str(compression_factor) + "'.")

        nslices = len(pointings)
        nsamples = tuple([p.size for p in pointings])

        # concatenate pointings
        pointing = Pointing.empty(sum(nsamples), dtype=self.POINTING_DTYPE)
        d = 0
        for p, n in zip(pointings, nsamples):
            names = (set(p.dtype.names) | set(p.__dict__.keys())) & \
                     set(pointing.dtype.names)
            for name in names:
                pointing[name][d:d+n] = p[name]
            d += n

        # store pointing information
        pointing.masked  |= \
            (policy_inscan     == 'mask')   & \
                (pointing.info == Pointing.INSCAN) | \
            (policy_turnaround == 'mask')   & \
                (pointing.info == Pointing.TURNAROUND) | \
            (policy_other      == 'mask')   & \
                (pointing.info == Pointing.OTHER)
        pointing.removed |= \
            (policy_inscan     == 'remove') & \
                (pointing.info == Pointing.INSCAN) | \
            (policy_turnaround == 'remove') & \
                (pointing.info == Pointing.TURNAROUND) | \
            (policy_other      == 'remove') & \
                (pointing.info == Pointing.OTHER)
        self.pointing = pointing

        self.slice = np.recarray(nslices, dtype=[
            ('filename', 'S256'),
            ('nsamples_all', int),
            ('start', int),
            ('stop', int),
            ('obsid', int),
            ('mode', 'S32'),
            ('compression_factor', int),
            ('delay', float),
            ('unit', 'S32'),
            ('ra', float),
            ('dec', float),
            ('instrument_angle', float),
            ('scan_angle', float),
            ('scan_length', float),
            ('scan_nlegs', int),
            ('scan_step', float),
            ('scan_speed', float)])

        self.slice.filename = 'simulation'
        self.slice.nsamples_all = nsamples
        self.slice.start[0] = 0
        self.slice.start[1:] = np.cumsum(nsamples)[:-1]
        self.slice.stop = np.cumsum(nsamples)
        self.slice.obsid = 0
        self.slice.mode = mode
        self.slice.compression_factor = compression_factor
        self.slice.delay = delay
        self.slice.unit = ''
        
        for s, p in zip(self.slice, pointings):
            if p.header is not None:
                for field in ('ra', 'dec', 'instrument_angle', 'scan_angle',
                              'scan_nlegs', 'scan_length', 'scan_step'):
                    s[field] = p.header[field] if field in p.header else 0
            s.scan_speed = _scan_speed(p)
        
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
            raise ValueError('The pointing and slice attribute are incompatible'
                             '. This should not happen.')

        status = np.recarray(p.size, dtype=[('BBID', np.int64), ('FINETIME',
            np.int64), ('BAND', 'S2'),('CHOPFPUANGLE', np.float64),
            ('RaArray', np.float64), ('DecArray', np.float64),
            ('PaArray', np.float64)])
        band = self.instrument.band
        status.BAND = {'blue':'BS', 'green':'BL', 'red':'R '}[band]
        status.FINETIME = np.round(p.time*1000000.)
        status.RaArray  = p.ra
        status.DecArray = p.dec
        status.PaArray  = p.pa
        status.CHOPFPUANGLE = 0. if not hasattr(p, 'chop') else p.chop
        status.BBID = 0
        status.BBID[p.info == Pointing.INSCAN] = 0xcd2 << 16
        status.BBID[p.info == Pointing.TURNAROUND] = 0x4000 << 16
        self._status = status

        return status


#-------------------------------------------------------------------------------


class PacsMultiplexing(Operator):
    """
    Performs the multiplexing of the PACS subarrays.
    The subarray columns are read one after the other, in a 0.025s cycle (40Hz).
    """
    def __init__(self, obs, shapein=None, description=None):
        Operator.__init__(self,
                          shapein=shapein,
                          classin=Tod,
                          description=description)
        self.fine_sampling_factor = obs.instrument.fine_sampling_factor
        self.ij = obs.instrument.ij

    def direct(self, input, output):
        tmf.pacs_multiplexing_direct(input.T, output.T,
                                     self.fine_sampling_factor, self.ij)

    def transpose(self, input, output):
        tmf.pacs_multiplexing_transpose(input.T, output.T,
                                        self.fine_sampling_factor, self.ij)

    def reshapein(self, shapein):
        return shapein[:-1] + (shapein[-1] // self.fine_sampling_factor,)

    def reshapeout(self, shapeout):
        return shapeout[:-1] + (shapeout[-1] * self.fine_sampling_factor,)

    def validatein(self, shapein):
        if shapein[-1] % self.fine_sampling_factor != 0:
            raise ValueError("The input timeline size '{0}' is incompatible wit"
                "h the fine sampling factor '{1}'.".format(shapein[-1],
                self.fine_sampling_factor))


#-------------------------------------------------------------------------------


def PacsConversionAdu(obs, gain='nominal', offset='direct'):
    """
    Returns the operator which converts PACS timelines from Volts to ADU.
    """
    if gain not in obs.instrument.adu_converter_gain:
        raise ValueError("Invalid gain. Expected values are 'low', 'high' or 'n"
                         "ominal'.")
    if offset not in obs.instrument.adu_converter_offset:
        raise ValueError("Invalid offset. Expected values are 'direct' or 'ddcs"
                         "'.")

    g = obs.instrument.adu_converter_gain[gain]
    o = obs.instrument.adu_converter_offset[offset]
    a = obs.instrument.adu_converter_min[offset]
    z = obs.instrument.adu_converter_max[offset]

    return np.product([
        ClipOperator(a, z), # ADU converter saturation
        RoundOperator(),    # ADU converter rounding
        I + o,              # ADU converter offset
        1 / g])             # ADU converter gain


#-------------------------------------------------------------------------------


def PacsConversionVolts(obs):
    """Jy-to-Volts conversion."""
    op = HomothetyOperator(obs.instrument.responsivity / \
                           obs.instrument.active_fraction)
    op.__name__ = 'Jy-to-Volts conversion'
    return op


#-------------------------------------------------------------------------------


def pacs_preprocess(obs, tod,
                    projection_method='sharp',
                    header=None,
                    downsampling=False,
                    npixels_per_sample=0,
                    deglitching_hf_length=20,
                    deglitching_nsigma=5.,
                    hf_length=30000,
                    transparent_mode_compression_factor=1):
    """
    deglitch, filter and potentially compress if the observation is in
    transparent mode
    """
    projection = ProjectionOperator(obs,
                            method='sharp',
                            header=header,
                            downsampling=True,
                            npixels_per_sample=npixels_per_sample)
    tod_filtered = filter_median(tod, deglitching_hf_length)
    tod.mask = deglitch_l2mad(tod_filtered,
                              projection,
                              nsigma=deglitching_nsigma)
    tod = filter_median(tod, hf_length)
    masking = MaskOperator(tod.mask)
    tod = masking(tod)
    
    # get the proper projector if necessary
    if projection is None or projection_method != 'sharp' or \
       not downsampling and np.any(obs.slice.compression_factor * \
       obs.instrument.fine_sampling_factor > 1):
        projection = ProjectionOperator(obs,
                                method=projection_method,
                                downsampling=downsampling,
                                header=header,
                                npixels_per_sample=npixels_per_sample)

    # bail out if not in transparent mode
    if all(obs.slice[0].compression_factor != 1) or \
       transparent_mode_compression_factor == 1:
        if not downsampling:
            compression = CompressionAverageOperator(
                obs.slice.compression_factor)
            model = compression * projection
        else:
            model = projection
        map_mask = model.T(tod.mask)
        model = masking * model
        return tod, model, mapper_naive(tod, model), map_mask

    # compress the transparent observation
    compression = CompressionAverageOperator(
        transparent_mode_compression_factor)
    todc = compression(tod)
    mask = compression(tod.mask)
    mask[mask != 0] = 1
    todc.mask = np.array(mask, dtype='uint8')
    maskingc = MaskOperator(todc.mask)

    model = compression * projection
    map_mask = model.T(tod.mask)
    model = masking * model

    return todc, model, mapper_naive(todc, model), map_mask


#-------------------------------------------------------------------------------


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
                         strenum(choices) + '.')

    kind = kind.lower()
    choices = ('airy', 'gaussian', 'calibration')
    if kind not in choices:
        raise ValueError("Invalid PSF kind '" + kind + \
            "'. Expected values are " + strenum(choices) + '.')

    if kind in ('airy', 'gaussian'):
        func = { 'airy': airy_disk, 'gaussian' : gaussian }[kind]
        fwhm = PacsBase.PSF_FWHM[band]
        size = int(np.round(10 * fwhm / resolution)) // 2 * 2 + 1
        return func((size,size), fwhm=fwhm, resolution=resolution)

    raise NotImplementedError()


#-------------------------------------------------------------------------------


def pacs_compute_delay(obs, tod, model, invntt=None, tol_delay=1.e-3,
                       tol_mapper=1.e-3, hyper=1., brack=(-60.,-30.,0.),
                       full_output=False):
    
    """
    Compute temporal offset between the PACS counter and the spacecraft clock.

    Parameters
    ----------
    obs : PacsObservation
    tod : Tod
    model : acquisition model
    invntt : N^-1 operator
    tol_delay : float
        Relative error in delay acceptable for convergence.
    tol_mapper : float
        Mapper tolerance.
    hyper : float
        Mapper smoothness hyperparameter.
    full_output : boolean
        if true, returns additional information.

    Returns
    -------

    delay : float
        Temporal offset between the PACS counter and the spacecraft clock.
    criterion : float
        Criterion of the map constructed with the returned delay.
    delays : list of floats
        Iterated delays.
    criteria : list of floats
        Criteria associated with the iterated delays.
    map : Map
        Map constructed with the returned delay.
    """
    global map_rls

    def func_model(obs, tod, model):

        def substitute_projection(model):
            for i, c in enumerate(model.operands):
                if isinstance(c, ProjectionOperator):
                    model.operands[i] = ProjectionOperator(obs, method=c.method,
                       header=c.header, npixels_per_sample=c.npixels_per_sample)
                    return True
                elif hasattr(c, 'operands'):
                    if substitute_projection(c):
                        return True
            return False

        if isinstance(model, ProjectionOperator):
            return ProjectionOperator(obs, method=model.method, header= \
                model.header, npixels_per_sample=model.npixels_per_sample)
        if hasattr(model, 'operands'):
            if substitute_projection(model):
                return model

        raise TypeError('There is no projection in this acquisition model.')
    
    delays = []
    criteria = []

    def criteria_delay(params, obs, tod, model, invntt, tol, hyper):
        global map_rls
        delay = params if np.rank(params) == 0 else params[0]
        obs.slice[0].delay = delay
        model = func_model(obs, tod, model)
        print('\nDelay: ' + str(delay) + 'ms...')
        map_rls = mapper_rls(tod, model, invntt=invntt, tol=tol, hyper=hyper)
        criterion = map_rls.header['criter']
        delays.append(delay)
        criteria.append(criterion)
        gc.collect()
        return criterion

    if obs.slice.size > 1:
        raise ValueError('Temporal delay is determined within a single obsid.')

    if weight is None:
        weight = IdentityOperator()

    method = scipy.optimize.brent
    result = method(criteria_delay, (obs, tod, model, weight, tol_mapper,
                    hyper), brack=brack, tol=tol_delay)

    map_rls.header.update('delay', result)

    if full_output:
        return result, map_rls.header['criter'], delays, criteria, map_rls

    return result


#-------------------------------------------------------------------------------


def _get_workload(band, slices, detector_masked, policy, transparent, comm):

    detector_removed = detector_masked.copy()
    if policy != 'remove':
        detector_removed[:] = False

    # mask the non-transmitting detectors in transparent mode
    if transparent:
        if band == 'red':
            detector_removed[0:8,0:8] = True
            detector_removed[0:8,16:] = True
            detector_removed[8:,:]    = True
        else:
            detector_removed[0:16,0:16] = True
            detector_removed[0:16,32:]  = True
            detector_removed[16:,:]     = True

    # distribute the workload over the processors
    return distribute_observation(detector_removed, slices, comm=comm)


#-------------------------------------------------------------------------------


def _files2tmf(filename):
    nfilenames = len(filename)
    length = max(len(f) for f in filename)
    filename_ = ''
    for f in filename:
        filename_ += f + (length-len(f))*' '
    return filename_, nfilenames


#-------------------------------------------------------------------------------


def _hcss_fits_keyword(header, keyword, *args):
    """ Return FITS keyword, according to HCSS handling. """
    if len(args) > 1:
        raise ValueError('Invalid number of arguments.')
    for k in header.keys():
        if not k.startswith('KEY.'):
            continue
        if header[k] == keyword:
            return header[k[4:]]
    if len(args) == 1:
        return args[0]
    raise KeyError("Keyword '{0}' not found.".format(keyword))


#-------------------------------------------------------------------------------


def _write_mask(obs, mask, filename):

    header = create_fitsheader(naxis=(), extname='Mask')
    header.update('DSETS___', 1)
    header.update('DS_0', 4)
    pyfits.append(filename, None, header)

    mask = mask.view()
    shape = obs.instrument.detector.shape
    mask.shape = (obs.instrument.detector.size, -1)
    bitmask = tmf.pacs_bitmask(mask.T.view(np.int8))
    bitmask = bitmask.T
    bitmask.shape = (shape[0], shape[1],-1)
    header = create_fitsheader(fromdata=bitmask, extname='Tamasis')
    header.update('INFO____', 'Tamasis mask')
    pyfits.append(filename, bitmask, header)


#-------------------------------------------------------------------------------


def _write_status(obs, filename, fitskw=None):

    s = obs.slice[0]

    if any(obs.slice.compression_factor != obs.slice[0].compression_factor):
        raise ValueError('Unable to save into a single file. The observations d'
                         'o not have the same compression factor.')
    compression_factor = s.compression_factor

    if any(obs.slice.mode != obs.slice[0].mode):
        raise ValueError('Unable to save into a single file. The observations d'
                         'o not have the same observing mode.')
    mode = s.mode
    channel = 'Red' if obs.instrument.band == 'red' else 'Blue'
    band_type = {'blue':'BS', 'green':'BS', 'red':'RS'}[obs.instrument.band]

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

    if obs.pointing.size != np.sum(obs.slice.nsamples_all):
        raise ValueError('The pointing and slice attribute are incompatible. Th'
                         'is should not happen.')

    fits = pyfits.HDUList()

    # Primary header
    cc = pyfits.create_card
    
    header = pyfits.Header([
        cc('simple', True), 
        cc('BITPIX', 32), 
        cc('NAXIS', 0), 
        cc('EXTEND', True, 'May contain datasets'), 
        cc('TYPE', 'HPPAVG'+band_type, 'Product Type Identification'), 
        cc('CREATOR', 'TAMASIS v' + var.__version__, 'Generator of this file'), 
        cc('INSTRUME', 'PACS', 'Instrument attached to this file'), 
        cc('TELESCOP', 'Herschel Space Observatory', 'Name of telescope'),
        cc('OBS_MODE', 'Scan map', 'Observation mode name'),
        cc('RA', s.ra),
        cc('DEC', s.dec),
        cc('EQUINOX', 2000., 'Equinox of celestial coordinate system'),
        cc('RADESYS', 'ICRS', 'Coordinate reference frame for the RA and DEC'),
        cc('CUSMODE', cusmode[mode], 'CUS observation mode'),
        cc('META_0', obs.instrument.detector.shape[0]),
        cc('META_1', obs.instrument.detector.shape[1]), 
        cc('META_2', channel + ' Photometer'),
        cc('META_3', ('Floating Average  : ' + str(compression_factor)) \
           if compression_factor > 1 else 'None'), 
        cc('META_4', comp_mode),
        cc('META_5', s.scan_angle),
        cc('META_6', s.scan_length),
        cc('META_7', s.scan_nlegs),
        cc('META_8', 'blue2' if obs.instrument.band == 'green' else 'blue1'),
        cc('HIERARCH key.META_0', 'detRow'), 
        cc('HIERARCH key.META_1', 'detCol'), 
        cc('HIERARCH key.META_2', 'camName'), 
        cc('HIERARCH key.META_3', 'algorithm'), 
        cc('HIERARCH key.META_4', 'compMode'),
        cc('HIERARCH key.META_5', 'mapScanAngle'),
        cc('HIERARCH key.META_6', 'mapScanLegLength'),
        cc('HIERARCH key.META_7', 'mapScanNumLegs'),
        cc('HIERARCH key.META_8', 'blue'),
    ])

    if fitskw is not None:
        for k, v in fitskw.items():
            if len(k) > 8:
                k = 'HIERARCH ' + k
            if isinstance(v, tuple):
                if len(v) != 2:
                    raise ValueError('A tuple input for fitskw must have two el'
                                     'ements (value, comment).')
                v, c = v
            else:
                c = None
            v = v if isinstance(v, (str, int, long, float, complex, bool,
                  np.floating, np.integer, np.complexfloating)) else repr(v)
            header.update(k, v, c)

    hdu = pyfits.PrimaryHDU(None, header)
    fits.append(hdu)
    
    # Status
    v = ~obs.pointing.removed
    status = obs.status[v]
    status.RaArray = obs.pointing.ra[v]
    status.DecArray = obs.pointing.dec[v]
    status.PaArray = obs.pointing.pa[v]
    status.CHOPFPUANGLE = obs.pointing.chop[v]
    hdu = pyfits.BinTableHDU(status, None, name='STATUS')
    fits.append(hdu)
    fits.writeto(filename, clobber=True)


#-------------------------------------------------------------------------------


def _scan_speed(pointing):
    if pointing.size == 1:
        return np.nan
    velocity = pointing.velocity
    inscan = (pointing.info == pointing.INSCAN) & ~pointing.removed
    return np.round(np.median(velocity[inscan]),1)


#-------------------------------------------------------------------------------
# processing steps

def step_scanline_masking(obs, n_repetition=0, n_scanline=None):
    """
    Mask everything up to the n_scanline scan line of the n_repetition
    repetition. If n_scanline is None, mask the whole repetition.

    Arguments
    ----------

    obs: PacsObservation
        The considered PACS observation instance.

    n_repetition: int
        Mask everything up to n_repetition.

    n_scanline: int or None (default: None)
        Mask everything up to n_scanline of n_repetition.
        If None mask n_repetition entirely.

    Returns
    -------

    Returns nothing, obs is altered inplace.

    Note
    ----

    scan line index starts at 1 !
    """
    if n_scanline == None:
        for i in xrange(n_repetition):
            obs.pointing.masked[obs.status.Repetition == i] = True
    else:
        for i in xrange(n_repetition - 1):
            obs.pointing.masked[obs.status.Repetition == i] = True
        for j in xrange(1, n_scanline):
            obs.pointing.masked[(obs.status.Repetition == n_repetition) *
                                (obs.status.ScanLineNumber == j)] = True

def step_deglitching(obs, tod, length=100, nsigma=25., method="mad"):
    """
    Include all the deglitching steps : 

       - define a sharp projector with downsampling and with
         default resolution.

       - perform tod filtering

       - second level deglitching with mad or std

    Arguments
    ---------
    obs : PacsObservation instance
    tod : Time Ordered Data
    length : int (default: 100)
        Median filtering window size.
    nsigma : float (default: 25.)
        Median filtering threshold
    method : "mad" or "std" (default: "mad")
        Filtering method (median absolute deviation or standard deviation).

    Returns
    -------
    Returns nothing. tod mask is updated in-place.
    """
    import tamasis as tm

    if method == "mad":
        method = tm.deglitch_l2mad
    elif method == "std":
        method = tm.deglitch_l2std
    else:
        raise ValueError("Unrecognized deglitching method.")

    # special deglitching projector
    proj_glitch = tm.ProjectionOperator(obs,
                                        method='sharp',
                                        downsampling=True,
                                        npixels_per_sample=6)
    # filter tod with narrow window
    tod_glitch = tm.filter_median(tod, length=length)
    # actual degltiching according to selected method (mad or std)
    tod.mask = method(tod_glitch, proj_glitch, nsigma=nsigma)

