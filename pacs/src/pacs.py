import glob
import kapteyn
import numpy
import os
import pyfits
import re
import tamasisfortran as tmf
import tempfile
import utils

from acquisitionmodels import AcquisitionModel, CompressionAverage, Masking, Projection, ValidationError
from config import __version__, tamasis_dir
from datatypes import *
from mappers import mapper_naive
from matplotlib import pyplot
from mpi4py import MPI
from processing import deglitch_l2mad, filter_median
from unit import Quantity

__all__ = [ 'PacsObservation', 'PacsSimulation', 'pacs_plot_scan', 'pacs_preprocess', 'pacs_status' ]


class _Pacs(object):
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
        ncorrelations, status = tmf.pacs_read_filter_calibration_ncorrelations(tamasis_dir, self.band)
        if status != 0: raise RuntimeError()

        data, status = tmf.pacs_read_filter_calibration(tamasis_dir, self.band, ncorrelations, self.ndetectors, numpy.asfortranarray(self.detector_mask))
        if status != 0: raise RuntimeError()

        return data.T

   
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

        band, transparent_mode, status = tmf.pacs_info_band(filename_, nfilenames)
        if status != 0: raise RuntimeError()
        band = band.strip()

        nrows, ncolumns = (16,32) if band == 'red' else (32,64)

        # get the detector mask, before distributing to the processors
        if detector_mask is not None:
            if detector_mask.shape != (nrows, ncolumns):
                raise ValueError('Invalid shape of the input detector mask: ' + str(detector_mask.shape) + ' for the ' + band + ' band.')
            detector_mask = numpy.array(detector_mask, dtype='int8', copy=False)

        else:
            detector_mask, status = tmf.pacs_info_bad_detector_mask(tamasis_dir, band, transparent_mode, reject_bad_line, nrows, ncolumns)
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
        frame_policy = utils.MaskPolicy('inscan,turnaround,other,invalid'.split(','), (frame_policy_inscan, frame_policy_turnaround, frame_policy_other, frame_policy_invalid), 'Frame Policy')

        # retrieve information from the observations
        nthreads = tmf.info_nthreads()
        print 'Info: ' + MPI.Get_processor_name() + ' (' + str(nthreads) + ' core' + ('s' if nthreads > 1 else '') + ' handling ' + str(ndetectors) + ' detector' + ('s' if ndetectors > 1 else '') + ')'
        compression_factor, nsamples, unit, responsivity, detector_area, dflat, oflat, status = tmf.pacs_info(tamasis_dir, filename_, nfilenames, transparent_mode, fine_sampling_factor, numpy.array(frame_policy), numpy.asfortranarray(detector_mask))
        if status != 0: raise RuntimeError()

        self.filename = filename
        self.band = band
        self.nrows = nrows
        self.ncolumns = ncolumns
        self.default_npixels_per_sample = 11 if band == 'red' else 6
        self.default_resolution = 6.4 if band == 'red' else 3.2
        self.nobservations = nfilenames
        self.nsamples = tuple(nsamples)
        self.nfinesamples = tuple(nsamples * compression_factor * fine_sampling_factor)
        self.ndetectors = ndetectors
        self.detector_mask = detector_mask
        self.reject_bad_line = reject_bad_line
        self.frame_policy = frame_policy
        self.fine_sampling_factor = fine_sampling_factor
        self.mode = '' #XXX FIX ME!
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
   

#-------------------------------------------------------------------------------


class PacsSimulation(_Pacs):
    def __init__(self, band, mode, ra0, dec0, pa0=0., scan_angle=0., scan_length=30., scan_nlegs=3, scan_step=20., scan_speed=10., fine_sampling_factor=1):
        band = band.lower()
        if band not in ('blue', 'green', 'red'):
            raise ValueError("Band is not 'blue', 'green', nor 'red'.")

        mode = mode.lower()
        if mode not in ('prime', 'parallel', 'transparent'):
            raise ValueError("Observing mode is not 'prime', 'parallel', nor 'transparent'.")
        
        compression_factor = 8 if mode == 'parallel' and band != 'red' else 1 if mode == 'transparent' else 4
        self._file = tempfile.NamedTemporaryFile('w', suffix='.fits')
        self.pointing = _generate_pointing(ra0, dec0, pa0, scan_angle, scan_length, scan_nlegs, scan_step, scan_speed, compression_factor)
        self.ra0 = ra0
        self.dec0 = dec0
        self.pa0 = pa0
        self.scan_angle = scan_angle
        self.scan_length = scan_length
        self.scan_nlegs = scan_nlegs
        self.scan_step = scan_step
        self.scan_speed = scan_speed

        self.filename = self._file.name
        self.band = band
        self.nrows, self.ncolumns = (16,32) if band == 'red' else (32,64)
        self.default_npixels_per_sample = 11 if band == 'red' else 6
        self.default_resolution = 6.4 if band == 'red' else 3.2
        self.nobservations = 1
        self.nsamples = (self.pointing.size,)
        self.nfinesamples = (self.pointing.size * compression_factor * fine_sampling_factor,)
        self.ndetectors = self.nrows * self.ncolumns
        self.detector_mask = numpy.zeros((self.nrows, self.ncolumns), dtype='int8')
        self.reject_bad_line = False
        self.frame_policy = frame_policy = utils.MaskPolicy('inscan,turnaround,other,invalid'.split(','), 4*('keep',), 'Frame Policy')
        self.fine_sampling_factor = fine_sampling_factor
        self.mode = mode
        self.transparent_mode = mode == 'transparent'
        self.compression_factor = compression_factor
#        self.unit = unit.strip()
#        if self.unit.find('/') == -1:
#            self.unit += ' / detector'
#        self.responsivity = Quantity(responsivity, 'V/Jy')
#        self.detector_area = Map(detector_area, unit='arcsec^2/detector')
#        self.flatfield = {
#            'total'   : Map(dflat*oflat),
#            'detector': Map(numpy.ascontiguousarray(dflat)),
#            'optical' : Map(numpy.ascontiguousarray(oflat))
#            }

        _write_status(self, self.filename)

    def plot(self, map=None, title=None, num=None, figsize=None, dpi=None):
        if map is None:
            crval = [utils.mean_degrees(self.pointing.ra), 
                     numpy.mean(self.pointing.dec)]
            ra_min,  ra_max  = utils.minmax_degrees(self.pointing.ra)
            dec_min, dec_max = numpy.min(self.pointing.dec), numpy.max(self.pointing.dec)
            cdelt = numpy.max((ra_max-ra_min)/1000., (dec_max-dec_min)/1000.)
            header = utils.create_fitsheader(None, naxis=[1,1], cdelt=cdelt, crval=crval, crpix=[1,1])
            proj = kapteyn.wcs.Projection(header)
            coords = numpy.array([self.pointing.ra, self.pointing.dec]).T
            xy = proj.topixel(coords)
            xmin = int(numpy.round(numpy.min(xy[:,0])))
            xmax = int(numpy.round(numpy.max(xy[:,0])))
            ymin = int(numpy.round(numpy.min(xy[:,1])))
            ymax = int(numpy.round(numpy.max(xy[:,1])))
            xmargin = int(numpy.ceil((xmax - xmin + 1) / 10.))
            ymargin = int(numpy.ceil((ymax - ymin + 1) / 10.))
            xmin -= xmargin
            xmax += xmargin
            ymin -= ymargin
            ymax += ymargin
            naxis = (xmax - xmin + 1, ymax - ymin + 1)
            crpix = (-xmin+2, -ymin+2)
            header = utils.create_fitsheader(None, naxis=naxis, cdelt=cdelt, crval=crval, crpix=crpix)
            fitsobj = kapteyn.maputils.FITSimage(externalheader=header)
            fig = pyplot.figure(num=num, figsize=figsize, dpi=dpi)
            frame = fig.add_axes([0.1,0.1,0.9,0.9])
            if title is not None:
                frame.set_title(title)
            image = fitsobj.Annotatedimage(frame, blankcolor='w')
            grat = image.Graticule()
            image.plot()
            image.interact_toolbarinfo()
            image.interact_writepos()
            pyplot.show()
        else:
            image = map.imshow(title=title, num=num, figsize=figsize, dpi=dpi)
        x, y = image.topixel(self.pointing.ra, self.pointing.dec)
        pyplot.plot(x, y, color='magenta', linewidth=2)
        pyplot.plot(x[0], y[0], 'o', color='magenta')

    def save(self, filename, tod):
        _write_status(self, filename)
        header = create_fitsheader(tod, extname='Signal')
        pyfits.append(filename, tod, header)
        if tod.mask is not None:
            mask = numpy.abs(self.mask).view('uint8')
            header = create_fitsheader(mask, extname='Mask')
            pyfits.append(filename, mask)
        
   
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
        pyplot.plot(status['ra'], status['dec'])


#-------------------------------------------------------------------------------


class pacs_status(object):
    def __init__(self, filename):
        hdu = pyfits.open(filename)['STATUS']
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


def _generate_pointing(ra0, dec0, pa0, scan_angle, scan_length=30., scan_nlegs=3, scan_step=20., scan_speed=10., compression_factor=4):
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
    flags         = numpy.zeros(nsamples, dtype=int)
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
        flag = 255
   
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
            flag  = 0
   
        # constant speed
        if working_time >=  extra_time1 and working_time < extra_time1+line_time:
            delta = signe*(-scan_length/2+ (working_time-extra_time1)*scan_speed)
            flag  = 1
 
        # Deceleration at then end of the scanline to stop
        if working_time >= extra_time1+line_time and working_time < extra_time1+line_time+extra_time1:
            dt = working_time - extra_time1 - line_time
            delta = signe * (scan_length/2 + scan_speed*dt - 0.5 * gamma*dt*dt)
            flag  = 2
  
        # Acceleration to go toward the next scan line
        if working_time >= 2*extra_time1+line_time and working_time < 2*extra_time1+line_time+extra_time2:
            dt = working_time-2*extra_time1-line_time
            alpha = alpha0 + 0.5*gamma*dt*dt
            flag  = 3
   
        # Deceleration to stop at the next scan line
        if working_time >= 2*extra_time1+line_time+extra_time2 and working_time < full_line_time:
            dt = working_time-2*extra_time1-line_time-extra_time2
            speed = gamma*extra_time2
            alpha = (alpha0+scan_step/2.) + speed*dt - 0.5*gamma*dt*dt
            flag  = 4

        time[i] = i / sampling_frequency
        flags[i] = flag
        latitude[i] = delta
        longitude[i] = alpha
        line_counters[i] = line_counter
        working_time = working_time + sampling_period

    # Convert the longitude and latitude *expressed in degrees) to ra and dec
    ra, dec = _change_coord(ra0, dec0, scan_angle, -longitude/3600., latitude/3600.)

    p = numpy.recarray(nsamples, dtype=utils.pointing)
    p.time = time
    p.ra = ra
    p.dec = dec
    p.pa = pa0
    p[flags == 1].flag = 0xcd2
    p[flags != 1].flag = 0x4000

    return p


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


def _write_status(obs, filename):

    band = 'BS' if obs.band == 'blue' else 'BL' if obs.band == 'green' else 'R '

    if obs.mode == 'prime':
        observing_mode = 'Photometry Default Mode'
    elif mode == 'parallel':
        observing_mode = 'Photometry Double Compression Mode'
    else:
        observing_mode = 'Photometry Lossless Compression Mode'

    fits = pyfits.HDUList()

    # Primary header
    cc = pyfits.createCard
    header = pyfits.Header([
            cc('simple', True), 
            cc('BITPIX', 32), 
            cc('NAXIS', 0), 
            cc('EXTEND', True), 
            cc('TYPE', 'HPPAVG'+band.strip()), 
            cc('CREATOR', 'TAMASIS v' + __version__), 
            cc('INSTRUME', 'PACS    '), 
            cc('SOURCE', 'largeScan'), 
            cc('RA', obs.ra0),
            cc('DEC', obs.dec0),
            cc('META_0', obs.nrows), 
            cc('META_1', obs.ncolumns), 
            cc('META_2', obs.band.title()+' Photometer'), 
            cc('META_3', ('Floating Average  : '+str(obs.compression_factor)) if obs.compression_factor > 1 else 'None'), 
            cc('META_4', observing_mode),
            cc('META_5', obs.scan_angle),
            cc('META_6', obs.scan_length),
            cc('META_7', obs.scan_nlegs),
            cc('META_8', obs.scan_speed),
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
    
    # write status
    table = numpy.recarray(obs.nsamples[0], dtype=[('BBID', numpy.int64), ('FINETIME', numpy.int64), ('BAND', 'S2'), ('CHOPFPUANGLE', numpy.float64), ('RaArray', numpy.float64), ('DecArray', numpy.float64), ('PaArray', numpy.float64)])
    table.BAND = band
    table.FINETIME = numpy.round(obs.pointing.time*1000000.)
    table.RaArray = obs.pointing.ra
    table.DecArray = obs.pointing.dec
    table.PaArray = obs.pointing.pa
    table.CHOPFPUANGLE = 0.
    table.BBID = obs.pointing.flag

    status = pyfits.BinTableHDU(table, None, 'STATUS')
    fits.append(status)
    fits.writeto(filename)


