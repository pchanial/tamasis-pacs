from __future__ import division

import gc
import numpy as np

from kapteyn import wcs
from matplotlib import pyplot
from pyoperators.utils import strenum, strshape
from pyoperators.utils.mpi import MPI
from . import tamasisfortran as tmf
from .acquisitionmodels import PointingMatrix
from .datatypes import Map
from .mpiutils import gather_fitsheader_if_needed
from .wcsutils import barycenter_lonlat, combine_fitsheader, create_fitsheader

__all__ = ['Instrument']


class Instrument(object):
    """
    Class storing information about the instrument.

    Attributes
    ----------
    name
    comm
    detector.removed
    detector.masked
    detector.center
    detector.corner
    """
    def __init__(self, name, shape, removed=None, masked=None,
                 detector_center=None, detector_corner=None,
                 default_resolution=None, dtype=None, comm=MPI.COMM_WORLD):

        self.name = str(name)
        shape = tuple(shape)
        self.default_resolution = default_resolution
        self.comm = comm

        dtype_default = [('masked', np.bool8), ('removed', np.bool8)]

        if removed is not None:
            if removed.shape != shape:
                raise ValueError('The input specifying the removed detectors ha'
                                 's an incompatible shape.')
        else:
            removed = False

        if masked is not None:
            if masked.shape != shape:
                raise ValueError('The input specifying the masked detectors has'
                                 ' an incompatible shape.')
        else:
            masked = False

        if detector_center is not None:
            if detector_center.shape != shape:
                raise ValueError('The detector centers have an incompatible sha'
                                 'pe.')
            dtype_default += [('center', float, 2)]

        if detector_corner is not None:
            if detector_corner.ndim != len(shape) + 2 or \
               detector_corner.shape[0:2] != shape:
                raise ValueError('The detector corners have an incompatible sha'
                                 'pe.')
            dtype_default += [('center', float, 2), ('corner', float, (4,2))]
            detector_center = np.mean(detector_corner, axis=-2)

        if dtype is None:
            dtype = dtype_default

        self.detector = Map.zeros(shape, dtype=dtype, origin='upper')
        
        self.detector.masked = masked
        self.detector.removed = removed
        if detector_center is not None:
            self.detector.center = detector_center
        if detector_corner is not None:
            self.detector.corner = detector_corner

    def get_ndetectors(self, masked=False):
        """
        Return the number of valid detectors.

        Parameters
        ----------
        masked : boolean
              If False, the valid detectors are those that are not removed.
              Otherwise, they are those that not removed and not masked.
        """
        if masked:
            return int(np.sum(~self.detector.removed & ~self.detector.masked))
        return int(np.sum(~self.detector.removed))

    def get_valid_detectors(self, masked=False):
        """
        Return valid detectors subscripts.

        Parameters
        ----------
        masked : boolean
              If False, the valid detectors are those that are not removed.
              Otherwise, they are those that not removed and not masked.
        """
        mask = ~self.detector.removed
        if masked:
            mask &= ~self.detector.masked
        if np.sum(mask) == self.detector.size:
            return Ellipsis
        return mask

    def pack(self, input, masked=False):
        """
        Convert an ndarray which only includes the valid detectors into 
        another ndarray which contains all the detectors under the control
        of the detector mask.

        Parameters
        ----------
        input : ndarray 
              Array to be packed, whose first dimensions are equal to those
              of the detector attribute.
        masked : boolean
              If False, the valid detectors are those that are not removed.
              Otherwise, they are those that not removed and not masked.

        Returns
        -------
        output : ndarray
              Packed array, whose first dimension is the number of valid
              detectors.

        See Also
        --------
        unpack : inverse method.

        Notes
        -----
        This method does not necessarily make a copy.
        """
        if isinstance(input, dict):
            output = input.copy()
            for k, v in output.items():
                try:
                    output[k] = self.pack(v, masked=masked)
                except (ValueError, TypeError):
                    output[k] = v
            return output
        if not isinstance(input, np.ndarray):
            raise TypeError('The input is not an ndarray.')
        if input.ndim < self.detector.ndim or \
           input.shape[:self.detector.ndim] != self.detector.shape:
            raise ValueError("The shape of the argument '{0}' is incompatible w"
                "ith that of the detectors '{1}'.".format(strshape(input.shape),
                strshape(self.detector.shape)))
        index = self.get_valid_detectors(masked=masked)
        if index is Ellipsis:
            new_shape = (-1,) + input.shape[self.detector.ndim:]
            output = input.reshape(new_shape)
        else:
            output = input[index]
        if type(input) != np.ndarray:
            output = output.view(type(input))
            for k, v in input.__dict__.items():
                try:
                    output.__dict__[k] = self.pack(v, masked=masked)
                except (ValueError, TypeError):
                    output.__dict__[k] = v
        return output

    def unpack(self, input, masked=False):
        """
        Convert an ndarray which only includes the valid detectors into 
        another ndarray which contains all the detectors under the control
        of the detector mask.

        Parameters
        ----------
        input : ndarray 
              Array to be unpacked, whose first dimension is equal to the
              number of valid detectors.
        masked : boolean
              If False, the valid detectors are those that are not removed.
              Otherwise, they are those that not removed and not masked.

        Returns
        -------
        output : ndarray
              Unpacked array, whose first dimensions are those of the detector
              attribute.

        See Also
        --------
        pack : inverse method.

        Notes
        -----
        This method does not necessarily make a copy.
        """
        if isinstance(input, dict):
            output = input.copy()
            for k, v in output.items():
                try:
                    output[k] = self.unpack(v, masked=masked)
                except (ValueError, TypeError):
                    output[k] = v
            return output
        if not isinstance(input, np.ndarray):
            raise TypeError('The input is not an ndarray.')

        n = self.get_ndetectors(masked=masked)
        if input.ndim == 0 or n != input.shape[0]:
            raise ValueError("The shape of the argument '{0}' is incompatible w"
                "ith the number of valid detectors '{1}'.".format(strshape(
                input.shape),n))
        index = self.get_valid_detectors(masked=masked)
        new_shape = self.detector.shape + input.shape[1:]
        if index is Ellipsis:
            return input.reshape(new_shape)
        output = np.zeros(new_shape, dtype=input.dtype)
        output[index] = input
        if type(input) != np.ndarray:
            output = output.view(type(input))
            for k, v in input.__dict__.items():
                try:
                    output.__dict__[k] = self.unpack(v, masked=masked)
                except (ValueError, TypeError):
                    output.__dict__[k] = v
        return output

    def get_map_header(self, pointing, resolution=None):
        """
        Return the FITS header of the smallest map that encompasses
        a set of pointings, by taking into account the instrument geometry.

        Parameters
        ----------
        pointing : array of flexible type
            Pointing directions.
        resolution : float
            Sky pixel increment, in arc seconds. Default is .default_resolution.

        Returns
        -------
        header : pyfits.Header
            The resulting FITS header.
        """
        if resolution is None:
            resolution = self.default_resolution
        if 'corner' in self.detector.dtype.names:
            coords = self.detector.corner
        else:
            coords = self.detector.center
        coords = self.pack(coords, masked=True)

        mask = ~pointing['removed'] & ~pointing['masked']
        if not np.any(mask):
            raise ValueError('The FITS header cannot be inferred: there is no v'
                             'alid pointing.')
        pointing = pointing[mask]

        # get a dummy header, with correct cd and crval
        ra0, dec0 = barycenter_lonlat(pointing['ra'], pointing['dec'])
        header = create_fitsheader((1,1), cdelt=resolution/3600,
                                   crval=(ra0,dec0), crpix=(1,1))

        # compute the coordinate boundaries according to the header's astrometry
        xmin, ymin, xmax, ymax = self.instrument2xy_minmax(
            coords, pointing, str(header).replace('\n',''))
        ixmin = int(np.round(xmin))
        ixmax = int(np.round(xmax))
        iymin = int(np.round(ymin))
        iymax = int(np.round(ymax))
        nx = ixmax-ixmin+1
        ny = iymax-iymin+1

        # move the reference pixel (not the reference value!)
        header = create_fitsheader((nx,ny), cdelt=resolution/3600,
                                   crval=(ra0,dec0), crpix=(-ixmin+2,-iymin+2))

        # gather and combine the FITS headers
        headers = self.comm.allgather(header)
        return combine_fitsheader(headers)

    def get_pointing_matrix(self, pointing, header, npixels_per_sample=0,
                            method=None, downsampling=False,
                            units=('/detector', '/pixel'),
                            derived_units=(None,None), comm=MPI.COMM_WORLD,
                            **keywords):
        """
        Return the pointing matrix for a given set of pointings.

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
            If True, return a pointing matrix downsampled by the instrument
            fine sampling factor. Otherwise return a pointing matrix sampled
            at the instrument fine sampling factor.
        comm : mpi4py.MPI.Comm
            MPI communicator, only used if the input FITS header is local.

        """
        if method is None:
            if 'corner' in self.detector.dtype.names:
                method = 'sharp'
            else:
                method = 'nearest'
        method = method.lower()
        choices = ('nearest', 'sharp')
        if method not in choices:
            raise ValueError("Invalid method '" + method + "'. Expected values "
                "are " + strenum(choices) + '.')

        if not downsampling and self.fine_sampling_factor > 1:
            raise NotImplementedError('Oversampling is not implemented.')

        header = gather_fitsheader_if_needed(header, comm=comm)
        shape_input = tuple(header['NAXIS' + str(i+1)]
                            for i in range(header['NAXIS']))[::-1]
        if np.product(shape_input) > np.iinfo(np.int32).max:
            raise RuntimeError('The map is too large: pixel indices cannot be s'
                'stored using 32 bits: {0}>{1}'.format(np.product(shape_input),
                np.iinfo(np.int32).max))

        mask = ~pointing['removed']
        pointing = pointing[mask]
        nvalids = pointing.size

        if method == 'nearest':
            npixels_per_sample = 1

        # Allocate memory for the pointing matrix
        ndetectors = self.get_ndetectors()
        shape = (ndetectors, nvalids, npixels_per_sample)
        info = {'header':header,
                'method':method,
                'units':units,
                'derived_units':derived_units,
                'outside':False,
                'npixels_per_sample_min':0}
        try:
            pmatrix = PointingMatrix.empty(shape, shape_input, info=info,
                                           verbose=False)
        except MemoryError:
            gc.collect()
            pmatrix = PointingMatrix.empty(shape, shape_input, info=info,
                                           verbose=False)

        # compute the pointing matrix
        if method == 'sharp':
            coords = self.pack(self.detector.corner)
            new_npixels_per_sample, outside = self. \
                instrument2pmatrix_sharp_edges(coords, pointing, header,
                                               pmatrix, npixels_per_sample)
        elif method == 'nearest':
            coords = self.pack(self.detector.center)
            new_npixels_per_sample = 1
            raise NotImplementedError()
        else:
            raise NotImplementedError()

        info['outside'] = bool(outside)
        info['npixels_per_sample_min'] = new_npixels_per_sample

        if new_npixels_per_sample <= npixels_per_sample:
            return pmatrix

        # if the actual number of pixels per sample is greater than
        # the specified one, redo the computation of the pointing matrix
        del pmatrix
        return self.get_pointing_matrix(pointing, header,
            new_npixels_per_sample, method, downsampling, comm, **keywords)

    @staticmethod
    def instrument2ad(coords, pointing):
        """
        Convert coordinates in the instrument frame into celestial coordinates,
        assuming a pointing direction and a position angle.
        
        The instrument frame is the frame used for the 'center' or 'corner'
        coordinates, which are fields of the 'detector' attribute.

        Parameters
        ----------
        coords : float array (last dimension is 2)
            Coordinates in the instrument frame in arc seconds.
        pointing : array with flexible dtype including 'ra', 'dec' and 'pa'
            Direction corresponding to (0,0) in the local frame, in degrees.

        Returns
        -------
        coords_converted : float array (last dimension is 2, for Ra and Dec)
            Converted coordinates, in degrees.

        Notes
        -----
        The routine is not accurate at the poles.
        """
        coords = np.array(coords, float, order='c', copy=False)
        shape = coords.shape
        coords = coords.reshape((-1,2))
        new_shape = (pointing.size,) + coords.shape
        result = np.empty(new_shape, float)
        for r, ra, dec, pa in zip(result, pointing['ra'].flat,
                                  pointing['dec'].flat, pointing['pa'].flat):
            tmf.pointing.instrument2ad(coords.T, r.T, ra, dec, pa)
        result = result.reshape(pointing.shape + shape)
        for dim in range(pointing.ndim):
            result = np.rollaxis(result, 0, -1)
        return result

    def instrument2xy(self, coords, pointing, header):
        """
        Convert coordinates in the instrument frame into sky pixel coordinates,
        assuming a pointing direction and a position angle.
        """
        coords = np.array(coords, float, order='c', copy=False)
        proj = wcs.Projection(header)
        return proj.topixel(self.instrument2ad(coords, pointing))
        

    @staticmethod
    def instrument2xy_minmax(coords, pointing, header):
        """
        Return the minimum and maximum sky pixel coordinate values given a set
        of coordinates specified in the instrument frame.
        """
        coords = np.array(coords, float, order='c', copy=False)
        xmin, ymin, xmax, ymax, status = tmf.pointing.instrument2xy_minmax(
            coords.reshape((-1,2)).T, pointing['ra'].ravel(),
            pointing['dec'].ravel(), pointing['pa'].ravel(),
            str(header).replace('\n',''))
        if status != 0:
            raise RuntimeError()
        return xmin, ymin, xmax, ymax

    @staticmethod
    def instrument2pmatrix_sharp_edges(coords, pointing, header, pmatrix,
                                       npixels_per_sample):
        """
        Return the dense pointing matrix whose values are intersection between
        detectors and map pixels.
        """
        coords = coords.reshape((-1,2))
        ra = pointing['ra'].ravel()
        dec = pointing['dec'].ravel()
        pa = pointing['pa'].ravel()
        masked = pointing['masked'].view(np.int8).ravel()
        if pmatrix.size == 0:
            # f2py doesn't accept zero-sized opaque arguments
            pmatrix = np.empty(1, np.int64)
        else:
            pmatrix = pmatrix.ravel().view(np.int64)
        header = str(header).replace('\n','')

        new_npixels_per_sample, out, status = tmf.pointing. \
            instrument2pmatrix_sharp_edges(coords.T, ra, dec, pa, masked,
            header, pmatrix, npixels_per_sample)
        if status != 0: raise RuntimeError()

        return new_npixels_per_sample, out

    @classmethod
    def create_grid(cls, nrows, ncolumns, size, active_fraction=1.,
                    focal_distance=None, inversion=False):
        """
        Return the physical positions of the corners of square detectors in a
        matrix of shape (nrows, ncolumns).
        If the focal distance is provided, the detector positions are then
        expressed in arc seconds, provided that the focal distance units are
        the same as that of the size argument.

        Parameters
        ----------
        nrows : int
            Number of rows of detectors.
        ncolumns : int
            Number of columns of detectors.
        size : float
            Detector size, in the same units as the focal distance.
        active_fraction : float
            Fraction of the detector surface that transmits light.
        focal_distance : float
            The focal distance, in the same units as the detector size.
        inversion : boolean
            If set the true, the detector positions are inversed through
            the detector center.
        """
        corners = np.empty((nrows, ncolumns, 4, 2), float)
        i, j = np.ogrid[0:nrows,0:ncolumns]
        corners[:,:,0,0] = size * j
        corners[:,:,0,1] = size * i
        corners[:,:,1,0] = size * (j+1)
        corners[:,:,1,1] = size * i
        corners[:,:,2,0] = size * (j+1)
        corners[:,:,2,1] = size * (i+1)
        corners[:,:,3,0] = size * j
        corners[:,:,3,1] = size * (i+1)
        corners[...,0] -= np.mean(corners[...,0])
        corners[...,1] -= np.mean(corners[...,1])
        if active_fraction != 1:
            coef = np.sqrt(active_fraction)
            for i in range(nrows):
                for j in range(ncolumns):
                    u0, v0 = np.mean(corners[i,j,:,0]),np.mean(corners[i,j,:,1])
                    corners[i,j,:,0] = (corners[i,j,:,0] - u0) * coef + u0
                    corners[i,j,:,1] = (corners[i,j,:,1] - v0) * coef + v0

        corners[...,0], corners[...,1] = \
            np.rad2deg(corners[...,0] / focal_distance), \
            np.rad2deg(corners[...,1] / focal_distance)
        corners *= 3600

        if inversion:
            corners = -corners

        return corners

    @classmethod
    def plot_grid(cls, corners, update_axes=True):
        """
        Plot a detector grid.

        Parameters
        ----------
        corners : float ndarray [...,4,2]
            The detector corners the penultimate dimension is for the vertices,
            and the last one is for x and y, as defined in the current graphic
            axis instance.
        update_axes : boolean
            If true, the axes of the plot will updated to match the boundaries
            of the input detector corners.

        Example
        -------
        # overlay the detector grid on the observation pointings
        obs = MyObservation(...)
        annim = plot_scan(obs.pointing)
        xy = obs.instrument.instrument2xy(obs.instrument.detector.corner,
                                          obs.pointing[0], annim.hdr)
        plot_detector_grid(xy)

        """
        a = pyplot.gca()
        corners = corners.view(float).reshape(-1,4,2)
        for i in range(corners.shape[0]):
            a.add_patch(pyplot.Polygon(corners[i], closed=True, fill=False))
        if update_axes:
            xlim = a.get_xlim()
            ylim = a.get_ylim()
            pyplot.xlim((min(xlim[0], np.min(corners[...,0])),
                         max(xlim[1], np.max(corners[...,0]))))
            pyplot.ylim((min(ylim[0], np.min(corners[...,1])),
                         max(ylim[1], np.max(corners[...,1]))))
        pyplot.show()
