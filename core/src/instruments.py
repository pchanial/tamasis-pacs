from __future__ import division
import gc
import numpy as np
import tamasisfortran as tmf
from kapteyn import wcs
from mpi4py import MPI
from pyoperators.utils import strenum, strshape
from .datatypes import Map
from .wcsutils import barycenter_lonlat, combine_fitsheader, create_fitsheader

class Instrument(object):
    """
    Class storing information about the instrument.

    Attributes
    ----------
    name
    comm_tod
    detector.removed
    detector.masked
    detector.center
    detector.corner
    """
    def __init__(self, name, shape, removed=None, masked=None,
                 detector_center=None, detector_corner=None,
                 default_resolution=None, dtype=None, comm_tod=MPI.COMM_WORLD):

        self.name = str(name)
        shape = tuple(shape)
        self.default_resolution = default_resolution
        self.comm_tod = comm_tod

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
            for k,v in input.__dict__.items():
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
            for k,v in input.__dict__.items():
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
        headers = self.comm_tod.allgather(header)
        return combine_fitsheader(headers)
        
    def get_pointing_matrix(self, pointing, header, npixels_per_sample=0,
                            method=None):
        if method is None:
            if 'corner' in self.detector.dtype.names:
                method = 'sharp'
            else:
                method = 'nearest'
        method = method.lower()
        choices = ('nearest', 'sharp')
        if method not in choices:
            raise ValueError("Invalid method '" + method + \
                "'. Expected values are " + strenum(choices) + '.')

        mask = ~pointing['removed']
        pointing = pointing[mask]
        nvalids = pointing.size

        if method == 'sharp':
            coords = self.pack(self.detector.corner)
        else:
            coords = self.pack(self.detector.center)
            npixels_per_sample = 1

        # Allocate memory for the pointing matrix
        ndetectors = coords.shape[0]
        shape = (ndetectors, nvalids, npixels_per_sample)
        dtype = [('weight', 'f4'), ('pixel', 'i4')]
        if npixels_per_sample != 0:
            print('Info: Allocating ' + str(np.product(shape) / 2**17) + ' MiB '
                  'for the pointing matrix.')
        try:
            pmatrix = np.empty(shape, dtype).view(np.recarray)
        except MemoryError:
            gc.collect()
            pmatrix = np.empty(shape, dtype).view(np.recarray)

        # compute the pointing matrix
        if method == 'sharp':
            new_npixels_per_sample, out = self. \
                instrument2pmatrix_sharp_edges(coords, pointing, header,
                                               pmatrix, npixels_per_sample)
        else:
            new_npixels_per_sample = 1
            raise NotImplementedError()

        if ndetectors > 0 and nvalids > 0:
            if new_npixels_per_sample == 0:
                print('Warning:  All detectors fall outside the map.')
            elif out:
                print('Warning: Some detectors fall outside the map.')

        if new_npixels_per_sample != npixels_per_sample:
            print("Warning: For this observation, you can set the keyword 'npix"
                  "els_per_sample' to {0} for better performances.".format(
                  new_npixels_per_sample))

        if new_npixels_per_sample <= npixels_per_sample:
            return pmatrix, method, ('/detector', '/pixel'), None

        # if the actual number of pixels per sample is greater than
        # the specified one, redo the computation of the pointing matrix
        del pmatrix
        return self.get_pointing_matrix(pointing, header,
                                        new_npixels_per_sample, method)

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
