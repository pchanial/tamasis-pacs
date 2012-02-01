# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
from __future__ import division

import numpy as np
import os
import pyfits
import sys
import tamasisfortran as tmf

from pyoperators.operators_mpi import distribute_shape, distribute_slice
from pyoperators.utils import openmp_num_threads, strshape
from . import MPI
from .wcsutils import create_fitsheader

__all__ = []

def distribute_observation(detectors, observations, rank=None, comm=None):

    if comm is None:
        comm = MPI.COMM_WORLD
    
    size  = comm.Get_size()
    if size == 1:
        return detectors.copy(), list(observations)

    if rank is None:
        rank = comm.Get_rank()
    nthreads = openmp_num_threads()
    ndetectors = np.sum(~detectors)
    nobservations = len(observations)

    # number of observations. They should approximatively be of the same length
    nx = nobservations

    # number of detectors, grouped by the number of cpu cores
    ny = int(np.ceil(ndetectors / nthreads))

    # we start with the minimum blocksize and increase it until we find a
    # configuration that covers all the observations
    blocksize = int(np.ceil(nx * ny / size))
    while True:
        # by looping over x first, we favor larger numbers of detectors and
        # fewer numbers of observations per processor, to minimise inter-
        # processor communication in case of correlations between
        # detectors
        for xblocksize in range(1, blocksize+1):
            if blocksize / xblocksize != blocksize // xblocksize:
                continue
            yblocksize = int(blocksize // xblocksize)
            nx_block = int(np.ceil(nx / xblocksize))
            ny_block = int(np.ceil(ny / yblocksize))
            if nx_block * ny_block <= size:
                break
        if nx_block * ny_block <= size:
            break
        blocksize += 1

    ix = rank // ny_block
    iy = rank %  ny_block

    # check that the processor has something to do
    if ix >= nx_block:
        idetector = slice(0,0)
        iobservation = slice(0,0)
    else:
        idetector    = slice(iy * yblocksize * nthreads, (iy+1) * yblocksize * \
                             nthreads)
        iobservation = slice(ix * xblocksize, (ix+1) * xblocksize)

    detectors_ = detectors.copy()
    igood = np.where(~detectors_.ravel())[0]
    detectors_.ravel()[igood[0:idetector.start]] = True
    detectors_.ravel()[igood[idetector.stop:]] = True
    observations_ = observations[iobservation]

    return detectors_, observations_


#-------------------------------------------------------------------------------


def read_fits(filename, extname, comm):
    """
    Read and distribute a FITS file into local arrays.

    Parameters
    ----------
    filename : str
        The FITS file name.
    extname : str
        The FITS extension name. Use None to read the first HDU with data.
    comm : mpi4py.Comm
        The MPI communicator of the local arrays.
    """

    # check if the file name is the same for all MPI jobs
    files = comm.allgather(filename+str(extname))
    all_equal = all([f == files[0] for f in files])
    if comm.size > 1 and not all_equal:
        raise ValueError('The file name is not the same for all MPI jobs.')

    # get primary hdu or extension
    fits = pyfits.open(filename)
    if extname is not None:
        hdu = fits[extname]
    else:
        ihdu = 0
        while True:
            try:
                hdu = fits[ihdu]
            except IndexError:
                raise IOError('The FITS file has no data.')
            if hdu.header['NAXIS'] == 0:
                ihdu += 1
                continue
            if hdu.data is not None:
                break

    header = hdu.header
    n = header['NAXIS' + str(header['NAXIS'])]
    s = distribute_slice(n, comm=comm)
    output = pyfits.Section(hdu)[s]

    if not output.dtype.isnative:
        output = output.byteswap().newbyteorder('=')

    # update the header
    header['NAXIS' + str(header['NAXIS'])] = s.stop - s.start
    try:
        if header['CTYPE1'] == 'RA---TAN' and header['CTYPE2'] == 'DEC--TAN':
            header['CRPIX2'] -= s.start
    except KeyError:
        pass
    comm.Barrier()

    return output, header


#-------------------------------------------------------------------------------


def write_fits(filename, data, header, extension, extname, comm):
    """
    Write and combine local arrays into a FITS file.

    Parameters
    ----------
    filename : str
        The FITS file name.
    data : ndarray
        The array to be written.
    header : pyfits.Header
        The data FITS header. None can be set, in which case a minimal FITS
        header will be inferred from the data.
    extension : boolean
        If True, the data will be written as an extension to an already
        existing FITS file.
    extname : str
        The FITS extension name. Use None to write the primary HDU.
    comm : mpi4py.Comm
        The MPI communicator of the local arrays. Use MPI.COMM_SELF if the data
        are not meant to be combined into a global array. Make sure that the MPI
        processes are not executing this routine with the same file name.
    """

    # check if the file name is the same for all MPI jobs
    files = comm.allgather(filename+str(extname))
    all_equal = all(f == files[0] for f in files)
    if comm.size > 1 and not all_equal:
        raise ValueError('The file name is not the same for all MPI jobs.')
    ndims = comm.allgather(data.ndim)
    if any(n != ndims[0] for n in ndims):
        raise ValueError("The arrays have an incompatible number of dimensions:"
                         " '{0}'.".format(', '.join(str(n) for n in ndims)))
    ndim = ndims[0]
    shapes = comm.allgather(data.shape)
    if any(s[1:] != shapes[0][1:] for s in shapes):
        raise ValueError("The arrays have incompatible shapes: '{0}'.".format(
                         strshape(shapes)))

    # get header
    if header is None:
        header = create_fitsheader(fromdata=data, extname=extname)
    else:
        header = header.copy()
    if extname is not None:
        header.update('extname', extname)

    # we remove the file first to avoid an annoying pyfits informative message
    if not extension:
        try:
            os.remove(filename)
        except:
            pass

    # case without MPI communication
    if comm.size == 1:
        if not extension:
            hdu = pyfits.PrimaryHDU(data, header)
            hdu.writeto(filename, clobber=True)
        else:
            pyfits.append(filename, data, header)
        return

    # get global/local parameters
    shape_global = (sum(s[0] for s in shapes),) + shapes[0][1:]
    nglobal = shape_global[0]
    s = distribute_slice(nglobal)
    nlocal = s.stop - s.start

    if not extension:
        try:
            os.remove(filename)
        except:
            pass

    if comm.rank == 0:
        header['NAXIS' + str(ndim)] = nglobal
        shdu = pyfits.StreamingHDU(filename, header)
        data_loc = shdu._datLoc
        shdu.close()
    else:
        data_loc = None

    data_loc = comm.bcast(data_loc)

    # get a communicator excluding the processes which have no work to do
    # (Create_subarray does not allow 0-sized subarrays)
    chunk = np.product(shape_global[1:])
    rank_nowork = min(comm.size, nglobal)
    group = comm.Get_group()
    group.Incl(range(rank_nowork))
    newcomm = comm.Create(group)
    
    if comm.rank < rank_nowork:
        mtype = {1:MPI.BYTE, 4: MPI.FLOAT, 8:MPI.DOUBLE}[data.dtype.itemsize]
        ftype = mtype.Create_subarray([nglobal*chunk], [nlocal*chunk],
                                      [s.start*chunk])
        ftype.Commit()
        f = MPI.File.Open(newcomm, filename, amode=MPI.MODE_APPEND + \
                          MPI.MODE_WRONLY + MPI.MODE_CREATE)
        f.Set_view(data_loc, mtype, ftype, 'native', MPI.INFO_NULL)
        # mpi4py 1.2.2: pb with viewing data as big endian KeyError '>d'
        if sys.byteorder == 'little' and data.dtype.byteorder == '=' or \
           data.dtype.byteorder == '<':
            data = data.byteswap().newbyteorder('=')
        f.Write_all(data[0:nlocal])
        f.Close()

    if comm.rank == 0:
        shdu._ffo = pyfits.file._File(filename, 'append')
        shdu._ffo.getfile().seek(0,2)
        pyfitstype = {8:'uint8', 16:'int16', 32:'int32', 64:'int64', -32:'float32', -64:'float64'}[header['BITPIX']]
        completed = shdu.write(np.empty(0, dtype=pyfitstype))
        shdu.close()
        if not completed:
            raise RuntimeError('File is not completely written')

    comm.Barrier()
