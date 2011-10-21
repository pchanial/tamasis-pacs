# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
from __future__ import division

import math
import numpy as np
import os
import pyfits
import sys
import tamasisfortran as tmf
from mpi4py import MPI
from .wcsutils import create_fitsheader

__all__ = []


def split_work(nglobal, rank=None, comm=None):
    comm = comm or MPI.COMM_WORLD
    size  = comm.Get_size()
    if rank is None:
        rank = comm.Get_rank()
    nlocal = int(np.ceil(nglobal / size))
    return slice(min(rank * nlocal, nglobal), min((rank + 1) * nlocal, nglobal))
    

#-------------------------------------------------------------------------------


def split_shape(shape, comm=None):
    comm = comm or MPI.COMM_WORLD
    if not isinstance(shape, (list, np.ndarray, tuple)):
        shape = (shape,)
    n = int(np.ceil(shape[0] / comm.Get_size()))
    return (n,) + tuple(shape[1:])

#-------------------------------------------------------------------------------


def split_observation(detectors, observations, rank=None, comm=None):

    if comm is None:
        comm = MPI.COMM_WORLD
    
    size  = comm.Get_size()
    if size == 1:
        return detectors.copy(), list(observations)

    if rank is None:
        rank = comm.Get_rank()
    nthreads = tmf.info_nthreads()
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

def read_fits(hdu, comm=None, sequential=False):
    """Read a FITS file into local images of the same shape"""

    if comm is None:
        comm = MPI.COMM_SELF
    header = hdu.header
    naxis = [header['NAXIS'+str(i)] for i in range(1, header['NAXIS']+1)]
    s = split_work(naxis[-1], comm=comm)
    n = s.stop - s.start
    nmax = int(np.ceil(naxis[-1] / comm.Get_size()))
    shape_extra = list(reversed(naxis))
    shape_extra[0] = nmax - n
    fitsdtype = {8:'uint8', 16:'>i2', 32:'>i4', 64:'>i8', -32:'>f4', -64:'>f8'}

    for iproc in range(comm.Get_size()):
        if iproc == comm.Get_rank():
            if n == 0:
                output = np.zeros(shape_extra,dtype=fitsdtype[header['BITPIX']])
            else:
                output = pyfits.Section(hdu)[s]
                if n < nmax:
                    extra = np.zeros(shape_extra, dtype=output.dtype)
                    output = np.concatenate([output, extra])
        if sequential:
            comm.Barrier()
    comm.Barrier()

    return output, header, tuple(naxis)[::-1]

def write_fits(filename, data, header, shape_global, extension, comm,
               extname=None):
    """Write local images into a FITS file"""

    if comm is None:
        comm = MPI.COMM_SELF

    if not extension:
        try:
            os.remove(filename)
        except:
            pass

    if header is None:
        header = create_fitsheader(shape_global[::-1], data.dtype,
                                   extname=extname)

    rank = comm.Get_rank()
    size = comm.Get_size()
    files = comm.allgather(filename)
    allsame = all([f == files[0] for f in files])
    alldiff = len(files) == len(np.unique(files))
    if not alldiff and not allsame:
        raise ValueError('Some target filenames are equal, but not all.')
    if alldiff or size == 1:
        if not extension:
            hdu = pyfits.PrimaryHDU(data, header)
            hdu.writeto(filename, clobber=True)
        else:
            pyfits.append(filename, data, header)
        return

    if rank == 0:
        shdu = pyfits.StreamingHDU(filename, header)
        data_loc = shdu._datLoc
        shdu.close()
    else:
        data_loc = None
    data_loc = comm.bcast(data_loc)

    # get a communicator excluding the processes which have no work to do
    # (Create_subarray does not allow 0-sized subarrays)
    nglobal = shape_global[0]
    chunk = np.product(shape_global[1:])
    s = split_work(nglobal)
    nlocal = s.stop - s.start
    nmax = int(np.ceil(nglobal / size))
    rank_nowork = int(np.ceil(nglobal / nmax))
    group = comm.Get_group()
    group.Incl(range(rank_nowork))
    newcomm = comm.Create(group)
    
    if rank < rank_nowork:
        mtype = {1:MPI.BYTE, 4: MPI.FLOAT, 8:MPI.DOUBLE}[data.dtype.itemsize]
        ftype = mtype.Create_subarray([nglobal*chunk], [nlocal*chunk],
                                      [s.start*chunk])
        ftype.Commit()
        f = MPI.File.Open(newcomm, filename, amode=MPI.MODE_APPEND + \
                          MPI.MODE_WRONLY + MPI.MODE_CREATE)
        f.Set_view(data_loc, mtype, ftype, 'native', MPI.INFO_NULL)
        # mpi4py 1.2.2: pb with viewing data as big endian KeyError '>d'
        if sys.byteorder == 'little' and data.dtype.byteorder in ('=', '<'):
            data = data.byteswap()
        else:
            data = data.newbyteorder('=')
        f.Write_all(data[0:nlocal])
        f.Close()

    if rank == 0:
        shdu._ffo = pyfits.core._File(filename, 'append')
        shdu._ffo.getfile().seek(0,2)
        pyfitstype = {8:'uint8', 16:'int16', 32:'int32', 64:'int64', -32:'float32', -64:'float64'}[header['BITPIX']]
        completed = shdu.write(np.empty(0, dtype=pyfitstype))
        shdu.close()
        if not completed:
            raise RuntimeError('File is not completely written')

    comm.Barrier()
