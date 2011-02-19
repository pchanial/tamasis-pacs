import numpy as np
import tamasisfortran as tmf

__all__ = []

def split_observation(comm, ndetectors, nobservations):

    rank = comm.Get_rank()
    nnodes  = comm.Get_size()
    nthreads = tmf.info_nthreads()

    # number of observations. They should approximatively be of the same length
    nx = nobservations

    # number of detectors, grouped by the number of cpu cores
    ny = int(np.ceil(float(ndetectors) / nthreads))

    # we start with the miminum blocksize and increase it until we find a
    # configuration that covers all the observations
    blocksize = int(np.ceil(float(nx * ny) / nnodes))
    while True:
        # by looping over x first, we favor larger numbers of detectors and
        # fewer numbers of observations per processor, to minimise inter-
        # processor communication in case of correlations between
        # detectors
        for xblocksize in range(1, blocksize+1):
            if float(blocksize) / xblocksize != blocksize // xblocksize:
                continue
            yblocksize = int(blocksize // xblocksize)
            nx_block = int(np.ceil(float(nx) / xblocksize))
            ny_block = int(np.ceil(float(ny) / yblocksize))
            if nx_block * ny_block <= nnodes:
                break
        if nx_block * ny_block <= nnodes:
            break
        blocksize += 1

    ix = rank // ny_block
    iy = rank %  ny_block

    # check that the processor has something to do
    if ix >= nx_block:
        iobservation = slice(0,0)
        idetector = slice(0,0)
    else:
        iobservation = slice(ix * xblocksize, (ix+1) * xblocksize)
        idetector    = slice(iy * yblocksize * nthreads, (iy+1) * yblocksize * \
                             nthreads)

    return iobservation, idetector
        

#-------------------------------------------------------------------------------
