# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
import numpy as np
import tamasisfortran as tmf

__all__ = []

def split_observation(comm, detectors, observations):

    size  = comm.Get_size()
    if size == 1:
        return detectors.copy(), list(observations)

    rank = comm.Get_rank()
    nthreads = tmf.info_nthreads()
    ndetectors = np.sum(~detectors)
    nobservations = len(observations)

    # number of observations. They should approximatively be of the same length
    nx = nobservations

    # number of detectors, grouped by the number of cpu cores
    ny = int(np.ceil(float(ndetectors) / nthreads))

    # we start with the minimum blocksize and increase it until we find a
    # configuration that covers all the observations
    blocksize = int(np.ceil(float(nx * ny) / size))
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
