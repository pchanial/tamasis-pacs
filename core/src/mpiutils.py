# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
import math
import numpy as np
import tamasisfortran as tmf
from . import var
from mpi4py import MPI
from scipy.sparse.linalg import aslinearoperator

__all__ = []

def gradient(f):
    """Parallelised three-point estimation of the gradient"""
    def df(x, out=None, comm=None):
        x = np.array(x, copy=False)
        if out is None:
            out = np.empty(x.size)
        elif len(out.shape) == 0:
            out.shape = (1,)
        if comm is None:
            comm = MPI.COMM_WORLD
        dx = np.sqrt(np.finfo(float).eps) * x
        rank = comm.Get_rank()
        size = comm.Get_size()
        for r in range(size):
            n = comm.bcast(x.size, root=r)
            if r == rank:
                x_dx = x.copy()
                for i in xrange(n):
                    x_dx[i] = x[i] + dx[i]
                    try:
                        f2 = f(x_dx, comm=comm)
                    except TypeError:
                        f2 = f(x_dx)
                    x_dx[i] = x[i] - dx[i]
                    try:
                        f1 = f(x_dx, comm=comm)
                    except TypeError:
                        f1 = f(x_dx)
                    out[i] = (f2-f1) / (2*dx[i])
                    x_dx[i] = x[i]
            else:
                for i in xrange(n):
                    try:
                        f(x, comm=comm)
                        f(x, comm=comm)
                    except TypeError:
                        f(x)
                        f(x)
        out.shape = x.shape
        return out
    return df

#-------------------------------------------------------------------------------


def sum(array, comm=None):
    """Parallelised sum"""
    array = np.ascontiguousarray(array)
    if comm is None:
        comm = MPI.COMM_WORLD
    output = np.array(tmf.sum(array.T))
    comm.Allreduce(MPI.IN_PLACE, [output, MPI.DOUBLE], op=MPI.SUM)
    return float(output)


#-------------------------------------------------------------------------------


def norm_l1(array, comm=None):
    """Parallelised L-1 norm"""
    array = np.ascontiguousarray(array)
    if comm is None:
        comm = MPI.COMM_WORLD
    return float(tmf.norm_l1(array.T, comm.py2f()))


#-------------------------------------------------------------------------------


def dnorm_l1(array, out=None, comm=None):
    """Gradient of the L-1 norm"""
    array = np.ascontiguousarray(array)
    if out is None:
        out = np.empty(array.shape)
    tmf.dnorm_l1(array.T, out.T)
    return out


#-------------------------------------------------------------------------------


def norm_l2(array, comm=None):
    """Parallelised euclidian norm"""
    array = np.ascontiguousarray(array)
    if comm is None:
        comm = MPI.COMM_WORLD
    return math.sqrt(float(tmf.norm2(array.T, comm.py2f())))


#-------------------------------------------------------------------------------


def dnorm_l2(array, out=None, comm=None):
    """Gradient of the euclidian norm"""
    array = np.ascontiguousarray(array)
    if out is None:
        out = np.empty(array.shape)
    if comm is None:
        comm = MPI.COMM_WORLD
    tmf.dnorm_l2(array.T, comm.py2f(), out.T)
    return out


#-------------------------------------------------------------------------------


def norm_lp(p):
    """Parallelised L-p norm"""
    def norm(array, comm=None):
        array = np.ascontiguousarray(array)
        if comm is None:
            comm = MPI.COMM_WORLD
        return float(tmf.normp(array.T, p, comm.py2f()))**(1./p)
    return norm


#-------------------------------------------------------------------------------


def dnorm_lp(p):
    """Gradient of the L-p norm"""
    def dnorm(array, out=None, comm=None):
        array = np.ascontiguousarray(array)
        if out is None:
            out = np.empty(array.shape)
        if comm is None:
            comm = MPI.COMM_WORLD
        tmf.dnorm_lp(array.T, p, comm.py2f(), out.T)
        return out
    return dnorm


#-------------------------------------------------------------------------------


def norm_linf(array, comm=None):
    """Parallelised L-inf norm"""
    array = np.ascontiguousarray(array)
    if comm is None:
        comm = MPI.COMM_WORLD
    return float(tmf.norm_linf(array.T, comm.py2f()))


#-------------------------------------------------------------------------------


def norm_huber(delta):
    """Parallelised Huber's norm"""
    if delta < 0:
        raise ValueError("Huber's norm delta must be positive.")
    if delta == 0:
        return norm2

    def norm(array, comm=None):
        array = np.ascontiguousarray(array)
        if comm is None:
            comm = MPI.COMM_WORLD
        return float(tmf.norm_huber(array.T, delta, comm.py2f()))
    return norm


#-------------------------------------------------------------------------------


def dnorm_huber(delta):
    """Gradient of the Huber's norm"""
    if delta < 0:
        raise ValueError("Huber's norm delta must be positive.")
    if delta == 0:
        return dnorm2

    def dnorm(array, out=None, comm=None):
        array = np.ascontiguousarray(array)
        if out is None:
            out = np.empty(array.shape)
        tmf.dnorm_huber(array.T, delta, out.T)
        return out
    return dnorm
        

#-------------------------------------------------------------------------------


def norm2(array, comm=None):
    """Parallelised squared euclidian norm"""
    array = np.ascontiguousarray(array)
    if comm is None:
        comm = MPI.COMM_WORLD    
    return float(tmf.norm2(array.T, comm.py2f()))


#-------------------------------------------------------------------------------


def dnorm2(array, out=None, comm=None):
    """Gradient of the squared euclidian norm"""
    array = np.ascontiguousarray(array)
    if out is None:
        out = np.empty(array.shape)
    tmf.dnorm2(array.T, out.T)
    return out


#-------------------------------------------------------------------------------


def normp(p):
    """Parallelised L-p norm to the power of p"""
    def norm(array, comm=None):
        array = np.ascontiguousarray(array)
        if comm is None:
            comm = MPI.COMM_WORLD    
        return float(tmf.normp(array.T, p, comm.py2f()))
    return norm


#-------------------------------------------------------------------------------


def dnormp(p):
    """Gradient of the L-p norm to the power of pm"""
    def dnorm(array, out=None, comm=None):
        array = np.ascontiguousarray(array)
        if out is None:
            out = np.empty(array.shape)
        tmf.dnormp(array.T, p, out.T)
        return out
    return dnorm


#-------------------------------------------------------------------------------


def norm2_ellipsoid(A):
    """Ellipsoid norm

    Returns x^T A x, where A is a definite positive symmetric matrix
    """
    A = aslinearoperator(A)
    def norm(x, comm=None):
        if comm is None:
            comm = MPI.COMM_WORLD
        return dot(x, A.matvec(x.ravel()), comm=comm)
    return norm


#-------------------------------------------------------------------------------


def dnorm2_ellipsoid(A):
    """Derivative of the ellipsoid norm

    Returns 2 * A x, where A is a definite positive symmetric matrix
    """
    A = aslinearoperator(A)
    def dnorm(x, out=None, comm=None):
        if comm is None:
            comm = MPI.COMM_WORLD
        result = 2 * A.matvec(x.ravel())
        if out is not None:
            out[:] = result
        return result
    return dnorm


#-------------------------------------------------------------------------------


def dot(array1, array2, comm=None):
    """Parallelised dot product"""
    array1 = np.ascontiguousarray(array1)
    array2 = np.ascontiguousarray(array2)
    if comm is None:
        comm = MPI.COMM_WORLD
    output = np.array(tmf.dot(array1.T, array2.T))
    comm.Allreduce(MPI.IN_PLACE, [output, MPI.DOUBLE], op=MPI.SUM)
    return float(output)


#-------------------------------------------------------------------------------


def split_work(nglobal, rank=None, comm=None):
    if comm is None:
        comm = MPI.COMM_WORLD
    size  = comm.Get_size()
    if rank is None:
        rank = comm.Get_rank()
    nlocal = int(np.ceil(float(nglobal) / size))
    return slice(min(rank * nlocal, nglobal), min((rank + 1) * nlocal, nglobal))
    

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
