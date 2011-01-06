import numpy
import os
import tamasisfortran as tmf
from mpi4py import MPI

mpi_comm = MPI.COMM_WORLD
path = os.path.abspath(os.path.dirname(__file__) + '/../../../../share/tamasis')
verbose = False

FLOAT_DTYPE = {
    4  : numpy.dtype(numpy.float32),
    8  : numpy.dtype(numpy.float64),
    16 : numpy.dtype(numpy.float128),
}[tmf.info_nbytes_real()]

COMPLEX_DTYPE = {
    4  : numpy.dtype(numpy.complex64),
    8  : numpy.dtype(numpy.complex128),
    16 : numpy.dtype(numpy.complex256),
}[tmf.info_nbytes_real()]

VERSION = tmf.info_version().strip()

def get_default_dtype(data):
    import numpy
    if numpy.iscomplexobj(data):
        return COMPLEX_DTYPE
    else:
        return FLOAT_DTYPE

del MPI, os, tmf
