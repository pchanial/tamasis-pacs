import numpy as np
import os
import tamasisfortran as tmf
from mpi4py import MPI

mpi_comm = MPI.COMM_WORLD
path = os.path.abspath(os.path.dirname(__file__) + '/../../../../share/tamasis')
verbose = False

FLOAT_DTYPE = {
    4  : np.dtype(np.float32),
    8  : np.dtype(np.float64),
    16 : np.dtype(np.float128),
}[tmf.info_nbytes_real()]

COMPLEX_DTYPE = {
    4  : np.dtype(np.complex64),
    8  : np.dtype(np.complex128),
    16 : np.dtype(np.complex256),
}[tmf.info_nbytes_real()]

VERSION = tmf.info_version().strip()

def get_default_dtype(data):
    import numpy
    if np.iscomplexobj(data):
        return COMPLEX_DTYPE
    else:
        return FLOAT_DTYPE

del MPI, os, tmf
