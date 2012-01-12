# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
import numpy as np
import os
import tamasisfortran as tmf
from . import MPI

comm_tod = MPI.COMM_WORLD
comm_map = MPI.COMM_SELF

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
    if np.iscomplexobj(data):
        return COMPLEX_DTYPE
    else:
        return FLOAT_DTYPE

del MPI, os, tmf
