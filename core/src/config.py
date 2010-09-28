import os
__verbose__ = False
tamasis_dir = os.path.dirname(__file__) + '/../' if os.path.dirname(__file__) != '' else '../'
del os

__version_info__ = (1, 4, 0)
__version__ = '.'.join((str(i) for i in __version_info__))

__all__ = [ 'get_default_dtype', 'get_default_dtype_complex', 'get_default_dtype_float', 'tamasis_dir', '__verbose__', 
            '__version__', '__version_info__' ]

def get_default_dtype(data):
    import numpy
    iscomplex = isinstance(data, numpy.ndarray) and data.dtype.type in (numpy.complex64, numpy.complex128, numpy.complex256) or \
        type(data) is complex
    if iscomplex:
        return get_default_dtype_complex()
    else:
        return get_default_dtype_float()

def get_default_dtype_complex():
    import numpy
    import tamasisfortran as tmf
    nbytes = tmf.info_nbytes_real()
    if nbytes == 4:
        return numpy.complex64
    elif nbytes == 8:
        return numpy.complex128
    elif nbytes == 16:
        return numpy.complex256
    else:
        raise ValueError("Invalid number of bytes per real '"+str(nbytes)+"'.")

def get_default_dtype_float():
    import numpy
    import tamasisfortran as tmf
    nbytes = tmf.info_nbytes_real()
    if nbytes == 4:
        return numpy.float32
    elif nbytes == 8:
        return numpy.float64
    elif nbytes == 16:
        return numpy.float128
    else:
        raise ValueError("Invalid number of bytes per real '"+str(nbytes)+"'.")
