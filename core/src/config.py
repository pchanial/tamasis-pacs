import os
__verbose__ = False
tamasis_dir = os.path.dirname(__file__) + '/../' if os.path.dirname(__file__) != '' else '../'
del os

__version_info__ = (1, 8, 0)
__version__ = '.'.join((str(i) for i in __version_info__))

__all__ = [ 'get_default_dtype', 'get_default_dtype_complex', 'get_default_dtype_float', 'tamasis_dir', '__verbose__', 
            '__version__', '__version_info__' ]

def get_default_dtype(data):
    import numpy
    if numpy.iscomplexobj(data):
        return get_default_dtype_complex()
    else:
        return get_default_dtype_float()

def get_default_dtype_complex():
    import numpy
    import tamasisfortran as tmf
    nbytes = tmf.info_nbytes_real()
    if nbytes == 4:
        return numpy.dtype(numpy.complex64)
    elif nbytes == 8:
        return numpy.dtype(numpy.complex128)
    elif nbytes == 16:
        return numpy.dtype(numpy.complex256)
    else:
        raise ValueError("Invalid number of bytes per real '"+str(nbytes)+"'.")

def get_default_dtype_float():
    import numpy
    import tamasisfortran as tmf
    nbytes = tmf.info_nbytes_real()
    if nbytes == 4:
        return numpy.dtype(numpy.float32)
    elif nbytes == 8:
        return numpy.dtype(numpy.float64)
    elif nbytes == 16:
        return numpy.dtype(numpy.float128)
    else:
        raise ValueError("Invalid number of bytes per real '"+str(nbytes)+"'.")
