import numpy

__all__ = [ 'any_neq' ]

def any_neq(a, b, rtol=None, atol=0.):
    """
    Returns True if two arrays are element-wise equal within a tolerance.
    Differs from numpy.allclose in two aspects: the default rtol values (10^-7 
    and 10^-14 for single and double floats or complex) and the treatment of 
    NaN values (do not return False if the two arrays have element-wise
    both NaN value)
    """
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    doubletype = ('float64', 'complex128')
    if rtol is None:
        if str(a.dtype) in doubletype or str(b.dtype) in doubletype:
            rtol = 1.e-14
        else:
            rtol = 1.e-7
    mask = numpy.isnan(a)
    if numpy.any(mask != numpy.isnan(b)):
        return True
    if numpy.all(mask):
        return False
    result = abs(a-b) > rtol * numpy.maximum(abs(a), abs(b)) + atol
    if numpy.isscalar(mask):
        return result
    else:
        return numpy.any(result[~mask])


#-------------------------------------------------------------------------------


def get_type(data):
    """Returns input's data type."""
    data_ = numpy.asarray(data)
    type_ = data_.dtype.type.__name__
    if type_[-1] == '_':
        type_ = type_[0:-1]
    if type_ != 'object':
        return type_
    return type(data).__name__

   
#-------------------------------------------------------------------------------


def _my_issctype(dtype):
    """Hack around numpy.issctype bug"""
    return numpy.issctype(dtype) and str(dtype)[0:2] != '|S'
   
   
#-------------------------------------------------------------------------------


def _my_isscalar(data):
    """Hack around numpy.isscalar bug"""
    return numpy.rank(data) == 0 if isinstance(data, numpy.ndarray) else numpy.isscalar(data)
