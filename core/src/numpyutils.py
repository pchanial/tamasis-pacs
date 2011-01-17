import numpy
from . import var

numpy.seterr(all='ignore')

__all__ = [ 'any_neq', 'minmax' ]

def any_neq(a, b, rtol=None, atol=0.):
    """
    Returns True if two arrays are element-wise equal within a tolerance.
    Differs from numpy.allclose in two aspects: the default rtol values (10^-7 
    and 10^-14 for single and double floats or complex) and the treatment of 
    NaN values (do not return False if the two arrays have element-wise
    both NaN value)
    """

    # for dictionaries, look up the items
    if isinstance(a, dict):
        if not isinstance(b, dict):
            if var.verbose: print('First argument is a dict and the second one is not.')
            return True
        if set([k for k in a]) != set([k for k in b]):
            if var.verbose: print('Input argument are dictionaries of different items.')
            return True
        for k in a:
            if any_neq(a[k], b[k]):
                if var.verbose: print('Input arguments are dictionaries of different values')
                return True
        return False

    a = numpy.asanyarray(a)
    b = numpy.asanyarray(b)

    # get common base class to give some slack
    cls = type(b)
    while True:
        if isinstance(a, cls):
            break
        cls = cls.__base__

    a = a.view(cls)
    b = b.view(cls)
    
    # first, lookup __slots__ and __dict__
    akind = a.dtype.kind
    bkind = b.dtype.kind

    # do not compare objects
    if akind == 'O' or bkind == 'O':
        return False

    aattr = get_attributes(a)
    battr = get_attributes(b)
    if set(aattr) != set(battr):
        if var.verbose: print('The arguments do not have the same attributes.')
        return True

    for k in aattr:
        if any_neq(getattr(a,k), getattr(b,k), rtol, atol):
            if var.verbose: print("The argument attributes '" + k + "' differ.")
            return True

    # then compare the values of the array
    if akind not in 'bifc' or bkind not in 'bifc':
        if akind != bkind:
            if var.verbose: print('The argument data types are incompatible.')
            return True
        if akind == 'S':
            result = numpy.any(a != b)
            if result and var.verbose: print('String arguments differ.')
            return result
        else:
            raise NotImplemented('Kind ' + akind + ' is not implemented.')
    
    if akind in 'bi' and bkind in 'bi':
        result = numpy.any(a != b)
        if result and var.verbose: print('Integer arguments differ.')
        return result
    
    if rtol is None:
        asize = a.dtype.itemsize // 2 if akind == 'c' else a.dtype.itemsize
        bsize = b.dtype.itemsize // 2 if bkind == 'c' else b.dtype.itemsize
        precision = min(asize, bsize)
        if precision == 8:
            rtol = 1.e-14
        elif precision == 4:
            rtol = 1.e-7
        else:
            raise NotImplementedError('The inputs are not in single or double precision.')

    mask = numpy.isnan(a)
    if numpy.any(mask != numpy.isnan(b)):
        if var.verbose: print('Argument NaNs differ.')
        return True
    if numpy.all(mask):
        return False
    result = abs(a-b) > rtol * numpy.maximum(abs(a), abs(b)) + atol
    if var.verbose: print('Argument data differ.')
    if numpy.isscalar(mask):
        return result
    else:
        return numpy.any(result[~mask])


#-------------------------------------------------------------------------------


def minmax(v):
    """Returns min and max values of an array, discarding NaN values."""
    v = numpy.asanyarray(v)
    return numpy.array((numpy.nanmin(v), numpy.nanmax(v)))


#-------------------------------------------------------------------------------


def get_attributes(obj):
    try:
        attributes = [ k for k in obj.__dict__ ]
    except AttributeError:
        attributes = []
    if hasattr(obj.__class__, '__mro__'):
        for cls in obj.__class__.__mro__:
            for slot in cls.__dict__.get('__slots__', ()):
                if hasattr(obj, slot):
                    attributes.append(slot)
    return attributes


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
