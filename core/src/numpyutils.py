# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#

from __future__ import division

import numpy as np

np.seterr(all='ignore')

__all__ = [ 'minmax' ]

def assert_all_eq(a, b, rtol=None, atol=0., msg=None):
    assert all_eq(a, b, rtol=rtol, atol=atol), msg

def all_eq(a, b, rtol=None, atol=0.):
    return not any_neq(a, b, rtol=rtol, atol=atol)

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
            print('First argument is a dict and the second one is not.')
            return True
        if set([k for k in a]) != set([k for k in b]):
            print('Arguments are dictionaries of different items.')
            return True
        for k in a:
            if any_neq(a[k], b[k]):
                print('Arguments are dictionaries of different values')
                return True
        return False

    # restrict comparisons to some classes
    class_ok = (np.ndarray, int, float, complex, list, tuple, str, unicode)
    if not isinstance(a, class_ok) or not isinstance(b, class_ok):
        return not isinstance(a, type(b)) and a is not None and \
               not isinstance(b, type(a)) and b is not None

    a = np.asanyarray(a)
    b = np.asanyarray(b)

    # get common base class to give some slack
    cls = type(b)
    while True:
        if isinstance(a, cls):
            break
        cls = cls.__base__

    a = a.view(cls)
    b = b.view(cls)
    
    akind = a.dtype.kind
    bkind = b.dtype.kind

    aattr = get_attributes(a)
    battr = get_attributes(b)
    if set(aattr) != set(battr):
        print('The arguments do not have the same attributes.')
        return True

    for k in aattr:
        if any_neq(getattr(a,k), getattr(b,k), rtol, atol):
            print("The argument attributes '" + k + "' differ.")
            return True

    # then compare the values of the array
    if akind not in 'biufc' or bkind not in 'biufc':
        if akind != bkind:
            print('The argument data types are incompatible.')
            return True
        if akind == 'S':
            result = np.any(a != b)
            if result:
                print('String arguments differ.')
            return result
        if akind == 'V':
            if set(a.dtype.names) != set(b.dtype.names):
                print('The names of the argument dtypes are different.')
                return True
            result = False
            for name in a.dtype.names:
                result_ = any_neq(a[name], b[name])
                if result_:
                    print("Values for '{0}' are different.".format(name))
                result = result or result_
            return result
        else:
            raise NotImplementedError('Kind ' + akind + ' is not implemented.')
    
    if akind in 'biu' and bkind in 'biu':
        result = np.any(a != b)
        if result:
            print('Integer arguments differ.')
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
            raise NotImplementedError('The inputs are not in single or double' \
                                      ' precision.')

    mask = np.isnan(a)
    if np.any(mask != np.isnan(b)):
        print('Argument NaNs differ.')
        return True
    if np.all(mask):
        return False

    result = abs(a-b) > rtol * np.maximum(abs(a), abs(b)) + atol
    if not np.isscalar(result):
        result = np.any(result[~mask])

    if result:
        factor = np.nanmax(abs(a-b) / (rtol * np.maximum(abs(a),abs(b)) + atol))
        print('Argument data difference exceeds tolerance by factor ' + \
              str(factor) + '.')

    return result


#-------------------------------------------------------------------------------


def minmax(v):
    """Returns min and max values of an array, discarding NaN values."""
    v = np.asanyarray(v)
    if np.all(np.isnan(v)):
        return np.array([np.nan, np.nan])
    return np.array((np.nanmin(v), np.nanmax(v)))
   
   
#-------------------------------------------------------------------------------


def isscalar(data):
    """Hack around np.isscalar oddity"""
    return data.ndim == 0 if isinstance(data, np.ndarray) else np.isscalar(data)


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
    return sorted(attributes)


#-------------------------------------------------------------------------------


def get_type(data):
    """Returns input's data type."""
    data_ = np.asarray(data)
    type_ = data_.dtype.type.__name__
    if type_[-1] == '_':
        type_ = type_[0:-1]
    if type_ != 'object':
        return type_
    return type(data).__name__


#-------------------------------------------------------------------------------


def _my_issctype(dtype):
    """Hack around np.issctype bug"""
    return np.issctype(dtype) and str(dtype)[0:2] != '|S'
