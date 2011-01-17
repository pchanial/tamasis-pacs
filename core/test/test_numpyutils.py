import numpy
from tamasis.numpyutils import *

class TestFailure(Exception):
    def __init__(self, *args):
        if args is not None:
            print('Test failure:\n' + ',\n'.join([repr(a) for a in args]))
        raise Exception()

NaN = numpy.nan
if any_neq(1, 1+1.e-15): raise TestFailure()
if any_neq([1], [1+1.e-15]): raise TestFailure()
if not any_neq(1, 1+1.e-13): raise TestFailure()
if not any_neq([1], [1+1.e-13]): raise TestFailure()
if any_neq(1, numpy.asarray(1+1.e-8, dtype='float32')): raise TestFailure()
if any_neq([1], [numpy.asarray(1+1.e-8, dtype='float32')]): raise TestFailure()
if not any_neq(1, numpy.asarray(1+1.e-6, dtype='float32')): raise TestFailure()
if not any_neq([1], [numpy.asarray(1+1.e-6, dtype='float32')]): raise TestFailure()
if any_neq(NaN, NaN): raise TestFailure()
if any_neq([NaN], [NaN]): raise TestFailure()
if any_neq([NaN,1], [NaN,1]): raise TestFailure()
if not any_neq([NaN,1,NaN], [NaN,1,3]): raise TestFailure()
if any_neq(minmax([NaN, 1., 4., NaN, 10.]), [1., 10.]): raise TestFailure()

class A(numpy.ndarray):
    __slots__ = ('__dict__', 'info1')
    def __new__(cls, data, dtype=None, info1=None, info2=None, info3=None):
        data = numpy.array(data, dtype=dtype).view(cls)
        data.info1 = info1
        data.info2 = info2
        data.info3 = info3
        return data
    def __array_finalize__(self, obj):
        from copy import deepcopy
        if not hasattr(self, 'info1'):
            self.info1 = deepcopy(getattr(obj, 'info1', None))
        if not hasattr(self, 'info2'):
            self.info2 = deepcopy(getattr(obj, 'info2', None))
        if not hasattr(self, 'info3'):
            self.info3 = deepcopy(getattr(obj, 'info3', None))
    def copy(self):
        import copy
        return copy.deepcopy(self)
    def __repr__(self):
        desc="""\
array(data=%s,
      info1=%s,
      info2=%s,
      info3=%s"""
        return desc % (numpy.ndarray.__repr__(self), repr(self.info1), repr(self.info2), repr(self.info3))

dtypes=(numpy.bool, numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64, numpy.complex64, numpy.complex128)

for dtype in dtypes:
    arr = numpy.ones((2,3), dtype=dtype)
    a = A(arr, info1='1', info2=True)
    b1 = A(arr, info1='1', info2=True)
    b1[0,1] = 0
    b2 = A(arr, info1='2', info2=True)
    b3 = A(arr, info1='1', info2=False)

    for b in (b1, b2, b3):
        if not any_neq(a, b): raise TestFailure(dtype, a, b)
    b = a.copy()
    if any_neq(a,b): raise TestFailure(dtype, a, b)

    a.info3 = b
    b1 = a.copy(); b1.info3[0,1] = 0
    b2 = a.copy(); b2.info3.info1 = '2'
    b3 = a.copy(); b3.info3.info2 = False
    for b in (b1, b2, b3):
        if not any_neq(a, b): raise TestFailure(dtype, a, b)
    b = a.copy()
    if any_neq(a,b): raise TestFailure(dtype, a, b)
