import numpy as np
from tamasis.numpyutils import all_eq, any_neq, minmax

def test_any_neq1():
    assert all_eq(1, 1+1.e-15)
    assert all_eq([1], [1+1.e-15])
    assert any_neq(1, 1+1.e-13)
    assert any_neq([1], [1+1.e-13])
    assert all_eq(1, np.asarray(1+1.e-8, dtype='float32'))
    assert all_eq([1], [np.asarray(1+1.e-8, dtype='float32')])
    assert any_neq(1, np.asarray(1+1.e-6, dtype='float32'))
    assert any_neq([1], [np.asarray(1+1.e-6, dtype='float32')])

def test_any_neq2():
    NaN = np.nan
    assert all_eq(NaN, NaN)
    assert all_eq([NaN], [NaN])
    assert all_eq([NaN,1], [NaN,1])
    assert any_neq([NaN,1,NaN], [NaN,1,3])
    assert all_eq(minmax([NaN, 1., 4., NaN, 10.]), [1., 10.])

class A(np.ndarray):
    __slots__ = ('__dict__', 'info1')
    def __new__(cls, data, dtype=None, info1=None, info2=None, info3=None):
        data = np.array(data, dtype=dtype).view(cls)
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
        return desc % (np.ndarray.__repr__(self), repr(self.info1), repr(self.info2), repr(self.info3))

dtypes=(np.bool, np.int8, np.int16, np.int32, np.int64, np.float32, np.float64, np.complex64, np.complex128)

def test_any_neq3():
    for dtype in dtypes:
        arr = np.ones((2,3), dtype=dtype)
        a = A(arr, info1='1', info2=True)
        b1 = A(arr, info1='1', info2=True)
        b1[0,1] = 0
        b2 = A(arr, info1='2', info2=True)
        b3 = A(arr, info1='1', info2=False)
        def func1(a, b):
            assert any_neq(a, b)
        for b in (b1, b2, b3):
            yield func1, a, b
        b = a.copy()
        assert all_eq(a,b)
        a.info3 = b
        b1 = a.copy(); b1.info3[0,1] = 0
        b2 = a.copy(); b2.info3.info1 = '2'
        b3 = a.copy(); b3.info3.info2 = False
        def func2(a, b):
            assert any_neq(a, b)
        for b in (b1, b2, b3):
            yield func2, a, b
        b = a.copy()
        assert all_eq(a,b)
