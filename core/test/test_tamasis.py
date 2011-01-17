import numpy
from tamasis import *

class TestFailure(Exception):
    pass

a = numpy.ones((10,10))
b = a.copy()

tmf.add_inplace(a.ravel(), b.ravel())
if any_neq(a, 2): raise TestFailure()

b[:] = -1
tmf.subtract_inplace(a.ravel(), b.ravel())
if any_neq(a, 3): raise TestFailure()

b[:] = 2
tmf.multiply_inplace(a.ravel(), b.ravel())
if any_neq(a, 6): raise TestFailure()

tmf.divide_inplace(a.ravel(), b.ravel())
if any_neq(a, 3): raise TestFailure()
