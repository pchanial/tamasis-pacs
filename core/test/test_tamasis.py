import numpy as np
import tamasis
from tamasis import *

class TestFailure(Exception):
    pass

a = np.ones((10,10))
b = a.copy()
tamasis.tmf.add_inplace(a.T, b.T)
if any_neq(a, 2): raise TestFailure()

b[:] = -1
tamasis.tmf.subtract_inplace(a.T, b.T)
if any_neq(a, 3): raise TestFailure()

b[:] = 2
tamasis.tmf.multiply_inplace(a.T, b.T)
if any_neq(a, 6): raise TestFailure()

tamasis.tmf.divide_inplace(a.T, b.T)
if any_neq(a, 3): raise TestFailure()

a1 = np.ones((2,3,4,5))
a2 = np.ones((2,3,4,5))
b = np.arange(2*3*4, dtype=tamasis.var.FLOAT_DTYPE).reshape((2,3,4))
tamasis.tmf.multiply_inplace(a1.T, b.T)
a2.T[:] *= b.T

if any_neq(a1, a2): raise TestFailure()
