import numpy
from tamasis.numpyutils import any_neq

class TestFailure(Exception):
    pass

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

print 'OK.'
