import numpy
from tamasis.utils import MaskPolicy, any_neq, mean_degrees

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

if any_neq(mean_degrees([1,2]), 1.5): raise TestFailure()
if any_neq(mean_degrees([1,359.1]), 0.05, rtol=1.e-12): raise TestFailure()
if any_neq(mean_degrees([0.1,359.1]), 359.6, rtol=1.e-12): raise TestFailure()

flags = ['bad', 'u1', 'u2']
good_policy = ['kEep', 'removE', 'MASK']
bad_policy = ['remove', 'KKeep']

mask_policy = MaskPolicy(flags, good_policy)
if numpy.any(numpy.array(mask_policy) != (0,2,1)): raise TestFailure()
try:
    junk = MaskPolicy(flags, bad_policy)
except ValueError:
    pass
try:
    junk = MaskPolicy(flags[0:2], bad_policy)
except KeyError:
    pass

print 'OK.'
