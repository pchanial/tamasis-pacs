import numpy
from tamasis.utils import mean_degrees
from tamasis.numpyutils import any_neq

class TestFailure(Exception):
    pass

if any_neq(mean_degrees([1,2]), 1.5): raise TestFailure()
if any_neq(mean_degrees([1,359.1]), 0.05, rtol=1.e-12): raise TestFailure()
if any_neq(mean_degrees([0.1,359.1]), 359.6, rtol=1.e-12): raise TestFailure()

print 'OK.'
