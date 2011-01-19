import numpy
from tamasis.utils import *
from tamasis.utils import _distance_slow
from tamasis.numpyutils import any_neq

class TestFailure(Exception):
    pass

if any_neq(mean_degrees([1,2]), 1.5): raise TestFailure()
if any_neq(mean_degrees([1,359.1]), 0.05, rtol=1.e-12): raise TestFailure()
if any_neq(mean_degrees([0.1,359.1]), 359.6, rtol=1.e-12): raise TestFailure()

origin = (1.,)
d0 = numpy.arange(5.)*0.5
d1 = distance(5, origin=origin, resolution=0.5)
d2 = _distance_slow((5,), origin, [0.5], None)
if any_neq(d0, d1): raise TestFailure()
if any_neq(d0, d2): raise TestFailure()

origin = (1.,2.)
d0 = numpy.array([[2., numpy.sqrt(5), numpy.sqrt(8)], [0,1,2]])
d1 = distance((2,3), origin=origin, resolution=(1.,2.))
d2 = _distance_slow((2,3), origin, [1.,2.], None)
if any_neq(d0, d1): raise TestFailure()
if any_neq(d0, d2): raise TestFailure()

m = gaussian((1000,1000),fwhm=10, resolution=0.1)
if numpy.sum(m[500,:] > numpy.max(m)/2) != 100: raise TestFailure()
m = gaussian((1000,1000),fwhm=10, resolution=1)
if numpy.sum(m[500,:] > numpy.max(m)/2) != 10: raise TestFailure()
m = gaussian((1000,1000),fwhm=100, resolution=10)
if numpy.sum(m[500,:] > numpy.max(m)/2) != 10: raise TestFailure()

m = airy_disk((1000,1000),fwhm=10, resolution=0.1)
if numpy.sum(m[500,:] > numpy.max(m)/2) != 100: raise TestFailure()
m = airy_disk((1000,1000),fwhm=10, resolution=1)
if numpy.sum(m[500,:] > numpy.max(m)/2) != 10: raise TestFailure()
m = airy_disk((1000,1000),fwhm=100, resolution=10)
if numpy.sum(m[500,:] > numpy.max(m)/2) != 10: raise TestFailure()
