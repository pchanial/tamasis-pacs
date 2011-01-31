import numpy
from tamasis import *
from tamasis.utils import _distance_slow, diff as udiff, diffT, diffTdiff

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

# diff
for axis in range(2):
    a=numpy.random.random_integers(1, 10, size=(8,9)).astype(float)
    ref = -numpy.diff(a, axis=axis)
    s = ref.shape
    udiff(a, axis=axis)
    a[:s[0],:s[1]] -= ref
    if any_neq(a, 0): raise TestFailure()

for axis in range(3):
    a=numpy.random.random_integers(1, 10, size=(8,9,10)).astype(float)
    ref = -numpy.diff(a, axis=axis)
    s = ref.shape
    udiff(a, axis=axis)
    a[:s[0],:s[1],:s[2]] -= ref
    if any_neq(a, 0): raise TestFailure()

for axis in range(4):
    a=numpy.random.random_integers(1, 10, size=(8,9,10,11)).astype(float)
    ref = -numpy.diff(a, axis=axis)
    s = ref.shape
    udiff(a, axis=axis)
    a[:s[0],:s[1],:s[2],:s[3]] -= ref
    if any_neq(a, 0): raise TestFailure()

# diff.T
for axis in range(4):
    dX = DiscreteDifference(axis=axis, shapein=(3,4,5,6))
    if any_neq(dX.dense().T, dX.T.dense()): raise TestFailure()

# diff.T diff
for axis in range(4):
    dX = DiscreteDifference(axis=axis, shapein=(3,4,5,6))
    dtd = DdTdd(axis=axis, shapein=(3,4,5,6))
    if any_neq(numpy.matrix(dX.T.dense()) * numpy.matrix(dX.dense()), dtd.dense()): raise TestFailure()

# profile

def profile_slow(input, origin=None, bin=1.):
    d = distance(input.shape, origin=origin)
    d /= bin
    d = d.astype(int)
    m = numpy.max(d)
    p = numpy.ndarray(int(m+1))
    n = numpy.zeros(m+1, int)
    for i in range(m+1):
        p[i] = numpy.mean(input[d == i])
        n[i] = numpy.sum(d==i)
    return p, n

d = distance((10,20))
x, y, n = profile(d, origin=(4,5), bin=2., histogram=True)
y2, n2 = profile_slow(d, origin=(4,5), bin=2.)
if any_neq(y, y2[0:y.size]): raise TestFailure()
if any_neq(n, n2[0:n.size]): raise TestFailure()
