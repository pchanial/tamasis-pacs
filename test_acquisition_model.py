import numpy
from tamasis import *

class TestFailure(Exception): pass

#---------
# Masking
#---------

try:
    mask = Masking('lkj')
except TypeError:
    pass
else: 
    raise TestFailure('mask str')
try:
    mask = Masking(3)
except TypeError:
    pass
else: 
    raise TestFailure('mask 3')
try:
    mask = Masking(object)
except TypeError:
    pass
else: 
    raise TestFailure('mask object')
try:
    mask = Masking(numpy.array(3))
except ValueError:
    pass
else: 
    raise TestFailure('mask ndarray rank=0')

mask = Masking(None)
if mask.shapein is not None or mask.shapeout is not None or mask.mask is not None:
    raise TestFailure('Masking(None)')

try:
    junk = mask.direct('str')
except TypeError:
    pass
else: 
    raise TestFailure('mask None wrong type str')
try:
    junk = mask.direct(32)
except TypeError:
    pass
else: 
    raise TestFailure('mask None wrong type int')
try:
    junk = mask.direct(object)
except TypeError:
    pass
else: 
    raise TestFailure('mask None wrong type object')

a = numpy.array([0, 0., 1., 3., 0., 2])
mask = Masking(a)
if mask.shapein != a.shape or mask.shapeout != a.shape:
    raise TestFailure('mask 1d')
b = numpy.array([3., 4., 1., 0., 3., 2.])
c = numpy.array([3., 4., 0., 0., 3., 0.])
if id(b) == id(c):
    raise TestFailure('mask 1d copy')
if numpy.any(mask.direct(b) != c):
    raise TestFailure('mask 1d direct')
if numpy.any(mask.transpose(b) != c):
    raise TestFailure('mask 1d transpose')

a = numpy.array([[0, 0.], [1., 3.], [0., 2]])
mask = Masking(a)
if mask.shapein != a.shape or mask.shapeout != a.shape:
    raise TestFailure('mask 2d')
b = numpy.array([[3., 4.], [1., 0.], [3., 2.]])
c = numpy.array([[3., 4.], [0., 0.], [3., 0.]])
if numpy.any(mask.direct(b) != c):
    raise TestFailure('mask 2d direct')
if numpy.any(mask.transpose(b) != c):
    raise TestFailure('mask 2d transpose')

a = numpy.array([[[0, 0.], [1., 3.]], [[0., 2], [1, 1]]])
mask = Masking(a)
if mask.shapein != a.shape or mask.shapeout != a.shape:
    raise TestFailure('mask 2d')
b = numpy.array([[[3, 4.], [1., 0.]], [[3., 2], [-1, 9]]])
c = numpy.array([[[3, 4.], [0., 0.]], [[3., 0], [0, 0]]])
if id(b) == id(c):
    raise TestFailure('mask 3d copy')
if numpy.any(mask.direct(b) != c):
    raise TestFailure('mask 3d direct')
if numpy.any(mask.transpose(b) != c):
    raise TestFailure('mask 3d transpose')

c = mask.direct(b, copyout=False)
if id(b) != id(c):
    raise TestFailure('mask no copy')

print 'OK.'
