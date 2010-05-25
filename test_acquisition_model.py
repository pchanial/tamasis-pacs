import numpy
from tamasis import *

class TestFailure(Exception): pass

#--------------------
# CompressionAverage
#--------------------

tod = Tod(numpy.ndarray((2,15)), nsamples=(10,5))
tod[0,:] = [1,1,1,1,1,3,3,3,3,3,4,4,4,4,4]
tod[1,:] = [1,2,1.,0.5,0.5,5,0,0,0,0,1,2,1,2,1.5]
compression = CompressionAverage(5)
tod2 = compression.direct(tod)
if tod2.shape != (2,3): raise TestFailure('compAvg1')
if tuple(tod2.ravel()) != (1.,3.,4.,1.,1.,1.5): raise TestFailure('compAvg2')
if tod2.nsamples != (2,1): raise TestFailure('compAvg3')
tod3 = compression.transpose(tod2)
if tod3.shape != (2,15): raise TestFailure('compAvg1t')
if tuple(tod3[0,:]) != (0.2,0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8): raise TestFailure('compAvg2t')
if tod3.nsamples != (10,5): raise TestFailure('compAvg3t')

tod = Tod(numpy.ndarray((2,15)), nsamples=(11,4))
try:
    tod2 = compression.direct(tod)
except ValidationError:
    pass

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

a = numpy.array([0, 0., 1., 1., 0., 1])
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

a = numpy.array([[0, 0.], [1., 1.], [0., 1.]])
mask = Masking(a)
if mask.shapein != a.shape or mask.shapeout != a.shape:
    raise TestFailure('mask 2d')
b = numpy.array([[3., 4.], [1., 0.], [3., 2.]])
c = numpy.array([[3., 4.], [0., 0.], [3., 0.]])
if numpy.any(mask.direct(b) != c):
    raise TestFailure('mask 2d direct')
if numpy.any(mask.transpose(b) != c):
    raise TestFailure('mask 2d transpose')

a = numpy.array([[[0, 0.], [1., 1.]], [[0., 1], [1, 1]]])
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

c = mask.direct(b, reusein=True)
if id(b) != id(c):
    raise TestFailure('mask no copy')


#---------
# Padding
#---------

padding = Padding(left=1,right=20)
tod = Tod(numpy.ndarray((10,15)))
tod2 = padding.direct(tod)
if tod2.shape != (10,36):          raise TestFailure('padding shapea')
if any_neq(tod2[:,0:1], 0, 15):    raise TestFailure('padding left')
if any_neq(tod2[:,1:16], tod, 15): raise TestFailure('padding middle')
if any_neq(tod2[:,16:], 0, 15):    raise TestFailure('padding right')
tod3 = padding.transpose(tod2)
if tod3.shape != tod.shape: raise TestFailure('padding shapeb')
if any_neq(tod3,tod,15): raise TestFailure('padding end')

padding = Padding(left=1,right=(4,20))
tod = Tod.ones((10,(12,3)))
tod2 = padding.direct(tod)
if tod2.shape != (10,41):          raise TestFailure('padding shapea')
if tod2.nsamples != (17,24):       raise TestFailure('padding nsamples')
if any_neq(tod2[:,0:1], 0, 15):    raise TestFailure('padding left2')
if any_neq(tod2[:,1:13], tod[:,0:12], 15): raise TestFailure('padding middle2')
if any_neq(tod2[:,13:17], 0, 15):  raise TestFailure('padding right2')
if any_neq(tod2[:,17:18], 0, 15):  raise TestFailure('padding left2')
if any_neq(tod2[:,18:21], tod[:,12:], 15): raise TestFailure('padding middle2')
if any_neq(tod2[:,21:], 0, 15):    raise TestFailure('padding right2')
tod3 = padding.transpose(tod2)
if tod3.shape != tod.shape: raise TestFailure('padding shapeb')
if any_neq(tod3,tod,15): raise TestFailure('padding end')


#-----
# Fft
#-----

n = 100
fft = Fft(n)
a = numpy.random.random(n)+1
b = fft.transpose(fft.direct(a))
if any_neq(a, b, 14): raise TestFailure('fft1')

a = numpy.random.random((3,n))+1
b = fft.transpose(fft.direct(a))
if any_neq(a, b, 14): raise TestFailure('fft2')

tod = Tod(numpy.random.random((10,1000))+1)
fft = Fft(tod.nsamples)
tod2 = fft.transpose(fft.direct(tod))
if any_neq(tod, tod2,14): raise TestFailure('fft3')

tod = Tod(numpy.random.random((10,1000))+1, nsamples=(100,300,5,1000-100-300-5))
fft = Fft(tod.nsamples)
tod2 = fft.transpose(fft.direct(tod))
if any_neq(tod, tod2,14): raise TestFailure('fft4')

print 'OK.'
