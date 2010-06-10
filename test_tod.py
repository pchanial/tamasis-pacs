import numpy
from tamasis import *

class TestFailure(Exception): pass

tod = Tod((2,))
tod = Tod([2])
tod = Tod(numpy.array([2]))
try:
    tod = Tod(2)
except ValueError:
    pass
try:
    tod = Tod(numpy.array(2))
except ValueError:
    pass
tod = Tod((2,), nsamples=1)
tod = Tod([2], nsamples=1)
tod = Tod(numpy.array([2]), nsamples=1)

a = numpy.ones((10,32))
tod = a.view(Tod)
if tod.nsamples != (32,): raise TestFailure()

tod = Tod.empty((10,(3,5)))
if tod.shape != (10,8): raise TestFailure('Tod.empty1')
tod = Tod.empty((10,(3,5)), nsamples=(3,5))
if tod.shape != (10,8): raise TestFailure('Tod.empty2')
try:
    tod = Tod.empty((10,(3,5)), nsamples=(3,4))
except ValidationError:
    pass
if tod.shape != (10,8): raise TestFailure('Tod.empty3')

tod = Tod.zeros((10,(3,5)))
if tod.shape != (10,8): raise TestFailure('Tod.zeros1')
tod = Tod.zeros((10,(3,5)), nsamples=(3,5))
if tod.shape != (10,8): raise TestFailure('Tod.zeros2')
try:
    tod = Tod.zeros((10,(3,5)), nsamples=(3,4))
except ValidationError:
    pass
if tod.shape != (10,8): raise TestFailure('Tod.zeros3')

tod = Tod.ones((10,(3,5)))
if tod.shape != (10,8): raise TestFailure('Tod.ones1')
tod = Tod.ones((10,(3,5)), nsamples=(3,5))
if tod.shape != (10,8): raise TestFailure('Tod.ones2')
try:
    tod = Tod.ones((10,(3,5)), nsamples=(3,4))
except ValidationError:
    pass
if tod.shape != (10,8): raise TestFailure('Tod.ones3')

tod2 = tod + 1
if tod2.shape != (10,8): raise TestFailure('Addition')
if tod2.nsamples != (3,5): raise TestFailure('Addition2')

a = Tod([10,20])
if a.dtype != numpy.float64: raise TestFailure()

a = Tod([10,20], mask=[True,False])
b = Tod(a, copy=True)
if id(a.mask) == id(b.mask): raise TestFailure()

b = Tod(a, copy=False)
if id(a) != id(b): raise TestFailure()

b = Tod(a, dtype='float32', copy=False)
if id(a) == id(b): raise TestFailure()
if id(a.mask) != id(b.mask): raise TestFailure()

header = create_fitsheader([20,10])
a = FitsArray([10,20], header=header, unit='m')
b = Tod(a)
if id(a.header) == id(b.header): raise TestFailure()
if id(a.unit) != id(b.unit): raise TestFailure()

b = Tod(a, copy=False)
if id(a.header) != id(b.header): raise TestFailure()
if id(a.unit) != id(b.unit): raise TestFailure()

a = Tod([20,10], header=header, unit='m')
b =  FitsArray(a)
if id(a.header) == id(b.header): raise TestFailure()
if id(a.unit) != id(b.unit): raise TestFailure()

b = FitsArray(a, copy=False)
if id(a.header) != id(b.header): raise TestFailure()
if id(a.unit) != id(b.unit): raise TestFailure()

b =  FitsArray(a, subok=True)
if id(a) == id(b): raise TestFailure()
if id(a.header) == id(b.header): raise TestFailure()
if id(a.unit) != id(b.unit): raise TestFailure()

b = FitsArray(a, copy=False, subok=True)
if id(a) != id(b): raise TestFailure()

print 'OK.'
