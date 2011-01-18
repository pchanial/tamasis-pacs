import numpy
import glob
import os
import pickle
import tamasis
from tamasis import *
from tamasis.numpyutils import get_attributes
from tamasis.datatypes import validate_sliced_shape
from uuid import uuid1

tamasis.var.verbose = True
filename = 'tamasistest-'+str(uuid1())
deftype = tamasis.var.FLOAT_DTYPE.type

def test_cleanup():
    files = glob.glob(filename+'*')
    for file in files:
        os.remove(file)

class TestFailure(Exception):
    test_cleanup()


# validate scalars
if validate_sliced_shape((),None) != ((),): raise TestFailure()
if validate_sliced_shape([],None) != ((),): raise TestFailure()
if validate_sliced_shape(numpy.array(()),None) != ((),): raise TestFailure()
if validate_sliced_shape((),()) != ((),): raise TestFailure()
if validate_sliced_shape([],()) != ((),): raise TestFailure()
if validate_sliced_shape(numpy.array(()),()) != ((),): raise TestFailure()

# validate arrays of size 0
if validate_sliced_shape((0,),None) != (0,): raise TestFailure()
if validate_sliced_shape((0,),(0,)) != (0,): raise TestFailure()
if validate_sliced_shape([0],None) != (0,): raise TestFailure()
if validate_sliced_shape([0],(0,)) != (0,): raise TestFailure()

# validate arrays with slices of size 0
a = numpy.ones((1,0,1))
if validate_sliced_shape(a.shape, None) != (1,0,1): raise TestFailure()
if validate_sliced_shape(a.shape, 1) != (1,0,1): raise TestFailure()
if validate_sliced_shape(a.shape, (1,)) != (1,0,1): raise TestFailure()
a = numpy.ones((1,1,0))
if validate_sliced_shape(a.shape, None) != (1,1,0): raise TestFailure()
if validate_sliced_shape(a.shape, 0) != (1,1,0): raise TestFailure()
if validate_sliced_shape(a.shape, (0,)) != (1,1,0): raise TestFailure()
a = numpy.ones((1,1,3))

if validate_sliced_shape((10,(10,3)), None) != (10, (10,3)): raise TestFailure()
if validate_sliced_shape((10,13), None) != (10, 13): raise TestFailure()

# check copies
d = numpy.array([1,2])
u = 'Jy'
du = {'Jy':Quantity(1,'Yj')}
h = create_fitsheader(d)
h.update('myk', 10)
o = 'upper'
n = (0, 2)
m = numpy.array([True, False])
c = numpy.array([1., 2.])
e = numpy.array([3., 4.])

a = Quantity(d, u, du)
b = FitsArray(a, header=h)
c = Map(b, origin=o, error=e, coverage=c)
d = Tod(b, nsamples=n, mask=m)
objs = [a, b, c, d]
for obj in objs:
    obj2 = obj.__class__(obj, copy=False)
    if id(obj) != id(obj2): raise TestFailure

for obj in objs:
    obj2 = obj.__class__(obj, copy=True)
    for k in get_attributes(obj):
        if (k in ['_unit', 'origin', 'nsamples']) is not \
           (id(getattr(obj, k)) == id(getattr(obj2, k))):
            print k, id(getattr(obj, k)), id(getattr(obj2, k))

b = FitsArray(a, unit=u, derived_units=du, header=h)
c = Map(a, unit=u, derived_units=du, header=h, origin=o, error=e, coverage=c)
d = Tod(a, unit=u, derived_units=du, header=h, nsamples=n, mask=m)

objs = [b, c, d]
for obj in objs:
    obj2 = obj.__class__(obj, copy=False)
    if id(obj) != id(obj2): raise TestFailure

for obj in objs:
    obj2 = obj.__class__(obj, copy=True)
    for k in get_attributes(obj):
        if (k in ['_unit', 'origin', 'nsamples']) is not \
           (id(getattr(obj, k)) == id(getattr(obj2, k))):
            print k, id(getattr(obj, k)), id(getattr(obj2, k))

# check errors
try:
    junk = validate_sliced_shape((), 0)
    raise TestFailure()
except ValueError:
    pass
try:
    junk = validate_sliced_shape((0,), ())
    raise TestFailure()
except ValueError:
    pass
try:
    junk = validate_sliced_shape((1,), ())
    raise TestFailure()
except ValueError:
    pass
try:
    junk = validate_sliced_shape((1,), (-1,0,2))
    raise TestFailure()
except ValueError:
    pass
try:
    junk = validate_sliced_shape((2,), (3,4))
    raise TestFailure()
except ValueError:
    pass


tod = Tod((2,))
tod = Tod([2])
tod = Tod(numpy.array([2]))
tod = Tod(2)
tod = Tod(numpy.array(2))
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
    raise TestFailure()
except ValueError:
    pass
if tod.shape != (10,8): raise TestFailure('Tod.empty3')

tod = Tod.zeros((10,(3,5)))
if tod.shape != (10,8): raise TestFailure('Tod.zeros1')
tod = Tod.zeros((10,(3,5)), nsamples=(3,5))
if tod.shape != (10,8): raise TestFailure('Tod.zeros2')
try:
    tod = Tod.zeros((10,(3,5)), nsamples=(3,4))
    raise TestFailure()
except ValueError:
    pass
if tod.shape != (10,8): raise TestFailure('Tod.zeros3')

tod = Tod.ones((10,(3,5)))
if tod.shape != (10,8): raise TestFailure('Tod.ones1')
tod = Tod.ones((10,(3,5)), nsamples=(3,5))
if tod.shape != (10,8): raise TestFailure('Tod.ones2')
try:
    tod = Tod.ones((10,(3,5)), nsamples=(3,4))
    raise TestFailure()
except ValueError:
    pass
if tod.shape != (10,8): raise TestFailure('Tod.ones3')

tod2 = tod + 1
if tod2.shape != (10,8): raise TestFailure('Addition')
if tod2.nsamples != (3,5): raise TestFailure('Addition2')

a = Tod([10,20])
if a.dtype.type != deftype: raise TestFailure()

a = Tod([10,20], mask=[True,False])
b = Tod(a, copy=True)
if id(a.mask) == id(b.mask): raise TestFailure()

b = Tod(a, copy=False)
if id(a) != id(b): raise TestFailure()

othertype = numpy.float32 if deftype is not numpy.float32 else numpy.float64
b = Tod(a, dtype=othertype, copy=False)
if id(a) == id(b): raise TestFailure()
if id(a.mask) != id(b.mask): raise TestFailure()

header = create_fitsheader(numpy.ones((20,10)))
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

b = FitsArray(a, subok=True)
if id(a) == id(b): raise TestFailure()
if id(a.header) == id(b.header): raise TestFailure()
if id(a.unit) != id(b.unit): raise TestFailure()

b = FitsArray(a, copy=False, subok=True)
if id(a) != id(b): raise TestFailure()

a = Tod([])
if a.shape != (0,): raise TestFailure()
if a.size != 0: raise TestFailure()

a = Tod([], nsamples = (0,))
if a.shape != (0,): raise TestFailure()
if a.size != 0: raise TestFailure()

a = Tod(1)
if a.shape != (): raise TestFailure()
if a.size != 1: raise TestFailure()

a = Tod(1, nsamples = ())
if a.shape != (): raise TestFailure()
if a.size != 1: raise TestFailure()

a = Tod.ones((10,20), nsamples=(5,15))
b = a[:,:12]
if b.nsamples != (12,): raise TestFailure()
b = a[:3,:]
if b.nsamples != a.nsamples: raise TestFailure()
b = a[3]
if b.nsamples != a.nsamples: raise TestFailure()
b = a[3,:]
if b.nsamples != a.nsamples: raise TestFailure()
b = a[:,2]
if not isinstance(b, FitsArray): raise TestFailure()
b = a[4,2]
if not isinstance(b, deftype): raise TestFailure()

m = numpy.ndarray((10,2,10), dtype='int8')
m.flat = numpy.random.random(m.size)*2
a = Tod(numpy.random.random_sample((10,2,10)), mask=m, nsamples=(2,8), unit='Jy')
a.save(filename+'_tod.fits')
b = Tod(filename+'_tod.fits')
tamasis.var.verbose = False
print 'start'
tamasis.var.verbose = False
a!=b
print 'stop'
tamasis.var.verbose = False
if numpy.any(a != b): raise TestFailure()
# or numpy.any(a.mask != b.mask) or a.nsamples != b.nsamples: raise TestFailure()

class MAP(Map):
    def __init__(self, data):
        self.info = 'info'

m = MAP(numpy.ones(3))
if get_attributes(m) != ['info', 'coverage', 'error', 'origin', '_header', '_unit', '_derived_units']: raise TestFailure()

# test pickling
a = numpy.ones((4,3))
a[1,2] = 4
q = Quantity(a, unit='myunit', derived_units={'myunit': Quantity(2., 'Jy')})
f = FitsArray(q, header=create_fitsheader(q, cdelt=0.5, crval=(4.,8.)))
m = Map(f, origin='upper', error=a*2, coverage=a*3)
mask = numpy.zeros((4,3), numpy.bool8)
mask[0,2] = True
t = Tod(f, mask=mask, nsamples=(2,1))

objs = (q,f,m,t)
for v in range(pickle.HIGHEST_PROTOCOL):
    for o in objs:
        o2 = pickle.loads(pickle.dumps(o,v))
        if any_neq(o, o2): raise TestFailure()

for o in objs[1:]:
    test_cleanup()
    o.save(filename+'_obj.fits')
    o2 = type(o)(filename+'_obj.fits')
    if any_neq(o, o2):
        print 'o', repr(o)
        print 'o2', repr(o2)
        raise Exception()

test_cleanup()
