import numpy as np
import glob
import os
import pickle
import tamasis
from tamasis import Quantity, FitsArray, Map, Tod, create_fitsheader
from tamasis.numpyutils import any_neq, get_attributes
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

# check copies
d = np.array([1,2])
u = 'Jy'
du = {'Jy':Quantity(1,'Yj')}
h = create_fitsheader(fromdata=d)
h.update('myk', 10)
o = 'upper'
n = (0, 2)
m = np.array([True, False])
c = np.array([1., 2.])
e = np.array([3., 4.])

a = Quantity(d, u, du)
b = FitsArray(a, header=h)
c = Map(b, origin=o, error=e, coverage=c)
d = Tod(b, mask=m)
objs = [a, b, c, d]
for obj in objs:
    obj2 = obj.__class__(obj, copy=False)
    if id(obj) != id(obj2): raise TestFailure()

for obj in objs:
    obj2 = obj.__class__(obj, copy=True)
    for k in get_attributes(obj):
        if (k in ['_unit', '_derived_units', 'comm', 'shape_global', 'origin'])\
            is not (id(getattr(obj, k)) == id(getattr(obj2, k))):
            print k, id(getattr(obj, k)), id(getattr(obj2, k))

b = FitsArray(a, unit=u, derived_units=du, header=h)
c = Map(a, unit=u, derived_units=du, header=h, origin=o, error=e, coverage=c)
d = Tod(a, unit=u, derived_units=du, header=h, mask=m)

objs = [b, c, d]
for obj in objs:
    obj2 = obj.__class__(obj, copy=False)
    if id(obj) != id(obj2): raise TestFailure()

for obj in objs:
    obj2 = obj.__class__(obj, copy=True)
    for k in get_attributes(obj):
        if (k in ['_unit', '_derived_units', 'comm', 'shape_global', 'origin'])\
            is not (id(getattr(obj, k)) == id(getattr(obj2, k))):
            print k, id(getattr(obj, k)), id(getattr(obj2, k))

tod = Tod((2,))
tod = Tod([2])
tod = Tod(np.array([2]))
tod = Tod(2)
tod = Tod(np.array(2))
tod = Tod((2,))
tod = Tod([2])

a = np.ones((10,32))
tod = a.view(Tod)
assert hasattr(tod, 'mask')

tod2 = tod + 1
assert isinstance(tod2, Tod)

a = Tod([10,20])
if a.dtype.type != deftype: raise TestFailure()

a = Tod([10,20], mask=[True,False])
b = Tod(a, copy=True)
if id(a.mask) == id(b.mask): raise TestFailure()

b = Tod(a, copy=False)
if id(a) != id(b): raise TestFailure()

othertype = np.float32 if deftype is not np.float32 else np.float64
b = Tod(a, dtype=othertype, copy=False)
if id(a) == id(b): raise TestFailure()
if id(a.mask) != id(b.mask): raise TestFailure()

header = create_fitsheader((10,20))
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

a = Tod(1)
if a.shape != (): raise TestFailure()
if a.size != 1: raise TestFailure()

m = np.ndarray((10,2,10), dtype='int8')
m.flat = np.random.random(m.size)*2
a = Tod(np.random.random_sample((10,2,10)), mask=m, unit='Jy')
a.save(filename+'_tod.fits')
b = Tod(filename+'_tod.fits')
if np.any(a != b): raise TestFailure()

class MAP(Map):
    def __init__(self, data):
        self.info = 'info'

m = MAP(np.ones(3))
if get_attributes(m) != ['_derived_units', '_header', '_unit','comm', 'coverage', 'error', 'info', 'origin', 'shape_global']: raise TestFailure()

# test pickling
a = np.ones((4,3))
a[1,2] = 4
q = Quantity(a, unit='myunit', derived_units={'myunit': Quantity(2., 'Jy')})
f = FitsArray(q, header=create_fitsheader(fromdata=q, cdelt=0.5, crval=(4.,8.)))
m = Map(f, origin='upper', error=a*2, coverage=a*3)
mask = np.zeros((4,3), np.bool8)
mask[0,2] = True
t = Tod(f, mask=mask)

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
