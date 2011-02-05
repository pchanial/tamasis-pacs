import numpy as np
import tamasis
from tamasis import *

class TestFailure(Exception): pass

def testFailure():
    if np.__version__ >= '1.4': raise TestFailure()

q = Quantity(1, 'km')
if any_neq(q.SI, Quantity(1000, 'm')): raise TestFailure()

q = Quantity(1, 'm')
if any_neq(q, q.tounit('m')): raise TestFailure()
q2 = q.copy()
q.inunit('m')
if any_neq(q, q2): raise TestFailure()

q = Quantity(1, 'km').tounit('m')
if q.magnitude != 1000 or q.unit != 'm': raise TestFailure()

q = Quantity(1, 'km')
q.inunit('m')
q2 = q.copy()
if q.magnitude != 1000 or q.unit != 'm': raise TestFailure()

def test_unit_add(x, y, v, u):
    c = x + y
    if not isinstance(c, Quantity): raise TestFailure(str(c))
    if c.magnitude != v: raise TestFailure(str(c))
    if c._unit != u: raise TestFailure(str(c))

q = Quantity(1.)
q2 = Quantity(1.)
a = np.array(1.)
i = 1
test_unit_add(q, q2, 2, {})
test_unit_add(q, a, 2, {})
test_unit_add(q, i, 2, {})

q = Quantity(1., 'm')
test_unit_add(q, q2, 2, {'m':1.0})
test_unit_add(q, a, 2, {'m':1.0})
test_unit_add(q, i, 2, {'m':1.0})

if np.__version__ >= '1.4':
    q = Quantity(1.)
    q2 = Quantity(1., 'km')
    test_unit_add(q, q2, 2, {'km':1.0})

    q = Quantity(1., 'm')
    q2 = Quantity(1., 'km')
    test_unit_add(q, q2, 1001, {'m':1.0})

    q = Quantity(1., 'km')
    q2 = Quantity(1., 'm')
    test_unit_add(q, q2, 1.001, {'km':1.0})

q = Quantity(1.)
a = np.array(1.)
i = 1
test_unit_add(a, q, 2, {})
test_unit_add(i, q, 2, {})

q = Quantity(1., 'm')
q2 = Quantity(1.)
test_unit_add(q2, q, 2, {'m':1.0})
test_unit_add(a, q, 2, {'m':1.0})
test_unit_add(i, q, 2, {'m':1.0})

try:
    a=Quantity(1, 'kloug')+Quantity(3., 'babar')
except UnitError:
    pass

# SUB

def test_unit_sub(x, y, v, u):
    if np.__version__ < '1.4':
        if getattr(x, '_unit', {}) is not {} and getattr(y, '_unit', {}) is not {}:
            if x._unit != y._unit:
                print('Disabling operation on Quantities for Numpy < 1.4')
                return
    c = x - y
    if not isinstance(c, Quantity): raise TestFailure(str(c))
    if c.magnitude != v: raise TestFailure(str(c))
    if c._unit != u: raise TestFailure(str(c))

q = Quantity(1.)
q2 = Quantity(1.)
a = np.array(1.)
i = 1
test_unit_sub(q, q2, 0, {})
test_unit_sub(q, a, 0, {})
test_unit_sub(q, i, 0, {})

q = Quantity(1., 'm')
test_unit_sub(q, q2, 0, {'m':1.0})
test_unit_sub(q, a, 0, {'m':1.0})
test_unit_sub(q, i, 0, {'m':1.0})

q = Quantity(1.)
q2 = Quantity(1., 'km')
test_unit_sub(q, q2, 0, {'km':1.0})

q = Quantity(1., 'm')
q2 = Quantity(1., 'km')
test_unit_sub(q, q2, -999, {'m':1.0})

q = Quantity(1., 'km')
q2 = Quantity(1., 'm')
test_unit_sub(q, q2, 0.999, {'km':1.0})

q = Quantity(1.)
a = np.array(1.)
i = 1
test_unit_sub(a, q, 0, {})
test_unit_sub(i, q, 0, {})

q = Quantity(1., 'm')
q2 = Quantity(1.)
test_unit_sub(q2, q, 0, {'m':1.0})
test_unit_sub(a, q, 0, {'m':1.0})
test_unit_sub(i, q, 0, {'m':1.0})

try:
    a=Quantity(1, 'kloug')-Quantity(1., 'babar')
except UnitError:
    pass

# unit conversion
a=Quantity(1., 'kloug')
try:
    a.inunit('notakloug')
except UnitError:
    pass
a=Quantity(1., 'kloug')

try:
    a.inunit('kloug^2')
except UnitError:
    pass

a = Quantity(1, 'MJy/sr').tounit('uJy/arcsec^2')
if not np.allclose(a, 23.5044305391): raise TestFailure()
if a.unit != 'uJy / arcsec^2': raise TestFailure()
a = (Quantity(1, 'MJy/sr')/Quantity(1, 'uJy/arcsec^2')).SI
if not np.allclose(a, 23.5044305391): raise TestFailure()
if a.unit != '': raise TestFailure()

# test __array_prepare__
if Quantity(10, 'm') >  Quantity(1,'km')  : testFailure()
if Quantity(10, 'm') >= Quantity(1,'km')  : testFailure()
if Quantity(1, 'km') <  Quantity(10,'m')  : testFailure()
if Quantity(1, 'km') <= Quantity(10,'m')  : testFailure()
if Quantity(1, 'km') == Quantity(1,'m')   : testFailure()
if Quantity(1, 'km') != Quantity(1000,'m'): testFailure()
if np.maximum(Quantity(10,'m'),Quantity(1,'km')) != 1000: testFailure()
if np.minimum(Quantity(10,'m'),Quantity(1,'km')) != 10  : testFailure()

# test constructor
a = np.array([10,20])
b = Quantity(a)
if a.dtype == b.dtype: raise TestFailure()

a = Quantity([10.],'m')
b = Quantity(a)
if a._unit != b._unit: raise TestFailure()
if id(a) == id(b): raise TestFailure()
if id(a._unit) != id(b._unit): raise TestFailure()

b = Quantity(a, copy=False)
b[0]=1
if a[0] != b[0]: raise TestFailure()


a = Tod([10.,10.], nsamples=(1,1), unit='m')
b = Quantity(a)
if a.unit != b.unit: raise TestFailure()

b = Quantity(a, copy=False)
b[0] = 0
if a[0] != b[0]: raise TestFailure()
if id(a) == id(b): raise TestFailure()

a = Tod([10.,10.], nsamples=(1,1))
b = Quantity(a, subok=True)
if a.__class__ is not b.__class__: raise TestFailure()
if a.nsamples is not b.nsamples: raise TestFailure()

a = Tod([10.,10.], nsamples=(1,1))
b = Quantity(a, subok=True, copy=False)
if id(a) != id(b): raise TestFailure()

# test mean/sum/std/var
a=Quantity([1.3,2,3], unit='Jy')
for func in (np.min, np.max, np.mean, np.ptp, np.round, np.sum, np.std):
    b = func(a)
    if any_neq(b, func(a.view(np.ndarray))): raise TestFailure(func)
    if b.unit != 'Jy': raise TestFailure()
    if func in (np.round,):
        continue
    b = func(a, axis=0)
    if any_neq(b, func(a.view(np.ndarray), axis=0)): raise TestFailure(func)
    if b.unit != 'Jy': raise TestFailure()

b = np.var(a)
if b != np.var(a.view(np.ndarray)): raise TestFailure(np.var)
if b.unit != 'Jy^2': raise TestFailure()
b = np.var(a, axis=0)
if b != np.var(a.view(np.ndarray), axis=0): raise TestFailure()
if b.unit != 'Jy^2': raise TestFailure()

# test upcasting
if Quantity(1).dtype is not tamasis.var.FLOAT_DTYPE: raise TestFailure()
if Quantity(1, dtype='float32').dtype.type is not np.float32: raise TestFailure()
if Quantity(1.).dtype is not tamasis.var.FLOAT_DTYPE: raise TestFailure()
if Quantity(complex(1,0)).dtype.type is not np.complex128: raise TestFailure()
if Quantity(1., dtype=np.complex64).dtype.type is not np.complex64: raise TestFailure()
if Quantity(1., dtype=np.complex128).dtype.type is not np.complex128: raise TestFailure()
if Quantity(1., dtype=np.complex256).dtype.type is not np.complex256: raise TestFailure()
if Quantity(np.array(complex(1,0))).dtype is not tamasis.var.COMPLEX_DTYPE: raise TestFailure()
if Quantity(np.array(np.complex64(1.))).dtype is not tamasis.var.COMPLEX_DTYPE: raise TestFailure()
if Quantity(np.array(np.complex128(1.))).dtype is not tamasis.var.COMPLEX_DTYPE: raise TestFailure()

# test custom derived units:
derived_units = {'detector':Quantity([1, 1/10.], 'm^2')}
a = Quantity([[1,2,3,4],[10,20,30,40]], 'detector', derived_units=derived_units)
if not np.allclose(a.SI, [[1,2,3,4],[1,2,3,4]]): raise TestFailure()
a = Quantity(1, 'detector C', derived_units)
if a.SI.shape != (2,): raise TestFailure()
a = Quantity(np.ones((1,10)), 'detector', derived_units=derived_units)
if a.SI.shape != (2,10): raise TestFailure()

a = Quantity(4., 'Jy/detector', {'detector':Quantity(2,'arcsec^2')})
a.inunit(a.unit + ' / arcsec^2 * detector')
if a.magnitude != 2 or a.unit != 'Jy / arcsec^2': raise TestFailure()

a = Quantity(2., 'brou', {'brou':Quantity(2.,'bra'),'bra':Quantity(2.,'bri'),'bri':Quantity(2.,'bro'), 'bro':Quantity(2, 'bru'), 'bru':Quantity(2.,'stop')})
b = a.tounit('bra')
if b.magnitude != 4: raise TestFailure()
b = a.tounit('bri')
if b.magnitude != 8: raise TestFailure()
b = a.tounit('bro')
if b.magnitude != 16: raise TestFailure()
b = a.tounit('bru')
if b.magnitude != 32: raise TestFailure()
b = a.tounit('stop')
if b.magnitude != 64: raise TestFailure()
b = a.SI
if b.magnitude != 64 or b.unit != 'stop': raise TestFailure()

# test pixels
if any_neq(Quantity(1,'pixel/sr/pixel_reference').SI, Quantity(1, '/sr')): TestFailure()
