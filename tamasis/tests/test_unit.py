import numpy
from tamasis import *

class TestFailure(Exception): pass

#test scalar, array and slices
na = numpy.array

def test_unit_add(x, y, v, u):
    c = x + y
    if not isinstance(c, Quantity): raise TestFailure(str(c))
    if na(c) != v: raise TestFailure(str(c))
    if c._unit != u: raise TestFailure(str(c))

q = Quantity(1.)
q2 = Quantity(1.)
a = numpy.array(1.)
i = 1
test_unit_add(q, q2, 2, None)
test_unit_add(q, a, 2, None)
test_unit_add(q, i, 2, None)

q = Quantity(1., 'm')
test_unit_add(q, q2, 2, {'m':1.0})
test_unit_add(q, a, 2, {'m':1.0})
test_unit_add(q, i, 2, {'m':1.0})

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
a = numpy.array(1.)
i = 1
test_unit_add(a, q, 2, None)
test_unit_add(i, q, 2, None)

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
    c = x - y
    if not isinstance(c, Quantity): raise TestFailure(str(c))
    if na(c) != v: raise TestFailure(str(c))
    if c._unit != u: raise TestFailure(str(c))

q = Quantity(1.)
q2 = Quantity(1.)
a = numpy.array(1.)
i = 1
test_unit_sub(q, q2, 0, None)
test_unit_sub(q, a, 0, None)
test_unit_sub(q, i, 0, None)

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
a = numpy.array(1.)
i = 1
test_unit_sub(a, q, 0, None)
test_unit_sub(i, q, 0, None)

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
    a.unit = 'notakloug'
except UnitError:
    pass
a=Quantity(1., 'kloug')

try:
    a.unit = 'kloug^2'
except UnitError:
    pass

# test __array_prepare__
if Quantity(10, 'm') >  Quantity(1,'km')  : raise TestFailure()
if Quantity(10, 'm') >= Quantity(1,'km')  : raise TestFailure()
if Quantity(1, 'km') <  Quantity(10,'m')  : raise TestFailure()
if Quantity(1, 'km') <= Quantity(10,'m')  : raise TestFailure()
if Quantity(1, 'km') == Quantity(1,'m')   : raise TestFailure()
if Quantity(1, 'km') != Quantity(1000,'m'): raise TestFailure()
if numpy.maximum(Quantity(10,'m'),Quantity(1,'km')) != 1000: raise TestFailure()
if numpy.minimum(Quantity(10,'m'),Quantity(1,'km')) != 10  : raise TestFailure()

# test constructor
a = numpy.array([10,20])
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


a = Tod([10,10], nsamples=(1,1), unit='m')
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
a=Quantity([1,2,3], unit='Jy')
for func in (numpy.min, numpy.max, numpy.mean, numpy.sum, numpy.std):
    b = func(a)
    if b != func(a.view(numpy.ndarray)): raise TestFailure(func)
    if b.unit != 'Jy': raise TestFailure()
    b = func(a, axis=0)
    if b != func(a.view(numpy.ndarray), axis=0): raise TestFailure(func)
    if b.unit != 'Jy': raise TestFailure()

b = numpy.var(a)
if b != numpy.var(a.view(numpy.ndarray)): raise TestFailure(numpy.var)
if b.unit != 'Jy^2': raise TestFailure()
b = numpy.var(a, axis=0)
if b != numpy.var(a.view(numpy.ndarray), axis=0): raise TestFailure()
if b.unit != 'Jy^2': raise TestFailure()

print 'OK.'
