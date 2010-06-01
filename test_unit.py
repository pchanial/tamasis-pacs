import numpy
from tamasis import *
from unit import *

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
test_unit_add(q, q2, 2, {})
test_unit_add(q, a, 2, {})
test_unit_add(q, i, 2, {})

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
    c = x - y
    if not isinstance(c, Quantity): raise TestFailure(str(c))
    if na(c) != v: raise TestFailure(str(c))
    if c._unit != u: raise TestFailure(str(c))

q = Quantity(1.)
q2 = Quantity(1.)
a = numpy.array(1.)
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
a = numpy.array(1.)
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
    a.unit = 'notakloug'
except UnitError:
    pass
a=Quantity(1., 'kloug')

try:
    a.unit = 'kloug^2'
except UnitError:
    pass

print 'OK.'
