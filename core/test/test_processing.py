import numpy
from tamasis import *

class TestFailure(Exception):
    pass

t = Tod([1,2,5,numpy.nan,numpy.nan,8])
remove_nan(t)
if any_neq(t, [1,2,5,0,0,8]): raise TestFailure()
if any_neq(t.mask, [False,False,False,True,True,False]): raise TestFailure()

t = Tod([1,2,5,numpy.nan,numpy.nan,8], mask=[False,False,False,True,False,True])
remove_nan(t)
if any_neq(t, [1,2,5,0,0,8]): raise TestFailure()
if any_neq(t.mask, [False,False,False,True,True,True]): raise TestFailure()

print 'OK.'
