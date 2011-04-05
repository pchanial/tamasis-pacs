import numpy as np
from tamasis import *

class TestFailure(Exception):
    pass

t = Tod([1,np.nan,5,np.inf,-np.inf,8])
remove_nonfinite(t)
if any_neq(t, [1,0,5,0,0,8]): raise TestFailure()
if any_neq(t.mask, [False,True,False,True,True,False]): raise TestFailure()

t = Tod([1,np.nan,5,np.inf,-np.inf,8], mask=[False,False,False,True,False,True])
remove_nonfinite(t)
if any_neq(t, [1,0,5,0,0,8]): raise TestFailure()
if any_neq(t.mask, [False,True,False,True,True,True]): raise TestFailure()
